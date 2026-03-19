#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samplesheet = "samplesheet.csv"
params.kmer_fasta = "/g/korbel/hain/centromere_genotyping_tool/centromere_genotyping/data/all_tagging_kmers.fasta"
params.reference_fasta_for_cram = "/scratch/hain/centromere/tool_data/GRCh38_full_analysis_set_plus_decoy_hla.fa"

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run genotype_centromere.nf --samplesheet <path>

    Required arguments:
      --samplesheet PATH    Path to the CSV samplesheet file

    Optional arguments:
      --help                Show this help message and exit

    The samplesheet should be a CSV file with columns:
      - NAME: Sample identifier (string)
      - FILE_PATH: Path to the input file (string)
      - TYPE: Type of the input file (string, e.g., 'bam', 'cram', 'fasta', 'fastq', 'dbg', gzipped file should be named as their filetype, e.g. 'fastq' for 'fastq.gz')
    """.stripIndent()
}

process convertToFasta {
    /*
    Convert input files to FASTA format if they are not already in FASTA. Supported input types are BAM, CRAM, FASTQ, and FASTA.
    */
    tag "${name}"

    input:
    tuple val(name), path(input_file), val(file_type)

    output:
    tuple val(name), path("${name}.fasta")

    script:
    if (file_type == 'bam') {
        """
        samtools fasta "${input_file}" > "${name}.fasta"
        """
    }
    else if (file_type == 'cram') {
        """
        samtools view -@ 4 -h -T ${params.reference_fasta_for_cram} "${input_file}" | samtools fasta - > "${name}.fasta"
        """
    }
    else if (file_type == 'fastq') {
        """
        seqtk seq -A "${input_file}" > "${name}.fasta"
        """
    }
    else {
        error "Unsupported file format in convertToFasta: ${file_type}"
    }
}

process splitReadsIntoFastaChunks {
    /*
    Split FASTA files into smaller chunks for downstream processing. 
    The number of splits is determined based on the total number of reads, with a target batch size of 1 million reads per chunk. 
    The process uses seqkit for splitting and ensures that the number of splits is between 2 and 999 to avoid creating too few or too many chunks.
    */
    tag "${name}"

    input:
    tuple val(name), path(fasta_file)

    output:
    tuple val(name), path("*.chunk.fasta"), emit: chunks_with_meta

    script:
    """
    # Count number of sequences
    READ_COUNT=\$(seqtk size ${fasta_file} | awk '{print \$1}')
    echo "Total reads: \$READ_COUNT" > ${name}.log
    
    BATCH_SIZE=1000000
    NUM_SPLITS=\$((READ_COUNT/BATCH_SIZE))
    
    # Ensure NUM_SPLITS is at least 2 and at most 999
    if [ \$NUM_SPLITS -lt 2 ]; then
        NUM_SPLITS=2
    elif [ \$NUM_SPLITS -gt 999 ]; then
        NUM_SPLITS=999
    fi
    echo "Number of splits: \$NUM_SPLITS" >> ${name}.log

    # Split sequence file using seqkit instead of fastp
    seqkit split2 -p \$NUM_SPLITS ${fasta_file} -O ./
    
    # Rename the split files to match expected pattern
    for file in ${fasta_file.baseName}.part_*.fasta; do
        if [ -f "\$file" ]; then
            basename_file=\$(basename "\$file")
            mv "\$file" "\${basename_file%.fasta}.chunk.fasta"
        fi
    done
    
    # Log the actual files created
    echo "Files created:" >> ${name}.log
    ls -la *.chunk.fasta >> ${name}.log
    """
}

process jellyfish {
    /*
    Count k-mers in the chunked FASTA files using Jellyfish. 
    The process creates a Jellyfish database for each chunk and then queries the database for the specified k-mers. 
    The results are saved to a text file for each chunk. 
    The Jellyfish database files are removed after querying to save disk space.
    */
    tag "${name}_${chunk_file.baseName}"

    input:
    tuple val(name), path(chunk_file)

    output:
    tuple val(name), path("${name}_${chunk_file.baseName}_query_results.txt")

    script:
    """
    jellyfish count \\
        -m 61 \\
        -s 2G \\
        -t ${task.cpus} \\
        -C \\
        ${chunk_file} \\
        -o ${chunk_file.baseName}.jf

    jellyfish query \\
        ${chunk_file.baseName}.jf \\
        -s ${params.kmer_fasta} \\
        -o ${name}_${chunk_file.baseName}_query_results.txt
    # remove jf database to save space
    rm ${chunk_file.baseName}.jf
    """
}

process DBG_TO_KMER_TABLE {
    tag "${name}"

    input:
    tuple val(name), path(file_path)

    output:
    tuple val(name), path("${name}_kmer_table.tsv")

    script:
    """
    python3 ${projectDir}/scripts/dbg_to_kmer_table.py --name "${name}" --file_path "${file_path}" --output "${name}_kmer_table.tsv"
    """
}

workflow {

    def input_ch = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ';')
        .map { row -> tuple(row.NAME, file(row.FILE_PATH), row.TYPE.toLowerCase()) }

    def branched_ch = input_ch
        .branch { name, file_path, type ->
            reads: type in ['bam', 'cram', 'fastq']
            fasta: type == 'fasta'
            dbg:   type == 'dbg'
        }

    // Process mapped reads and fastq into fasta and combine with direct fasta files
    def fasta_ch = convertToFasta(branched_ch.reads).mix(branched_ch.fasta.map { name, file_path, type -> tuple(name, file_path) })

    splitReadsIntoFastaChunks(fasta_ch)
        .transpose() // Flatten the list of chunk files into separate emissions
        .set { chunks_ch } // Count kmers in each chunk 
    
    jellyfish(chunks_ch) 

    // Process dbg files
    DBG_TO_KMER_TABLE(branched_ch.dbg.map { name, file_path, type -> tuple(name, file_path) })
}