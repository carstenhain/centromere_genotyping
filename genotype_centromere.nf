#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samplesheet = "samplesheet.csv"
params.outdir = null
params.kmer_fasta = "/g/korbel/hain/centromere_genotyping_tool/centromere_genotyping/data/all_tagging_kmers.fasta"
params.reference_fasta_for_cram = "/scratch/hain/centromere/tool_data/GRCh38_full_analysis_set_plus_decoy_hla.fa"

// Help message
def helpMessage() {
    log.info"""
    Usage:
        nextflow run genotype_centromere.nf --samplesheet <path> --outdir <path>

    Required arguments:
        --samplesheet PATH    Path to the CSV samplesheet file
        --outdir PATH         Directory where result files are stored

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
    # rm ${chunk_file.baseName}.jf
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
    ### rewrite fasta to a list of kmers
    awk 'NR % 2 == 0' "${params.kmer_fasta}" > kmer_list.txt

    ### create local uncompressed copy of the dbg file in task work directory
    gzip -dc "${file_path}" > dbg.uncompressed.fa

    ### subset dbg to those contigs with one matching kmers and write to new file
    pv dbg.uncompressed.fa | grep -B 1 -f kmer_list.txt > "${name}".kmers.zgrep.txt
    ### pv -f dbg.uncompressed.fa 2> pv.log | grep -B 1 -f kmer_list.txt > ${name}.kmers.zgrep.txt
    ### pv -f dbg.uncompressed.fa | stdbuf -oL grep -B 1 -f kmer_list.txt > ${name}.kmers.zgrep.txt
    ### cat dbg.uncompressed.fa | stdbuf -i0 -o0 -e0 grep -B 1 -f kmer_list.txt > ${name}.kmers.zgrep.txt
    

    ### call python file to convert the zgrep output into a kmer table
    python3 ${projectDir}/scripts/dbg_to_kmer_table.py \\
        --name "${name}" \\
        --zgrep_file "${name}.kmers.zgrep.txt"  \\
        --kmer_list kmer_list.txt \\
        --output "${name}_kmer_table.tsv"

    ### remove kmer list and uncompressed dbg file to save space
    # rm kmer_list.txt
    # dbg.uncompressed.fa
    """
}

process MERGE_KMER_COUNTS {
    input:
    path(samplesheet)
    path(jellyfish_output_list)

    output:
    path("reads_kmer_counts.tsv")

    script:
    """
    python3 ${projectDir}/scripts/merge_kmer_counts.py \\
        --samplesheet ${samplesheet} \\
        --jellyfish_output_list ${jellyfish_output_list} > reads_kmer_counts.tsv
    """
}

process FINAL_KMER_MERGE {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(kmer_subtable_list)

    output:
    path("final_kmer_merged.tsv")

    script:
    """
    python3 ${projectDir}/scripts/final_kmer_merge.py \
        --kmer_subtable_list ${kmer_subtable_list} > final_kmer_merged.tsv
    """
}

process GENOTYPE {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(samplesheet)
    path(final_kmer_merged)

    output:
    path("genotype_input_manifest.txt")

    script:
    """
    {
        echo "samplesheet=${samplesheet}"
        echo "final_kmer_merged=${final_kmer_merged}"
    } > genotype_input_manifest.txt
    """
}

workflow {

    if (!params.outdir) {
        error "Missing required parameter: --outdir"
    }
    if (!params.samplesheet) {
        error "Missing required parameter: --samplesheet"
    }

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
    
    // run jellyfish count and query and save the output files into a list
    def jellyfish_results_ch = jellyfish(chunks_ch)
    def jellyfish_output_list_ch = jellyfish_results_ch
        .map { name, output_file -> output_file.toString() }
        .collectFile(name: 'jellyfish_outputs.txt', newLine: true)

    // merge the jellyfish output files into a single table of kmer counts for the reads
    def reads_kmer_counts_ch = MERGE_KMER_COUNTS(
        channel.fromPath(params.samplesheet),
        jellyfish_output_list_ch
    )

    // Process dbg files into kmer tables
    def dbg_kmer_tables_ch = DBG_TO_KMER_TABLE(branched_ch.dbg.map { name, file_path, type -> tuple(name, file_path) })

    // write paths of kmer tables from reads and dbg files into a single file for merging 
    def kmer_table_list_ch = reads_kmer_counts_ch
        .map { reads_table -> reads_table.toString() }
        .mix(dbg_kmer_tables_ch.map { name, dbg_table -> dbg_table.toString() })
        .collectFile(name: 'kmer_table_list.txt', newLine: true)

    // merge kmer tables from different data types together into a single table for genotyping
    def final_kmer_merged_ch = FINAL_KMER_MERGE(kmer_table_list_ch)


}