#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samplesheet = "samplesheet.tsv"

process convertToFasta {
    tag "${name}"

    input:
    tuple val(name), path(file_path)

    output:
    tuple val(name), path("${name}.fasta")

    script:
    """
    samtools fasta "${file_path}" > "${name}.fasta"
    """
}

process splitReadsIntoFastaChunks {
    tag "${name}"

    input:
    tuple val(name), path(fasta_file)

    output:
    tuple val(name), path("${name}_chunk_*.fasta")

    script:
    """
    python3 ${projectDir}/scripts/split_fasta.py --input "${fasta_file}" --prefix "${name}_chunk_"
    """
}

process jellyfishCount {
    tag "${name}"

    input:
    tuple val(name), path(fasta_chunk)

    output:
    tuple val(name), path("${fasta_chunk.baseName}.jf")

    script:
    """
    jellyfish count -m 31 -s 100M -t 4 -C -o "${fasta_chunk.baseName}.jf" "${fasta_chunk}"
    """
}

process jellyfishQuery {
    tag "${name}"

    input:
    tuple val(name), path(jf_file)

    output:
    tuple val(name), path("${jf_file.baseName}_counts.tsv")

    script:
    """
    jellyfish query -s ${projectDir}/data/kmers.tsv "${jf_file}" > "${jf_file.baseName}_counts.tsv"
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
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.NAME, row.FILE_PATH, row.TYPE) }
        .branch {
            reads: it[2] in ['bam', 'cram', 'fastq']
            dbg:   it[2] == 'dbg'
        }
        .set { branched_ch }

    reads_ch = branched_ch.reads.map { name, file_path, type -> tuple(name, file(file_path)) }
    dbg_ch   = branched_ch.dbg.map   { name, file_path, type -> tuple(name, file(file_path)) }

    convertToFasta(reads_ch) \
        | splitReadsIntoFastaChunks \
        | transpose \
        | jellyfishCount \
        | jellyfishQuery

    DBG_TO_KMER_TABLE(dbg_ch)
}
