#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samplesheet = "samplesheet.tsv"

process RUN_PROCESS {
    tag "${name}"

    input:
    tuple val(name), val(file_path)

    script:
    """
    python3 ${projectDir}/process.py --name "${name}" --file_path "${file_path}"
    """
}

workflow {
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.NAME, row.FILE_PATH) }
        .set { samples_ch }

    RUN_PROCESS(samples_ch)
}
