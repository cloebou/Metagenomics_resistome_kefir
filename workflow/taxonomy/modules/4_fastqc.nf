#!/usr/bin/env nextflow

/*
 * Quality checking with fastQC
 */

process FASTQC2 {
    tag "${sample_id}" // Sert à affichier le sample id dans les logs
    publishDir "results/fastqc2/${sample_id}", mode: 'copy'
    conda "bioconda::fastqc=0.12.1"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*.zip')

    script:
    """
    fastqc \\
        ${reads.join(" ")} \\
        -t ${task.cpus} \\
        -o . || true

    """
}
