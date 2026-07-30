#!/usr/bin/env nextflow

/*
 * Removal of adapters and low quality bases with Trimmomatic
 */
process TRIMMOMATIC {

    tag "${sample_id}"
    publishDir "results/trimmomatic/${sample_id}", mode: 'copy'

    conda "bioconda::trimmomatic=0.40"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("1_paired.fastq"),
          path("2_paired.fastq"),
          path("1_unpaired.fastq"),
          path("2_unpaired.fastq")

    script:
    """
    trimmomatic PE \\
        -threads ${task.cpus} \\
        ${reads[0]} ${reads[1]} \\
        1_paired.fastq 1_unpaired.fastq \\
        2_paired.fastq 2_unpaired.fastq \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
