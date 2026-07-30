#!/usr/bin/env nextflow

/*
 * Quality checking with kraken2
 */

process KRAKEN {

    tag "${sample_id}"
    publishDir "results/kraken/${sample_id}", mode: 'copy'

    conda "bioconda::kraken2=2.1.3 bioconda::bracken=2.8 bioconda::krona=2.8.1"
    
    input:
        tuple val(sample_id), path(reads)
        val kraken_db

    output:
        tuple val(sample_id),
              path("${sample_id}.kraken"),
              path("${sample_id}.kraken.report"),
              path("${sample_id}.bracken", optional: true)

    script:
    """
    kraken2 \
        --db ${kraken_db} \
        --paired ${reads[0]} ${reads[1]} \
        --threads ${task.cpus} \
        --confidence 0.1  \
        --gzip-compressed \
        --output ${sample_id}.kraken \
        --report ${sample_id}.kraken.report

    #ktImportTaxonomy \
     #   ${sample_id}.kraken.report \
      #  -o ${sample_id}.kraken.html \

    bracken \
        -d ${kraken_db} \
        -i ${sample_id}.kraken.report \
        -o ${sample_id}.bracken \
        -r 150 \
        -t ${task.cpus} || true
    """
}
