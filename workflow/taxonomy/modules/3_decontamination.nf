#!/usr/bin/env nextflow

/*
 * Removal of host contamination with Bowtie2
 */

process BOWTIE2_DECONTAM {

   tag "${sample_id}"
   publishDir "results/decontamination/${sample_id}", mode: 'copy'


   conda "bioconda::bowtie2=2.5.1"
   
   input:
        tuple val(sample_id), path(read1), path(read2), path(host_dir)

   output:
        tuple val(sample_id),
              path("cleaned_1.fastq.gz", optional: true),
              path("cleaned_2.fastq.gz", optional: true),
              path("mapped_and_unmapped.sam"),
	      path("bowtie2.log")

   script:
   """
   bowtie2 \
        -p ${task.cpus} \
        --very-sensitive-local \
        -x ${host_dir}/${host_dir.getName()} \
        -1 ${read1} \
        -2 ${read2} \
        --un-conc-gz cleaned \
        -S mapped_and_unmapped.sam \
	2> bowtie2.log
   mv cleaned.1 cleaned_1.fastq.gz
   mv cleaned.2 cleaned_2.fastq.gz
   """
}
