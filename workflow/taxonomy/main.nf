#!/usr/bin/env nextflow
// Running the taxonomy workflow

nextflow.enable.dsl = 2

// Parameters
params.input_raw = "../../data/fastq_samples/*_{1,2}.fastq"

params.host_link_db = [
    "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip",
    "https://genome-idx.s3.amazonaws.com/bt/ARS-UCD1.2.zip"
]
params.kraken_db = "${projectDir}/data/kraken"

// Modules
include { FASTQC }          from './modules/1_fastqc.nf'
include { TRIMMOMATIC }     from './modules/2_trimming.nf'
include { HOST_DB }         from './modules/host.nf'
include { BOWTIE2_DECONTAM } from './modules/3_decontamination.nf'
include { FASTQC2 }          from './modules/4_fastqc.nf'
include { KRAKEN } from './modules/kraken.nf'

// Workflow
workflow {

    // 1. RAW READS
    raw_reads = channel
        .fromFilePairs(params.input_raw)
        .map { sample_id, reads -> tuple(sample_id, reads) }

    // 2. QC + TRIMMING
    fastqc_results = FASTQC(raw_reads)
    trimmed = TRIMMOMATIC(raw_reads)

    // HOST INDEX (OUTPUT = path to folders: data/host/species1/, host/species2/, ...)
    host_db = HOST_DB(params.host_link_db)
	.flatten()

    // 3. COMBINAISON DES ÉCHANTILLONS AVEC CHAQUE INDEX (CARTESIEN)
    paired_reads = trimmed.map { sample_id, r1_paired, r2_paired, r1_unpaired, r2_unpaired ->
    	tuple(sample_id, r1_paired, r2_paired)
    }

    decontam_input = paired_reads
        .combine(host_db)
    
    decontam = BOWTIE2_DECONTAM(decontam_input)
    clean_reads = decontam.map { sample_id, clean1, clean2, sam, log -> tuple(sample_id, [clean1, clean2]) }

    // 4. QC
    fastqc_results2 = FASTQC2(clean_reads)

    // 5. TAXONOMY
    kraken = KRAKEN(clean_reads, params.kraken_db)
    
    // OUTPUTS
    publish_fq     = fastqc_results
    publish_trim   = trimmed
    publish_deco   = decontam
    publish_fq2    = fastqc_results2
    publish_krak   = kraken
}
