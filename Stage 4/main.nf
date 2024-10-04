#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
    N A N O P O R E - NF    P I P E L I N E
    ========================================
    fastq paired-end reads        : ${params.reads}
    reference genome file : ${params.reference}
    outdir                   : ${params.outdir}
    """

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromFilePairs( params.reads )
    .view()

reference_ch = Channel
    .fromPath( params.reference )
    .view()

    
process fastqc {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    
    input:
    tuple val(sample), path(read) 
    
    script:
    """
    fastqc ${read}
    """
}

process fastp {
    input:
    tuple val(sample), path(reads)

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o out.R1.fq.gz -O out.R2.fq.gz
    """
}

/**
 * Workflow Execution
 */
workflow {
    fastqc(reads_ch)                     // Run FastQC on reads
    fastp(reads_ch)      // Run Fastp for trimming
}


workflow.onComplete {
    log.info (workflow.success ? "\n Your job is done :) ! You can check out your files here -> $params.outdir\n" : "Check your documents, the pipeline didn't work :(" )
}  





#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
    N A N O P O R E - NF    P I P E L I N E
    ========================================
    fastq paired-end reads        : ${params.reads}
    reference genome file         : ${params.reference}
    outdir                        : ${params.outdir}
    """

/**
 * Input channels
 */
reads_ch = Channel
    .fromFilePairs( params.reads )
    .view()

reference_ch = Channel
    .fromPath( params.reference )
    .view()

/**
 * FastQC for quality control
 */
process fastqc {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    
    input:
    tuple val(sample), path(reads)
    
    script:
    """
    fastqc ${reads[0]} ${reads[1]}
    """
}

/**
 * Fastp for quality control and trimming
 */
process fastp {
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("out.R1.fq.gz"), path("out.R2.fq.gz")

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o out.R1.fq.gz -O out.R2.fq.gz
    """
}

/**
 * BWA index process (index reference genome)
 */
process bwa_index {
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'
    
    input:
    path(reference)

    output:
    path("${reference}.bwt")
    
    script:
    """
    bwa index ${reference}
    """
}

/**
 * BWA mem for read alignment
 */
process bwa_mem {
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'
    
    input:
    tuple val(sample), path(reference)
    
    output:
    path "out.bam"

    script:
    """
    bwa mem ${reference} ${sample[1]} ${sample[2]} | samtools view -Sb - > out.bam
    """
}

/**
 * FreeBayes for variant calling
 */
process freebayes {
    container 'quay.io/biocontainers/freebayes:1.3.6--hd4865b3_0'
    
    input:
    path(reference)
    path "out.bam"
    
    output:
    path "var.vcf"

    script:
    """
    freebayes -f ${reference} ${out.bam} > var.vcf
    """
}

/**
 * Workflow Execution
 */
workflow {
    // Step 1: Run FastQC
    fastqc(reads_ch)

    // Step 2: Run Fastp for trimming
    trimmed_reads = fastp(reads_ch)

    // Step 3: Index the reference genome
    indexed_reference = bwa_index(reference_ch)

    // Step 4: Run BWA mem for alignment (correctly unpacking trimmed reads)
    aligned_bam = bwa_mem(trimmed_reads.collect())

    // Step 5: Run FreeBayes for variant calling
    freebayes(indexed_reference, aligned_bam)
}

workflow.onComplete {
    log.info (workflow.success ? "\n Your job is done :) ! You can check out your files here -> $params.outdir\n" : "Check your documents, the pipeline didn't work :(" )
}
