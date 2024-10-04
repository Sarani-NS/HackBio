#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq1 = "$baseDir/data/ERR8774458_1.fastq.gz"
params.fastq2 = "$baseDir/data/ERR8774458_2.fastq.gz"
params.reference = "$baseDir/data/Reference.fasta"

workflow {
    // Create channels for input files
    fastq_files = Channel.from(params.fastq1, params.fastq2)
    reference = file(params.reference)

    // Create a directory for FastQC results
    file("fastqc_results").mkdirs()

    // FastQC for both files
    fastqc_results = fastq_files.map { fastqc(it) }

    // Fastp for quality control and trimming
    trimmed_files = fastp(params.fastq1, params.fastq2)

    // BWA Index
    bwa_index(reference)

    // BWA MEM to align reads
    bam_file = bwa_mem(reference, trimmed_files)

    // Sort the BAM file
    sorted_bam_file = samtools_sort(bam_file)

    // FreeBayes variant calling
    freebayes(reference, sorted_bam_file)
}

process fastqc {
    input:
    path fastq_file

    output:
    path "fastqc_results/${fastq_file.baseName}_fastqc.html"
    path "fastqc_results/${fastq_file.baseName}_fastqc.zip"

    script:
    """
    fastqc $fastq_file -o fastqc_results/
    """
}

process fastp {
    input:
    path r1
    path r2

    output:
    tuple path('out.R1.fq.gz'), path('out.R2.fq.gz')

    script:
    """
    fastp -i $r1 -I $r2 -o out.R1.fq.gz -O out.R2.fq.gz
    """
}

process bwa_index {
    input:
    path reference

    script:
    """
    bwa index $reference
    """
}

process bwa_mem {
    input:
    path reference
    tuple path(r1), path(r2)

    output:
    path 'out.bam'

    script:
    """
    bwa mem $reference $r1 $r2 | samtools view -b -o out.bam -
    """
}

process samtools_sort {
    input:
    path bam_file

    output:
    path 'out_sorted.bam'

    script:
    """
    samtools sort $bam_file -o out_sorted.bam
    """
}

process freebayes {
    input:
    path reference
    path sorted_bam_file

    output:
    path 'var.vcf'

    script:
    """
    freebayes -f $reference $sorted_bam_file > var.vcf
    """
}
