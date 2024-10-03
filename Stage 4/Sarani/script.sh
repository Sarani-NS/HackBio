#!/usr/bin/env nextflow

Current folder (in which I run all the commands): /home/sarani_admin/ngs_pipeline

fastqc data/ACBarrie_R1.fastq.gz
fastqc data/ACBarrie_R2.fastq.gz
fastp -i data/ERR8774458_1.fastq.gz -I data/ERR8774458_2.fastq.gz -o out.R1.fq.gz -O
 out.R2.fq.gz
bwa index data/Reference.fasta
bwa mem data/Reference.fasta data/ERR8774458_1.fastq.gz data/ERR8774458_2.fastq.gz | samtools view -o out.bam
freebayes -f data/Reference.fasta out_sorted.bam > var.vcf
