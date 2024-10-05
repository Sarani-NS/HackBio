#!/bin/bash

# Download files and save them directly to the desired directory
wget -O /home/sarani_admin/ngs_pipeline/data/ERR8774458_1.fastq.gz "https://zenodo.org/record/10426436/files/ERR8774458_1.fastq.gz"
wget -O /home/sarani_admin/ngs_pipeline/data/ERR8774458_2.fastq.gz "https://zenodo.org/record/10426436/files/ERR8774458_2.fastq.gz"
wget -O /home/sarani_admin/ngs_pipeline/data/Reference.fasta "https://zenodo.org/records/10886725/files/Reference.fasta" 

# Define parameters
READ1="/home/sarani_admin/ngs_pipeline/data/ERR8774458_1.fastq.gz"
READ2="/home/sarani_admin/ngs_pipeline/data/ERR8774458_2.fastq.gz"
REFERENCE="/home/sarani_admin/ngs_pipeline/data/Reference.fasta"
OUTDIR="/home/sarani_admin/ngs_pipeline/data"

# Step 1: FastQC for quality control
echo "Running FastQC..."
fastqc $READ1 -o $OUTDIR
fastqc $READ2 -o $OUTDIR

# Step 2: Fastp for trimming and quality control
echo "Running Fastp..."
fastp -i $READ1 -I $READ2 -o $OUTDIR/out.R1.fq.gz -O $OUTDIR/out.R2.fq.gz

# Step 3: BWA index for reference genome
echo "Indexing reference genome with BWA..."
bwa index $REFERENCE

# Step 4: BWA mem for read alignment and Samtools for BAM conversion
echo "Aligning reads with BWA mem and converting to BAM..."
bwa mem $REFERENCE $OUTDIR/out.R1.fq.gz $OUTDIR/out.R2.fq.gz | samtools view -Sb - > $OUTDIR/out.bam

# Step 5: Sort BAM file using Samtools
echo "Sorting BAM file..."
samtools sort $OUTDIR/out.bam -o $OUTDIR/out_sorted.bam

# Step 6: FreeBayes for variant calling
echo "Calling variants with FreeBayes..."
freebayes -f $REFERENCE $OUTDIR/out_sorted.bam > $OUTDIR/var.vcf

echo "Pipeline completed successfully. Output files are in $OUTDIR."
