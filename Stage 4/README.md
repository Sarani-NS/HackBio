# Simple NGS Analysis Pipeline

This pipeline is designed to process Next-Generation Sequencing (NGS) data, perform quality control, trim the reads, align them to a reference genome, and call variants using common bioinformatics tools.

![image](https://github.com/user-attachments/assets/c515bed3-0fa7-4e90-81ca-3096a05a0557)

---

## Pipeline Overview

**Data Download (wget)**: Downloads the paired-end reads & the reference genome files in the precised directory.

**Quality Control (FastQC)**: Assess the quality of raw reads.

**Trimming (FastP)**: Trim low-quality bases and adapters from the reads.

**Genome Mapping (BWA)**: Align the trimmed reads to the reference genome.

**Variant Calling (FreeBayes)**: Identify variants (SNPs and indels) from the aligned reads.

### Requirements

For this pipeline, you will need the following tools:

- fastqc
- fastp
- bwa
- samtools
- bcftools
- freebayes

Install them with setup.sh:

```bash
bash setup.sh
```

### Running the Pipeline

The pipeline can be executed by running the script.sh file. The script downloads the data then processes them before performing the full NGS analysis. 

```bash
bash script.sh
```

After running the pipeline, you will get the following files of interest:

**The downloaded data**

*Forward reads*: ERR8774458_1.fastq.gz

*Reverse reads*: ERR8774458_2.fastq.gz

*Reference genome*: Reference.fasta

**The output**

*out_R1.fq.gz, out_R2.fq.gz*: Trimmed paired-end reads

*out.bam*: Aligned reads in BAM format

*var.vcf*: Variant call format (VCF) file with identified variants

---

## More dataset to try the pipeline on 

### Reference
https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta
### ACBarrie
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz
### Alsen
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R1.fastq.gz
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Alsen_R2.fastq.gz
### Baxter
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R1.fastq.gz
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Baxter_R2.fastq.gz
### Chara
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R1.fastq.gz
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Chara_R2.fastq.gz
### Drysdale
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R1.fastq.gz
https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/Drysdale_R2.fastq.gz
