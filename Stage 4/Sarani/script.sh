#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.outdir = './results'

process download {
    output:
    path 'data/*'

    script:
    """
    mkdir -p data
    wget -P data/ https://raw.githubusercontent.com/josoga2/yt-dataset/main/dataset/raw_reads/reference.fasta
    wget -P data/ https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R1.fastq.gz
    wget -P data/ https://github.com/josoga2/yt-dataset/raw/main/dataset/raw_reads/ACBarrie_R2.fastq.gz
    # Add more downloads for Alsen, Baxter, Chara, Drysdale
    """
}

process quality_control {
    input:
    path fastq_files

    output:
    path "qc/*"

    script:
    """
    mkdir -p qc
    fastqc -o qc ${fastq_files}
    """
}

process trimming {
    input:
    path fastq_file

    output:
    path "trimmed/*"

    script:
    """
    mkdir -p trimmed
    fastp -i $fastq_file -o trimmed/trimmed_${fastq_file.baseName}
    """
}

process mapping {
    input:
    path trimmed_fastq

    output:
    path "mapped/*"

    script:
    """
    mkdir -p mapped
    bwa mem reference.fasta $trimmed_fastq > mapped/${trimmed_fastq.baseName}.sam
    """
}

process variant_calling {
    input:
    path mapped_sam

    output:
    path "variants/*"

    script:
    """
    mkdir -p variants
    samtools view -bS $mapped_sam > ${mapped_sam.baseName}.bam
    samtools sort ${mapped_sam.baseName}.bam -o sorted_${mapped_sam.baseName}.bam
    bcftools mpileup -f reference.fasta sorted_${mapped_sam.baseName}.bam | bcftools call -mv -o variants/${mapped_sam.baseName}.vcf
    """
}

workflow {
    download()
        .view { fastq_files -> fastq_files }
        .set { fastq_files }

    fastq_files.each { fq ->
        quality_control(fq)
            .set { qc_results }
        trimming(qc_results)
            .set { trimmed_fastq }
        mapping(trimmed_fastq)
            .set { mapped_sam }
        variant_calling(mapped_sam)
    }
}
