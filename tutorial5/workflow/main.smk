#!/usr/bin/env python3

# -*- python -*-
from os.path import join

SAMPLES = ["A", "B", "C"]
FDIR = "data/"

rule all:
    input:
        "result/out.html"

rule bwa_map:
    input:
        join(FDIR, "genome.fa"),
        "data/samples/{sample}.fastq"
    output:
        "result/mapped_reads/{sample}.bam"
    benchmark:
        "result/benchmarks/{sample}.bwa.benchmark.txt"
    threads: 8
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "(bwa mem  -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "result/mapped_reads/{sample}.bam"
    output:
        "result/sorted_reads/{sample}.bam"
    message: "Executing samtools_sort on the following files {input}..."
    log:
        "logs/mapped_reads/{sample}.log"
    shell:
        "samtools sort -T result/sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output} >{log} 2> {log}"

rule samtools_index:
    input:
        "result/sorted_reads/{sample}.bam"
    output:
        "result/sorted_reads/{sample}.bam.bai"
    message: "Executing samtools_index..."
    log:
        "logs/sorted_reads/{sample}.log"
    shell:
        "samtools index {input} >{log} 2> {log}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("result/sorted_reads/{sample}.bam",sample=SAMPLES),
        bai=expand("result/sorted_reads/{sample}.bam.bai",sample=SAMPLES)
    output:
        "result/calls/all.vcf"
    message: "Executing bcftools_call..."

    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

rule report:
    input:
        T1="result/calls/all.vcf",
        T2=expand("result/benchmarks/{sample}.bwa.benchmark.txt", sample=SAMPLES)
    output:
        "result/out.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as f:
            n_calls = sum(1 for line in f if not line.startswith("#"))

        report("""
        An example workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        Benchmark results for BWA can be found in the tables T2_.
        """, output[0], **input)

