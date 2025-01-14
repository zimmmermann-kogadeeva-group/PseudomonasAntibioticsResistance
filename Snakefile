
from pathlib import Path
import re


configfile: "config.json"


rule all:
    input:
        "Output/QC/multiqc_report.html",
        "Output/DESeq2/dds.rds",
        "Output/Proteomics/normalised.csv",
        "Output/Snippy/.all",

rule fastqc:
    input:
        "Data/Transcriptomics/{filename}.fastq.gz"
    output:
        "Output/QC/{filename}_fastqc.html"
    group:
        "small"
    conda:
        "Envs/qc.yaml"
    shell:
        "fastqc {input} -o Output/QC/"

rule multiqc:
    input:
        expand(
            "Output/QC/{f}_L{lane}_fastqc.html", 
            f=config["transcriptomics"], 
            lane=[1, 2]
        )
    output:
        "Output/QC/multiqc_report.html",
        directory("Output/QC/multiqc_data/")
    conda:
        "Envs/qc.yaml"
    shell:
        "multiqc Output/QC/ -o Output/QC/"

rule trimmomatic:
    input:
        adapter="Data/illumina_adapters.fasta",
        fastq="Data/{path_to_file}.fastq.gz"
    output:
        temp("Output/TrimmedReads/{path_to_file}.fastq.gz")
    params:
        trim_opts = lambda wildcards, input :  
            f"ILLUMINACLIP:{input.adapter}:2:30:10 "
            "LEADING:3 "
            "TRAILING:3 "
            "SLIDINGWINDOW:4:15 "
            "MINLEN:36"
    conda:
        "Envs/trimmomatic.yaml"
    resources:
        runtime = 90
    shell:
        "trimmomatic SE -phred33 {input.fastq} {output} {params.trim_opts}"

rule merge:
    input:
        l1="Output/TrimmedReads/{path_to_file}_L1.fastq.gz",
        l2="Output/TrimmedReads/{path_to_file}_L2.fastq.gz"
    output:
        merged=temp("Output/Merged/{path_to_file}.fastq.gz")
    shell:
        "cat {input.l1} {input.l2} > {output.merged}"
        
rule salmon_index:
    input:
        "Data/reference.fasta.gz"
    output:
        directory("Output/SalmonIndex")
    conda:
        "Envs/salmon.yaml"
    shell:
        "salmon index -t {input} -i {output}"

rule salmon_quant:
    input:
        index="Output/SalmonIndex",
        transcript="Output/Merged/Transcriptomics/{filename}.fastq.gz"
    output:
        "Output/Quants/{filename}/quant.sf"
    conda:
        "Envs/salmon.yaml"
    params:
        outdir = lambda wildcards, output : Path(output[0]).parent
    shell:
        "salmon quant -i {input.index} -l A -p 1 --gcBias --validateMappings "
        "-r {input.transcript} -o {params.outdir}"

rule run_deseq2:
    input:
        quant_files = expand("Output/Quants/{f}/quant.sf", f=config["transcriptomics"]),
        annot = "Data/annotation.tsv"
    output:
        "Output/DESeq2/dds.rds",
        "Output/DESeq2/dds_results.csv",
        "Output/DESeq2/counts.csv",
        "Output/DESeq2/counts_normalised.csv",
        "Output/DESeq2/vst.csv",
        "Output/DESeq2/rlog.csv"
    conda:
        "Envs/deseq2.yaml"
    params:
        output_dir = "Output/DESeq2/"
    resources:
        mem_mb = "5G"
    script:
        "Scripts/run_deseq2.R"

rule get_proteomics:
    input:
        expand("Data/Proteomics/{f}.txt", f=config["proteomics"])
    output:
        "Output/Proteomics/raw.csv",
        "Output/Proteomics/normalised.csv",
        "Output/Proteomics/lfc.csv"
    params:
        output_dir="Output/Proteomics/"
    script:
        "Scripts/get_proteomics.py"

rule prokka:
    input:
        "Data/{path_to_file}.fasta",
        "Data/PAO1.gbk"
    output:
        directory("Output/Prokka/{path_to_file}/"),
        "Output/Prokka/{path_to_file}/prokka.gbk"
    container:
        # prokka 1.14.6
        "docker://staphb/prokka@sha256:6bb2522e077ef08a8be7a3856fe80372ede3b832becba0728e1ddbe83d89a042"
    params:
        opts = (
            "--evalue 0.000001 "
            "--force "
            "--increment 1 "
            "--gffver 3 "
            "--mincontig 200 "
            "--kingdom Bacteria "
            "--gcode 11 "
            "--genus Pseudomonas "
            "--species aeruginosa "
            "--prefix prokka"
        )
    resources:
        mem_mb = "5G",
        runtime = lambda wildcards, attempt: 60 * attempt
    shell:
        "prokka --proteins '{input[1]}' {params.opts} --outdir {output[0]} {input[0]}"

rule run_snippy_pe:
    input:
        r1 = "Data/WGS/TrimmedReads/{strain}_{cond}_R1.fastq.gz",
        r2 = "Data/WGS/TrimmedReads/{strain}_{cond}_R2.fastq.gz",
        ref="Output/Prokka/WGS/Assemblies/{ref_strain}_P/prokka.gbk"
    output:
        directory("Output/Snippy/{strain}_{cond}_vs_{ref_strain}_P")
    resources:
        cpus=16,
        mem_mb=lambda wildcards, attempt: f"{10 + 5 * attempt}G",
        runtime=lambda wildcards, attempt: 360 * attempt
    conda:
        "Envs/snippy.yaml"
    shadow: 
        "full"
    params:
        opts = "--force --report --cpus 16"
    shell:
        "snippy {params.opts} --outdir {output[0]} "
        "--ref {input.ref} --R1 {input.r1} --R2 {input.r2}"

rule run_snippy_se:
    input:
        se = "Data/WGS/TrimmedReads/{strain}_{cond}.fastq.gz",
        ref="Output/Prokka/WGS/Assemblies/{ref_strain}_P/prokka.gbk"
    output:
        directory("Output/Snippy/{strain}_{cond}_vs_{ref_strain}_P")
    resources:
        cpus=16,
        mem_mb=lambda wildcards, attempt: f"{10 + 5 * attempt}G",
        runtime=lambda wildcards, attempt: 360 * attempt
    conda:
        "Envs/snippy.yaml"
    shadow: 
        "full"
    params:
        opts = "--force --report --cpus 16"
    shell:
        "snippy {params.opts} --outdir {output[0]} --ref {input.ref} --se {input.se}"

rule all_snippy:
    input:
        expand(
            "Output/Snippy/{strain}_{cond}_vs_{ref_strain}_P", 
            strain=config["snippy"],
            cond=["P", "M", "C"],
            ref_strain=config["snippy"]
        )
    output:
        touch("Output/Snippy/.all")

