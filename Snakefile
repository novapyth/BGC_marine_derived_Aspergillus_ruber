# Snakemake pipeline for A_ruber BGC analysis
import os
from snakemake.utils import min_version
min_version("6.0")

configfile: "config.yaml"

WORKDIR = config["workdir"].format(**config)
GENOME_DIR = config["genome_dir"].format(**config)
GENOME = config["genome_fasta"].format(**config)
PROTEINS = config["proteins"].format(**config)
GFF = config["gff"].format(**config)
QUERY = config["query_proteins"].format(**config)
PKS_QUERIES = config["pks_nrps_queries"].format(**config)
CORE = config["core_proteins"].format(**config)
THREADS = config.get("threads", 8)

rule all:
    input:
        expand(WORKDIR + "/blast/blast_results.tsv"),
        expand(WORKDIR + "/antismash/antismash_output/index.html"),
        expand(WORKDIR + "/hmmer/pfam.domtblout"),
        expand(WORKDIR + "/pannzer/pannzer_output.tsv"),
        expand(WORKDIR + "/core/core_proteins.tree")

#############################################
# Step 1: Download with ncbi-datasets
#############################################
rule download_genome:
    output:
        genome_tar=WORKDIR + "/genome/ncbi_dataset.zip"
    conda: "envs/tools.yaml"
    shell:
        """
        mkdir -p {WORKDIR}/genome
        cd {WORKDIR}/genome
        datasets download genome accession {config[assembly_accession]} --include genome,protein,gff3,seq-report --filename ncbi_dataset.zip
        """

rule extract_genome:
    input:
        WORKDIR + "/genome/ncbi_dataset.zip"
    output:
        GENOME, PROTEINS, GFF = WORKDIR + "/genome/A_ruber_genome.fna", WORKDIR + "/genome/A_ruber_proteins.faa", WORKDIR + "/genome/A_ruber.gff3"
    shell:
        """
        unzip -o {input} -d {WORKDIR}/genome/dataset
        cp {WORKDIR}/genome/dataset/ncbi_dataset/data/{config[assembly_accession]}/**/*.fna {output.GENOME} 2>/dev/null || true
        cp {WORKDIR}/genome/dataset/ncbi_dataset/data/{config[assembly_accession]}/**/*.faa {output.PROTEINS} 2>/dev/null || true
        cp {WORKDIR}/genome/dataset/ncbi_dataset/data/{config[assembly_accession]}/**/*.gff* {output.GFF} 2>/dev/null || true
        """

#############################################
# Step 2: Make BLAST DB and run BLASTP
#############################################
rule make_blast_db:
    input:
        PROTEINS
    output:
        WORKDIR + "/blast/A_ruber_prot_db.pin"
    conda: "envs/tools.yaml"
    shell:
        """
        mkdir -p {WORKDIR}/blast
        makeblastdb -in {input} -dbtype prot -out {WORKDIR}/blast/A_ruber_prot_db
        """

rule blastp:
    input:
        query=QUERY,
        db=WORKDIR + "/blast/A_ruber_prot_db"
    output:
        WORKDIR + "/blast/blast_results.tsv"
    params:
        evalue=config.get("blast_evalue", "1e-10")
    threads: THREADS
    conda: "envs/tools.yaml"
    shell:
        """
        mkdir -p {WORKDIR}/blast
        blastp -query {input.query} -db {input.db} -evalue {params.evalue} -out {output} -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs" -num_threads {threads}
        """

#############################################
# Step 3: antiSMASH
#############################################
rule run_antismash:
    input:
        genome=GENOME,
        gff=GFF
    output:
        WORKDIR + "/antismash/antismash_output/index.html"
    conda: "envs/tools.yaml"
    threads: THREADS
    shell:
        """
        mkdir -p {WORKDIR}/antismash
        antismash --input-type nucl {input.genome} --gff3 {input.gff} --taxon fungi --strict --output-dir {WORKDIR}/antismash/antismash_output
        """

#########################################
