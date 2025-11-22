##############################################
# Aspergillus ruber BGC + Genome QC Pipeline
# Snakemake DSL2 — Reproducible Bioinformatics
##############################################

import os

configfile: "config.yaml"

WORKDIR      = config["workdir"]
ASSEMBLY     = WORKDIR + "/genome/A_ruber_genome.fna"
PROTEINS     = WORKDIR + "/genome/A_ruber_proteins.faa"
GFF          = WORKDIR + "/genome/A_ruber.gff3"
BUSCO_OUT    = WORKDIR + "/qc/busco"
BLOB_OUT     = WORKDIR + "/qc/blobtoolkit"
BLASTDB      = WORKDIR + "/blastdb/A_ruber_db"
BLAST_OUT    = WORKDIR + "/blast/results.tsv"
ANTISMASH_OUT= WORKDIR + "/antismash"
HMMER_OUT    = WORKDIR + "/hmmer/pfam_scan.tsv"
PANNZER_OUT  = WORKDIR + "/pannzer"
CLUST_CORE   = WORKDIR + "/core_proteins/core.faa"
ALIGN_OUT    = WORKDIR + "/msa/core.aln.fasta"
TREE_OUT     = WORKDIR + "/phylogeny/core.tree"

##############################################
# Rule: Final Output
##############################################

rule all:
    input:
        ASSEMBLY,
        PROTEINS,
        BUSCO_OUT + "/short_summary.json",
        BLOB_OUT + "/plots/blobplot.png",
        BLAST_OUT,
        ANTISMASH_OUT + "/index.html",
        HMMER_OUT,
        PANNZER_OUT + "/functional_annotations.tsv",
        TREE_OUT

##############################################
#Step 1: Genome Download
##############################################

rule download_genome:
    output:
        WORKDIR + "/genome/ncbi_dataset.zip"
    params:
        accession=config["assembly_id"]
    shell:
        """
        mkdir -p {WORKDIR}/genome
        datasets download genome accession {params.accession} \
        --include genome,gff3,protein --filename {output}
        """

##############################################
#Step 2: Extract Genome & Proteins
##############################################

rule extract_genome:
    input:
        WORKDIR + "/genome/ncbi_dataset.zip"
    output:
        ASSEMBLY,
        PROTEINS,
        GFF
    shell:
        """
        mkdir -p {WORKDIR}/genome/extracted
        unzip -o {input} -d {WORKDIR}/genome/extracted

        find {WORKDIR}/genome/extracted -name "*.fna" | head -1 > tmp_fna
        find {WORKDIR}/genome/extracted -name "*.faa" | head -1 > tmp_faa
        find {WORKDIR}/genome/extracted -name "*.gff" | head -1 > tmp_gff

        cp $(cat tmp_fna) {output[0]}
        cp $(cat tmp_faa) {output[1]}
        cp $(cat tmp_gff) {output[2]}

        rm tmp_fna tmp_faa tmp_gff
        """

##############################################
#Step 3: BUSCO — Genome Completeness
##############################################

rule busco:
    input:
        ASSEMBLY
    output:
        BUSCO_OUT + "/short_summary.json"
    conda:
        "envs/busco.yaml"
    params:
        lineage=config["busco_lineage"],
        threads=config["threads"]
    shell:
        """
        mkdir -p {BUSCO_OUT}
        busco -i {input} \
              -l {params.lineage} \
              -m genome \
              -c {params.threads} \
              -o A_ruber_busco \
              --out_path {BUSCO_OUT}
        """

##############################################
#Step 4: BlobToolKit — Taxonomic Profiling
##############################################

rule blobtoolkit:
    input:
        ASSEMBLY
    output:
        BLOB_OUT + "/plots/blobplot.png"
    conda:
        "envs/blobtoolkit.yaml"
    params:
        threads=config["threads"],
        nr_db=config["nr_db"],
        nodes=config["tax_nodes"],
        names=config["tax_names"]
    shell:
        """
        mkdir -p {BLOB_OUT}

        blobtools create --fasta {input} \
                         --out {BLOB_OUT}/A_ruber_blobdb

        diamond blastx -q {input} \
                       -d {params.nr_db} \
                       -o {BLOB_OUT}/hits.tsv \
                       --threads {params.threads}

        blobtools taxify --hits {BLOB_OUT}/hits.tsv \
                         --nodes {params.nodes} \
                         --names {params.names} \
                         --out {BLOB_OUT}/taxified

        blobtools plot --out {BLOB_OUT}/plots \
                       {BLOB_OUT}/A_ruber_blobdb
        """

##############################################
#Step 5: Make Local BLAST Database
##############################################

rule make_blast_db:
    input:
        PROTEINS
    output:
        BLASTDB + ".pin"
    params:
        db=BLASTDB
    shell:
        """
        mkdir -p {WORKDIR}/blastdb
        makeblastdb -in {input} -dbtype prot -out {params.db}
        """

##############################################
#Step 6: BLASTP — Identify Secondary Metabolite Proteins
##############################################

rule blastp:
    input:
        query=config["sm_protein_queries"],
        db=BLASTDB + ".pin"
    output:
        BLAST_OUT
    conda:
        "envs/blast.yaml"
    params:
        db=BLASTDB,
        threads=config["threads"]
    shell:
        """
        mkdir -p {WORKDIR}/blast
        blastp -query {input.query} \
               -db {params.db} \
               -out {output} \
               -evalue 1e-10 \
               -num_threads {params.threads} \
               -outfmt 6
        """

##############################################
#Step 7: antiSMASH — BGC Detection
##############################################

rule antismash:
    input:
        ASSEMBLY,
        PROTEINS,
        GFF
    output:
        ANTISMASH_OUT + "/index.html"
    conda:
        "envs/antismash.yaml"
    params:
        threads=config["threads"]
    shell:
        """
        mkdir -p {ANTISMASH_OUT}
        antismash {input[0]} \
            --genefinding-gff3 {input[2]} \
            --cb-general --cb-knownclusters \
            --cb-subclusters --pfam2go \
            --asf --tta --taxon fungi \
            --output-dir {ANTISMASH_OUT} \
            --cpus {params.threads}
        """

##############################################
#Step 8: HMMER — Pfam Domain Search
##############################################

rule hmmer_scan:
    input:
        PROTEINS
    output:
        HMMER_OUT
    conda:
        "envs/hmmer.yaml"
    params:
        pfam=config["pfam_db"],
        threads=config["threads"]
    shell:
        """
        mkdir -p {WORKDIR}/hmmer
        hmmscan --cpu {params.threads} \
                --domtblout {output} \
                {params.pfam} \
                {input}
        """

##############################################
#Step 9: PANNZER2 Functional Annotation
##############################################

rule pannzer:
    input:
        PROTEINS
    output:
        PANNZER_OUT + "/functional_annotations.tsv"
    shell:
        """
        mkdir -p {PANNZER_OUT}
        pannzer2.py -i {input} \
                    -o {PANNZER_OUT}
        """

##############################################
#Step 10: Extract Core Proteins for Phylogeny
##############################################

rule extract_core:
    input:
        ANTISMASH_OUT + "/index.html"
    output:
        CLUST_CORE
    shell:
        """
        python scripts/extract_core_proteins.py \
            --antismash {ANTISMASH_OUT} \
            --out {output}
        """

##############################################
#Step 11: COBALT Alignment
##############################################

rule cobalt_align:
    input:
        CLUST_CORE
    output:
        ALIGN_OUT
    shell:
        """
        cobalt -i {input} -o {output}
        """

##############################################
#Step 12: Phylogeny — Distance Tree
##############################################

rule build_tree:
    input:
        ALIGN_OUT
    output:
        TREE_OUT
    shell:
        """
        fasttree -wag {input} > {output}
        """
