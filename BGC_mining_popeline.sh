#!/bin/bash

#############################################
# 0. PREPARE WORKING DIRECTORY
#############################################
mkdir -p A_ruber_pipeline && cd A_ruber_pipeline


#############################################
# 1. DOWNLOAD GENOME & PROTEINS USING NCBI DATASETS
#############################################
echo "[STEP 1] Downloading Aspergillus ruber genome + proteins..."
mkdir -p genome && cd genome

datasets download genome accession GCA_000600275.1 \
    --include genome,protein,gff3,seq-report

unzip ncbi_dataset.zip -d dataset

cp dataset/ncbi_dataset/data/GCA_000600275.1/*.fna A_ruber_genome.fna
cp dataset/ncbi_dataset/data/GCA_000600275.1/*.faa A_ruber_proteins.faa
cp dataset/ncbi_dataset/data/GCA_000600275.1/*.gff A_ruber.gff3

cd ..


#############################################
# 2. BUILD LOCAL BLAST DATABASE + BLASTP SCREENING
#############################################
echo "[STEP 2] Making BLAST database and screening proteins..."
mkdir -p blast_analysis && cd blast_analysis

makeblastdb \
  -in ../genome/A_ruber_proteins.faa \
  -dbtype prot \
  -out A_ruber_prot_db

blastp \
  -query query.faa \
  -db A_ruber_prot_db \
  -evalue 1e-10 \
  -out blast_results.tsv \
  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs"

cd ..


#############################################
# 3. RUN ANTISMASH 6 IN FUNGAL MODE
#############################################
echo "[STEP 3] Running antiSMASH for BGC prediction..."
mkdir -p antismash_run && cd antismash_run

antismash \
  --input-type nucl \
  ../genome/A_ruber_genome.fna \
  --gff3 ../genome/A_ruber.gff3 \
  --taxon fungi \
  --strict \
  --cb-general \
  --cb-subclusters \
  --clusterblast \
  --subclusterblast \
  --knownclusterblast \
  --active-site-finder \
  --mibig \
  --pfam2go \
  --tigrfam \
  --smcog-trees \
  --cassis \
  --tfbs \
  --genefinding-tool prodigal \
  --output-dir antismash_output

cd ..


#############################################
# 4. DETECT PKS/NRPS DOMAINS — BLAST + HMMER
#############################################
echo "[STEP 4] Detecting PKS/NRPS domains..."
mkdir -p pks_nrps && cd pks_nrps

blastp \
  -query pks_nrps_queries.faa \
  -db ../blast_analysis/A_ruber_prot_db \
  -evalue 1e-10 \
  -out pks_nrps_hits.tsv \
  -outfmt 6

hmmscan \
  --cpu 8 \
  --domtblout pfam.domtblout \
  /path/to/Pfam-A.hmm \
  ../genome/A_ruber_proteins.faa

cd ..


#############################################
# 5. FUNCTIONAL REANNOTATION USING PANNZER2
#############################################
echo "[STEP 5] Running PANNZER2 annotation..."
mkdir -p pannzer && cd pannzer

pannzer2.py \
  -i ../genome/A_ruber_proteins.faa \
  -o pannzer_output \
  --uniprot /path/to/uniprot_sprot.fasta

cd ..


#############################################
# 6. RE-ANNOTATE CORE BGC PROTEINS
#############################################
echo "[STEP 6] Reannotating core proteins..."
mkdir -p core_reannotation && cd core_reannotation

blastp \
  -query core_proteins.faa \
  -db nr \
  -remote \
  -qcov_hsp_perc 70 \
  -perc_identity 50 \
  -evalue 1e-5 \
  -out core_proteins_reannot.tsv \
  -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore stitle"


#############################################
# 7. MULTIPLE SEQUENCE ALIGNMENT + PHYLOGENY
#############################################
echo "[STEP 7] Multiple sequence alignment + phylogeny..."
cobalt \
  -i core_proteins.faa \
  -o core_proteins.aln.fasta

FastTree core_proteins.aln.fasta > core_proteins.tree

echo "Pipeline completed successfully!"
