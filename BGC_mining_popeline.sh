#!/usr/bin/env bash
#Name: BGC_mining_popeline.sh
# Version: V2
set -euo pipefail
IFS=$'\n\t'

####### BGC_mining_popeline.sh
# Complete end-to-end driver script for A_ruber BGC analysis
# - Downloads genome via ncbi-datasets
# - Runs BUSCO completeness (fungi_odb10)
# - Runs BlobToolKit taxonomic profiling (diamond + blobtools)
# - Builds BLAST DB and runs BLASTP queries
# - Runs antiSMASH (fungal mode strict)
# - Runs hmmscan (Pfam)
# - Optionally calls PANNZER2 via pannzer_client.py
# - Extracts core proteins from GFF (core_annotator.py)
# - Performs alignment + phylogeny (MAFFT/MAFFT or MUSCLE + FastTree)
#
# REQUIREMENTS (examples):
#   conda create -n bgc-tools -c conda-forge -c bioconda ncbi-datasets-cli ncbi-blast+ hmmer fasttree mafft muscle diamond busco blobtools antismash python=3.10
#
# Usage:
#   ./BGC_mining_popeline.sh
#
# Edit variables below as needed.

########### User-editable variables ############
ASSEMBLY="GCA_000600275.1"
WORKDIR="$(pwd)/A_ruber_pipeline"
THREADS=8
BUSCO_LINEAGE="fungi_odb10"
PFAM_HMM="/usr/local/db/Pfam-A.hmm"      # edit to point to local Pfam HMM database
DIAMOND_DB="/db/nr.dmnd"                 # edit to your diamond DB for blobtoolkit taxify step
NAMES_DMP="/db/names.dmp"                # taxdump files for blobtools taxify (if used)
NODES_DMP="/db/nodes.dmp"
PANNZER_API_ENDPOINT=""                  # optional: set to PANNZER2 API endpoint if using pannzer_client.py
QUERY_PROTEINS="queries/query.faa"       # path to your query proteins for BLASTP (create if missing)
PKS_QUERIES="queries/pks_nrps_queries.faa"
CORE_GFF_FILTER="core_genes_list.txt"    # optional list of gene IDs for core annotation (one per line)
##############################################

LOG="${WORKDIR}/pipeline.log"
mkdir -p "${WORKDIR}"
echo "Pipeline started at $(date)" | tee "${LOG}"

# Optional: activate conda env
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate bgc-tools

# helper
run_step(){
  echo -e "\n==== $1 ====" | tee -a "${LOG}"
  shift
  "$@" 2>&1 | tee -a "${LOG}"
}

# 1. Download genome using ncbi-datasets
run_step "Downloading genome ${ASSEMBLY}" bash -lc "
mkdir -p ${WORKDIR}/genome &&
cd ${WORKDIR}/genome &&
datasets download genome accession ${ASSEMBLY} --include genome,protein,gff3,seq-report --filename ncbi_dataset.zip
"

# 2. Extract genome/proteins/gff
run_step "Extracting files from ncbi dataset" bash -lc "
cd ${WORKDIR}/genome &&
unzip -o ncbi_dataset.zip -d extracted &&
# find first matching files and copy to stable names
fna=\$(find extracted -type f -name '*.fna' | head -1) &&
faa=\$(find extracted -type f -name '*.faa' | head -1) &&
gff=\$(find extracted -type f -iname '*.gff*' | head -1) &&
if [[ -z \"\$fna\" || -z \"\$faa\" ]]; then
  echo 'ERROR: genome FASTA or protein FASTA not found in downloaded dataset' >&2
  exit 2
fi
cp \"\$fna\" ${WORKDIR}/genome/A_ruber_genome.fna &&
cp \"\$faa\" ${WORKDIR}/genome/A_ruber_proteins.faa &&
if [[ -n \"\$gff\" ]]; then cp \"\$gff\" ${WORKDIR}/genome/A_ruber.gff3 || true; fi
"

# 3. BUSCO completeness
run_step "Running BUSCO (completeness)" bash -lc "
mkdir -p ${WORKDIR}/qc/busco &&
busco -i ${WORKDIR}/genome/A_ruber_genome.fna -l ${BUSCO_LINEAGE} -m genome -c ${THREADS} -o A_ruber_busco -f --out_path ${WORKDIR}/qc/busco
"

# 4. BlobToolKit pipeline (diamond + blobtools)
run_step "Running BlobToolKit (diamond + blobtools)" bash -lc "
mkdir -p ${WORKDIR}/qc/blobtoolkit &&
# create blobtools DB (sequence + basic scaffold info)
blobtools create --fasta ${WORKDIR}/genome/A_ruber_genome.fna --out ${WORKDIR}/qc/blobtoolkit/A_ruber_blobdb || true

# run DIAMOND BLASTX (recommended) - ensure DIAMOND DB path is set
if [[ -f \"${DIAMOND_DB}\" ]]; then
  diamond blastx -d ${DIAMOND_DB} -q ${WORKDIR}/genome/A_ruber_genome.fna -o ${WORKDIR}/qc/blobtoolkit/A_ruber.diamond.tsv -p ${THREADS} --sensitive
  blobtools add --db ${WORKDIR}/qc/blobtoolkit/A_ruber_blobdb --hits ${WORKDIR}/qc/blobtoolkit/A_ruber.diamond.tsv || true
else
  echo 'WARNING: DIAMOND DB not found at ${DIAMOND_DB}. Skipping diamond step.' >&2
fi

# taxify only if taxonomy dump files available
if [[ -f \"${NAMES_DMP}\" && -f \"${NODES_DMP}\" ]]; then
  blobtools taxify --hits ${WORKDIR}/qc/blobtoolkit/A_ruber.diamond.tsv --db ${NAMES_DMP} --nodes ${NODES_DMP} --out ${WORKDIR}/qc/blobtoolkit/taxify || true
fi

# generate view & plots (requires R and blobtools plotting deps)
blobtools view --out ${WORKDIR}/qc/blobtoolkit/view ${WORKDIR}/qc/blobtoolkit/A_ruber_blobdb || true
blobtools plot --out ${WORKDIR}/qc/blobtoolkit/plots ${WORKDIR}/qc/blobtoolkit/A_ruber_blobdb || true
"

# 5. Build BLAST DB of proteome
run_step "Make BLAST DB (proteins)" bash -lc "
mkdir -p ${WORKDIR}/blast &&
makeblastdb -in ${WORKDIR}/genome/A_ruber_proteins.faa -dbtype prot -out ${WORKDIR}/blast/A_ruber_prot_db
"

# 6. Run BLASTP for query proteins (user-supplied)
if [[ -f "${QUERY_PROTEINS}" ]]; then
  run_step "Running BLASTP of query proteins vs local proteome" bash -lc "
  blastp -query ${QUERY_PROTEINS} -db ${WORKDIR}/blast/A_ruber_prot_db -out ${WORKDIR}/blast/blast_results.tsv -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' -num_threads ${THREADS}
  "
else
  echo "NOTE: query proteins not found at ${QUERY_PROTEINS} — skipping BLASTP" | tee -a "${LOG}"
fi

# 7. Run antiSMASH (fungal mode, strict)
run_step "Running antiSMASH (fungal mode)" bash -lc "
mkdir -p ${WORKDIR}/antismash &&
antismash --input-type nucl ${WORKDIR}/genome/A_ruber_genome.fna --gff3 ${WORKDIR}/genome/A_ruber.gff3 --taxon fungi --strict --output-dir ${WORKDIR}/antismash -c ${THREADS} || true
"

# 8. HMMER scan vs Pfam
if [[ -f "${PFAM_HMM}" ]]; then
  run_step "Running hmmscan (Pfam)" bash -lc "
  mkdir -p ${WORKDIR}/hmmer &&
  hmmscan --cpu ${THREADS} --domtblout ${WORKDIR}/hmmer/pfam.domtblout ${PFAM_HMM} ${WORKDIR}/genome/A_ruber_proteins.faa
  "
else
  echo "WARNING: PFAM HMM not found at ${PFAM_HMM}. Skipping hmmscan" | tee -a "${LOG}"
fi

# 9. PKS/NRPS query BLAST (optional)
if [[ -f "${PKS_QUERIES}" ]]; then
  run_step "Searching for PKS/NRPS homologs (BLASTP)" bash -lc "
  blastp -query ${PKS_QUERIES} -db ${WORKDIR}/blast/A_ruber_prot_db -out ${WORKDIR}/blast/pks_nrps_hits.tsv -evalue 1e-10 -outfmt '6 qseqid sseqid pident length qcovs evalue bitscore' -num_threads ${THREADS}
  "
else
  echo "NOTE: PKS/NRPS queries not found at ${PKS_QUERIES} — skipping PKS/NRPS BLAST" | tee -a "${LOG}"
fi

# 10. PANNZER2 functional re-annotation (optional)
if [[ -n "${PANNZER_API_ENDPOINT}" ]]; then
  run_step "Running PANNZER2 client (remote API)" bash -lc "
  mkdir -p ${WORKDIR}/pannzer
  python3 scripts/pannzer_client.py --input ${WORKDIR}/genome/A_ruber_proteins.faa --endpoint '${PANNZER_API_ENDPOINT}' --out ${WORKDIR}/pannzer/pannzer_output.tsv
  "
else
  echo "PANNZER API endpoint not set. If you have pannzer2 locally installed, you can run it manually or set PANNZER_API_ENDPOINT variable." | tee -a "${LOG}"
fi

# 11. Extract core proteins from GFF using core_annotator.py
if [[ -f "${WORKDIR}/genome/A_ruber.gff3" ]]; then
  run_step "Extracting core proteins from GFF (core_annotator.py)" bash -lc "
  mkdir -p ${WORKDIR}/core
  python3 scripts/core_annotator.py --gff ${WORKDIR}/genome/A_ruber.gff3 --faa ${WORKDIR}/genome/A_ruber_proteins.faa --out ${WORKDIR}/core/core_proteins.faa --genes ${CORE_GFF_FILTER}
  "
else
  echo "No GFF found - skipping core extraction" | tee -a "${LOG}"
fi

# 12. Align and tree (MAFFT + FastTree or MUSCLE + FastTree)
if [[ -f "${WORKDIR}/core/core_proteins.faa" ]]; then
  run_step "Alignment (MAFFT) + Tree (FastTree)" bash -lc "
  mkdir -p ${WORKDIR}/tree
  mafft --auto ${WORKDIR}/core/core_proteins.faa > ${WORKDIR}/tree/core_alignment.faa
  FastTree -wag ${WORKDIR}/tree/core_alignment.faa > ${WORKDIR}/tree/core_tree.nwk
  "
else
  echo "No core protein FASTA found to align. Skipping phylogeny." | tee -a "${LOG}"
fi

echo -e "\nPipeline finished at $(date)" | tee -a "${LOG}"
