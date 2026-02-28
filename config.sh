# =============================================================================
# config.sh — Pipeline configuration
# Edit this file before running any script.
# =============================================================================

# --- Directories ---
BASE_DIR="/path/to/rrn_operon_EU_datasets"   # raw nanopore reads (subdirs per country)
OUT_TRIMMED="/path/to/out_trimmed"           # output of 01_trimm.sh
MAPPINGS_DIR="/path/to/mappings"             # output of 02_map_filter.sh
RESULTS_DIR="/path/to/results"               # output of 03 and onwards

# --- Reference files ---
DB_REF="/path/to/refseq207nr_nrRep.fna.gz"  # RefSeq207nr FASTA
DB_MMI="/path/to/refseq207nr_nrRep.mmi"     # minimap2 index (created if missing)
TAX_FILE="/path/to/refseq207nr_taxRep.tsv"  # ID -> taxonomic lineage dictionary
META_FILE="/path/to/sample_metadata.tsv"    # sample metadata (must contain "Sample" column)

# --- Trimming ---
TRIM5=100         # nucleotides to remove from 5' end

# --- Mapping filters ---
THREADS=40
MIN_ALN=3000      # minimum alignment length (PAF col11)
MIN_COV=0.70      # minimum reference coverage (col11/col7)
MIN_PID=0.85      # minimum sequence identity  (col10/col11)
