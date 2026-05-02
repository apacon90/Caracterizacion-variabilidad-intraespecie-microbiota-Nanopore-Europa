# =============================================================================
# config.sh — Pipeline configuration
# Edit paths before running any script.
# =============================================================================

# --- Base directories ---
BASE_DIR="/path/to/rrn_operon_EU_datasets"
PROJECT_DIR="/path/to/GROND_nr_Rep"

# --- Phase 0 ---
OUT_TRIMMED="${PROJECT_DIR}/phase0_data/trimmed_reads"

# --- Phase 1 ---
MAPPINGS_DIR="${PROJECT_DIR}/phase1_data/mappings"
SPECIES_READS_DIR="${PROJECT_DIR}/phase1_data/species_reads"
RRN_CONSENSUS_DIR="${PROJECT_DIR}/phase1_data/Flye"
BLAST_DB="${PROJECT_DIR}/phase1_data/blast_db/refseq207nr_nrRep_db"
REF_DIR="${PROJECT_DIR}/phase1_data/results_phase1/ref_per_species"
POOLS_DIR="${PROJECT_DIR}/phase1_data/results_phase1/pool_mappings"
VCF_RAW_DIR="${PROJECT_DIR}/phase1_data/results_phase1/vcf_raw"
VCF_FILTERED_DIR="${PROJECT_DIR}/phase1_data/results_phase1/vcf_filtered"

# --- Phase 2 ---
FASTA_ALL="${PROJECT_DIR}/phase2_data/fasta_all"
FASTA_N="${PROJECT_DIR}/phase2_data/fasta_n100_or_max"
AUX_DIR="${PROJECT_DIR}/phase2_data/aux"
IND_BAM_DIR="${PROJECT_DIR}/phase2_data/results_phase2/individual_mappings"
IND_VCF_DIR="${PROJECT_DIR}/phase2_data/results_phase2/vcf_raw"
IND_VCF_FILT_DIR="${PROJECT_DIR}/phase2_data/results_phase2/vcf_filtered"
TSV_DIR="${PROJECT_DIR}/phase2_data/variants_tsv"

# --- Reference files ---
DB_REF="${PROJECT_DIR}/reference/refseq207nr_nrRep.fna"
DB_REF_FASTA="${PROJECT_DIR}/reference/refseq207nr_nrRep.fna"
DB_MMI="${PROJECT_DIR}/reference/refseq207nr_nrRep.mmi"
TAX_FILE="${PROJECT_DIR}/reference/refseq207nr_taxRep.tsv"
CODE2NAME_FILE="${PROJECT_DIR}/reference/code2name.tsv"

# --- Metadata and results ---
META_FILE="${PROJECT_DIR}/metadata/sample_metadata.tsv"
TOP50_TSV="${PROJECT_DIR}/metadata/Top50_species_abundance_median.tsv"
SELECT="${PROJECT_DIR}/metadata/select_sp.tsv"
RESULTS_DIR="${PROJECT_DIR}/metadata"
PAIRS_FILE="${PROJECT_DIR}/phase2_data/pairs_14bacteria_samples_ge100.tsv"

# --- Tools ---
BARRNAP="/path/to/barrnap/bin/barrnap"

# --- Trimming ---
TRIM5=100
READS_SUFFIX=".trim100.fa"

# --- Mapping filters ---
THREADS=40
MIN_ALN=3000
MIN_COV=0.70
MIN_PID=0.85

# --- bcftools parameters ---
MAPQ_MIN=20
BQ_MIN=10
MAX_BQ=35
MAX_DEPTH=500
PLOIDY=1
P_PRIOR="1e-3"
ANNOT="FORMAT/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR"

# --- VCF quality filters ---
QUAL_SNP=100
QUAL_INDEL=150

# --- Parallelization ---
JOBS=32
THREADS_PER_JOB=8
THREADS_PER_MAP=8
MAX_JOBS=4
