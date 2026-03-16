# 05_select_14bacteria.r — Select bacteria with coverage >= 100 reads in > 50% of individuals
# Input:  OTUtable_Top50_with_metadata.tsv (from 03_otu_meta_top50.r)
# Output: pairs_14bacteria_samples_ge100.tsv (sample-bacteria pairs with >= 100 reads)
#
# Usage: Rscript 05_select_14bacteria.r config.sh

# --- Parse config.sh ---
args     <- commandArgs(trailingOnly = TRUE)
cfg_file <- ifelse(length(args) > 0, args[1],
                   file.path(dirname(commandArgs()[4]), "config.sh"))
cfg_lines <- readLines(cfg_file)
cfg_lines <- cfg_lines[grepl("^[A-Z_]+=", cfg_lines)]
for (line in cfg_lines) {
    parts <- strsplit(line, "=")[[1]]
    assign(parts[1], gsub('"', '', paste(parts[-1], collapse = "=")))
}

# --- Paths from config ---
otu_top50_file <- file.path(RESULTS_DIR, "OTUtable_Top50_with_metadata.tsv")
pairs_out <- PAIRS_FILE

MIN_READS    <- 100    # minimum reads per individual-bacteria pair
MIN_PREV     <- 50     # minimum prevalence (% of individuals) to select a bacteria
N_META_ROWS  <- 6      # number of metadata rows at the top of the OTU table

# --- Load OTU table ---
otu <- read.delim(otu_top50_file, sep = "\t", header = TRUE, check.names = FALSE)
otu <- otu[-(1:N_META_ROWS), ]   # remove metadata rows, keep species only

total_individuals <- ncol(otu) - 1
total_bacteria    <- nrow(otu)

# --- Convert counts to numeric ---
otu_num <- apply(otu[, -1], 2, function(x) as.numeric(as.character(x)))
rownames(otu_num) <- otu[, 1]
colnames(otu_num) <- colnames(otu)[-1]

# --- Per-bacteria: count individuals with >= MIN_READS ---
n_ind_per_bact <- rowSums(otu_num >= MIN_READS, na.rm = TRUE)
df_bacteria <- data.frame(
    Bacteria           = rownames(otu_num),
    n_individuals      = as.numeric(n_ind_per_bact),
    pct_individuals    = 100 * n_ind_per_bact / total_individuals,
    row.names          = NULL
)

# --- Select bacteria present in > MIN_PREV% of individuals ---
selected <- df_bacteria$Bacteria[df_bacteria$pct_individuals > MIN_PREV]
cat("Bacteria selected (prevalence >", MIN_PREV, "%):\n")
print(sort(selected))
cat("Total:", length(selected), "\n\n")

# --- Build sample-bacteria pairs with >= MIN_READS ---
sub <- otu_num[selected, , drop = FALSE]
idx <- which(sub >= MIN_READS, arr.ind = TRUE)

pairs <- data.frame(
    Bacteria = rownames(sub)[idx[, 1]],
    Sample   = colnames(sub)[idx[, 2]],
    CountOTU = sub[idx],
    row.names = NULL
)

write.table(pairs, pairs_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Pairs saved:", pairs_out, "\n")
cat("Total pairs:", nrow(pairs), "\n")
