# 03_otu_meta_top50.r — Build OTU table and select Top50 species by median abundance
# Input:  *.keep.paf files (from 02_map_filter.sh), taxonomy file, sample metadata
# Output: OTU table with metadata, Top50 species table, OTU table filtered to Top50
#
# Usage: Rscript 03_otu_meta_top50.r config.sh

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
tax_file      <- TAX_FILE
paf_dir       <- MAPPINGS_DIR
meta_file     <- META_FILE
otu_out       <- file.path(RESULTS_DIR, "OTUtable_with_metadata.tsv")
top50_out     <- file.path(RESULTS_DIR, "Top50_species_abundance_median.tsv")
otu_top50_out <- file.path(RESULTS_DIR, "OTUtable_Top50_with_metadata.tsv")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Helper: normalize species name ---
norm_sp <- function(x) {
    x <- gsub("^s__", "", x)
    x <- gsub("['\"]", "", x)
    x <- gsub("\\[|\\]", "", x)
    x <- gsub("[[:space:]]+|/", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
}

# --- 1. Load taxonomy dictionary ID -> species ---
tax         <- read.table(tax_file, sep = "\t", header = FALSE, quote = "",
                          stringsAsFactors = FALSE)
colnames(tax) <- c("ID", "Lineage")
tax$Species <- sapply(tax$Lineage, function(lin) {
    parts <- strsplit(lin, "\\|")[[1]]
    sp    <- parts[grep("^s__", parts)]
    if (length(sp) == 0) return(NA_character_)
    norm_sp(sp)
})
id2sp <- setNames(tax$Species, tax$ID)

# --- 2. Build OTU table from PAF files ---
paf_files <- list.files(paf_dir, pattern = "\\.keep\\.paf$", full.names = TRUE)
if (length(paf_files) == 0) stop("No .keep.paf files found in: ", paf_dir)

otu <- NULL

for (f in paf_files) {
    sample <- sub("\\.keep\\.paf$", "", basename(f))
    cat("Processing:", sample, "\n")

    paf    <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                         quote = "", comment.char = "")
    qname  <- paf[[1]];  tname  <- paf[[6]]
    nmatch <- as.numeric(paf[[10]]); alen <- as.numeric(paf[[11]]); tlen <- as.numeric(paf[[7]])

    sp   <- id2sp[tname]
    keep <- !is.na(sp)
    if (!any(keep)) next

    df <- data.frame(q = qname[keep], sp = sp[keep],
                     pid = ifelse(alen[keep] > 0, nmatch[keep] / alen[keep], 0),
                     cov = ifelse(tlen[keep] > 0, alen[keep]  / tlen[keep], 0),
                     stringsAsFactors = FALSE)

    # Keep best alignment per read (highest pid, tiebreak by cov)
    df   <- df[order(df$q, -df$pid, -df$cov), ]
    best <- df[!duplicated(df$q), ]

    counts    <- table(best$sp)
    df_counts <- data.frame(Species = names(counts), Count = as.integer(counts))

    if (is.null(otu)) {
        otu <- df_counts; colnames(otu)[2] <- sample
    } else {
        otu <- merge(otu, df_counts, by = "Species", all = TRUE)
        colnames(otu)[ncol(otu)] <- sample
    }
}

otu[is.na(otu)] <- 0
otu <- otu[order(otu$Species), ]

# --- 3. Add sample metadata ---
meta        <- read.table(meta_file, sep = "\t", header = TRUE, quote = "",
                          stringsAsFactors = FALSE, check.names = FALSE)
samples_otu <- colnames(otu)[-1]
meta        <- meta[match(samples_otu, meta$Sample), , drop = FALSE]
stopifnot(identical(meta$Sample, samples_otu))

meta_rows           <- data.frame(`#Name` = colnames(meta)[-1],
                                  t(meta[, -1, drop = FALSE]), check.names = FALSE)
colnames(meta_rows) <- c("#Name", samples_otu)

otu2              <- otu; colnames(otu2)[1] <- "#Name"
otu_with_meta     <- rbind(meta_rows, otu2)

write.table(otu_with_meta, file = otu_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("OTU table saved:", otu_out, "\n")

# --- 4. Select Top50 species by median abundance ---
n_meta      <- nrow(meta_rows)
species_tbl <- otu_with_meta[-(1:n_meta), , drop = FALSE]
species_tbl[, -1] <- lapply(species_tbl[, -1], as.numeric)

medians     <- apply(species_tbl[, -1], 1, median, na.rm = TRUE)
tabla_med   <- data.frame(Species = species_tbl[[1]], Median = medians)
tabla_med   <- tabla_med[!grepl("sp\\.", tabla_med$Species), ]   # remove unresolved species
top50       <- head(tabla_med[order(-tabla_med$Median), ], 50)

write.table(top50, file = top50_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Top50 saved:", top50_out, "\n")

# --- 5. Export OTU table filtered to Top50 ---
sp_top50          <- species_tbl[match(top50$Species, species_tbl[[1]]), , drop = FALSE]
sp_top50          <- sp_top50[!is.na(sp_top50[[1]]), , drop = FALSE]
otu_top50_meta    <- rbind(meta_rows, sp_top50)

write.table(otu_top50_meta, file = otu_top50_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("OTU Top50 saved:", otu_top50_out, "\n")
