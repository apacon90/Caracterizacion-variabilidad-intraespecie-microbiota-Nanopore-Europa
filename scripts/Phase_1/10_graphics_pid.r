# 10_graphics_pid.r — Select best reference per species from BLAST results and plot % identity
# Input:  RESULTS_DIR/blast_top5_taxrep.tsv  (from 09_blast_match_top5.sh)
# Output: RESULTS_DIR/select_sp.tsv
#         RESULTS_DIR/plot_sp_ordered_bypid.png
#         RESULTS_DIR/plot_sp_by_family.png
#
# Usage: Rscript 11_graphics_pid.r config.sh

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

# --- Build derived paths ---
TABLES_DIR  <- file.path(RESULTS_DIR, "tables")
FIGURES_P1  <- file.path(RESULTS_DIR, "figures", "phase1")
FIGURES_P2  <- file.path(RESULTS_DIR, "figures", "phase2")

# --- Paths from config ---
infile  <- file.path(TABLES_DIR, "blast_top5_taxrep.tsv")
outfile <- file.path(RESULTS_DIR, "select_sp.tsv")

dir.create(FIGURES_P1, showWarnings = FALSE, recursive = TRUE)

# --- Load packages ---
for (pkg in c("ggplot2", "viridis")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(ggplot2)
library(viridis)

# --- Load BLAST results ---
df <- read.table(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df$species <- gsub(" ", "_", sub(".*\\|s__", "", gsub("\\[|\\]", "", df$taxon_full)))

# --- Select reference per species: lowest pident among same-species hits ---
# (tiebreak by longest alignment — ensures variant detection over perfect reference matching)
select_sp <- data.frame()
for (q in unique(df$query)) {
    d <- df[df$query == q & df$species == q, ]
    if (nrow(d) == 0) next
    d <- d[order(d$pident, -d$length), ][1, ]
    select_sp <- rbind(select_sp, d)
}
write.table(select_sp, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Reference table saved:", outfile, "\n")

# --- Prepare data for plots ---
get_family <- function(tx) {
    f <- sub(".*\\|f__", "", tx)
    f <- sub("\\|.*", "", f)
    f <- gsub("\\[|\\]| ", "", f)
    ifelse(f == "", "NA", f)
}

select_sp$pident <- as.numeric(select_sp$pident)
select_sp$qcovs  <- as.numeric(select_sp$qcovs)
select_sp$length <- as.numeric(select_sp$length)
select_sp$family <- get_family(select_sp$taxon_full)

# --- Plot 1: % identity per species ordered by pident (descending) ---
sp_ordered         <- select_sp[order(-select_sp$pident), ]
sp_ordered$query_f <- factor(sp_ordered$query, levels = sp_ordered$query)

p1 <- ggplot(sp_ordered, aes(x = query_f, y = pident, fill = pident)) +
    geom_col() +
    coord_cartesian(ylim = c(95, 100)) +
    scale_fill_viridis_c(option = "viridis", name = "% identidad") +
    labs(
        title    = "Porcentaje de identidad por especie",
        subtitle = "Ordenado de mayor a menor identidad",
        x        = "Especie",
        y        = "% identidad"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title         = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(size = 14, hjust = 0.5),
        axis.title.x       = element_text(size = 16, face = "bold"),
        axis.title.y       = element_text(size = 16, face = "bold"),
	axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        panel.grid.major.x = element_blank()
    )

ggsave(file.path(FIGURES_P1, "plot_sp_ordered_bypid.png"), p1, width = 15, height = 10, dpi = 200)
cat("Plot 1 saved\n")

# --- Plot 2: % identity per species grouped by family ---
fam_levels  <- sort(unique(sp_ordered$family))
fam_ordered <- select_sp[order(factor(select_sp$family, levels = fam_levels),
                               -select_sp$pident, select_sp$query), ]
fam_ordered$query_f  <- factor(fam_ordered$query,  levels = fam_ordered$query)
fam_ordered$family_f <- factor(fam_ordered$family, levels = fam_levels)

p2 <- ggplot(fam_ordered, aes(x = query_f, y = pident, fill = family_f)) +
    geom_col(color = "grey40", linewidth = 0.2) +
    coord_cartesian(ylim = c(95, 100)) +
    scale_fill_viridis_d(option = "turbo", name = "Familia") +
    labs(
        title = "Porcentaje de identidad por especie agrupado por familia",
        x     = "Especie",
        y     = "% identidad"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title         = element_text(size = 22, face = "bold", hjust = 0.5),
        legend.title       = element_text(size = 14, face = "bold"),
        axis.title.x       = element_text(size = 20, face = "bold"),
        axis.title.y       = element_text(size = 17, face = "bold"),
        axis.text.x        = element_text(size = 17, angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
        axis.text.y        = element_text(size = 20),
        panel.grid.major.x = element_blank()
    )

ggsave(file.path(FIGURES_P1, "plot_sp_by_family.png"), p2, width = 20, height = 11, dpi = 200)
cat("Plot 2 saved\n")
