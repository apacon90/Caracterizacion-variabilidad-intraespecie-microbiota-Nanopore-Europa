# 18_rrn_density.r — Variant density (VKb) per ribosomal region, species and country
# Input:  RESULTS_DIR/species_pool3_AF_val.csv  (validated variants, 16_export_variants.sh + manual curation)
#         RESULTS_DIR/master_rrn_coords_completo.tsv  (from 17_rrn_coords_full.sh)
# Output: RESULTS_DIR/p_box_VKb_validated.png
#         RESULTS_DIR/p_box_VKb_all.png
#         RESULTS_DIR/p_heat_VKb_validated.png
#         RESULTS_DIR/p_heat_VKb_all.png
#
# Usage: Rscript 18_rrn_density.r config.sh

# --- Parse config.sh ---
args     <- commandArgs(trailingOnly = TRUE)
cfg_file <- ifelse(length(args) > 0, args[1],
                   file.path(dirname(commandArgs()[4]), "../config.sh"))
cfg_lines <- readLines(cfg_file)
cfg_lines <- cfg_lines[grepl("^[A-Z_]+=", cfg_lines)]
for (line in cfg_lines) {
    parts <- strsplit(line, "=")[[1]]
    assign(parts[1], gsub('"', '', paste(parts[-1], collapse = "=")))
}

# --- Paths from config ---
variants_file <- file.path(RESULTS_DIR, "species_pool3_AF_val.csv")
coords_file   <- file.path(RESULTS_DIR, "master_rrn_coords_completo.tsv")
out_dir       <- RESULTS_DIR

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load packages ---
for (pkg in c("dplyr", "tidyr", "stringr", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# --- Load data ---
variants_raw <- read.table(variants_file, header = TRUE, sep = "\t")
coords_raw   <- read.table(coords_file,   header = TRUE, sep = "\t")

# --- Validation summary ---
cat("Total variantes candidatas:", nrow(variants_raw), "\n")
cat("Validadas en core_nt (SI):", sum(variants_raw$VALIDACION_CORE == "SI"), "\n")
cat("Validadas en core_nt o WGS:", sum(variants_raw$VALIDACION_CORE == "SI" | variants_raw$VALIDACION_WGS == "SI"), "\n")
cat("core_nt sin referencia, validadas en WGS:", sum(variants_raw$VALIDACION_CORE == "CORE_NO_REF" & variants_raw$VALIDACION_WGS == "SI"), "\n")
cat("core_nt sin referencia, no validadas en WGS:", sum(variants_raw$VALIDACION_CORE == "CORE_NO_REF" & variants_raw$VALIDACION_WGS == "NO"), "\n")
cat("No validadas en ninguna BD:", sum(variants_raw$VALIDACION_CORE == "NO" & variants_raw$VALIDACION_WGS == "NO"), "\n\n")

# --- Prepare coords (deduplicate) ---
coords <- coords_raw %>%
    distinct(Species, Contig, Region, Start, End, .keep_all = TRUE)

# --- Prepare variants ---
variants <- variants_raw %>%
    rename(Contig = CHROM) %>%
    mutate(VALIDADA = (VALIDACION_CORE == "SI") | (VALIDACION_WGS == "SI"))

cat("Tabla VALIDACION_CORE:\n"); print(table(variants$VALIDACION_CORE))
cat("Tabla VALIDACION_WGS:\n");  print(table(variants$VALIDACION_WGS))
cat("Tabla VALIDADA:\n");        print(table(variants$VALIDADA))

# --- Convert to long format ---
var_long <- variants %>%
    pivot_longer(
        cols      = starts_with("AF_"),
        names_to  = "Country",
        values_to = "AF"
    ) %>%
    mutate(Country = str_remove(Country, "^AF_"))

# --- Assign region to each variant ---
var_region <- var_long %>%
    left_join(coords, by = c("Species", "Contig"),
              relationship = "many-to-many") %>%
    filter(POS >= Start & POS <= End) %>%
    filter(Region != "5S_rRNA")

var_region_all <- var_region
var_region_val <- subset(var_region, VALIDADA)

# --- Density function (VKb) ---
calc_density <- function(vr, coords_df) {

    all_species   <- unique(coords_df$Species)
    all_regions   <- setdiff(unique(coords_df$Region), "5S_rRNA")
    all_countries <- unique(vr$Country)

    grid <- expand_grid(
        Species = all_species,
        Country = all_countries,
        Region  = all_regions
    )

    counts <- vr %>%
        filter(AF > 0, Region != "5S_rRNA") %>%
        distinct(Species, Country, Region, Contig, POS, REF, ALT, TYPE) %>%
        count(Species, Country, Region, name = "Nvar")

    lengths <- coords_df %>%
        filter(Region != "5S_rRNA") %>%
        mutate(kb = (End - Start + 1) / 1000) %>%
        group_by(Species, Region) %>%
        summarise(kb = sum(kb), .groups = "drop")

    grid %>%
        left_join(counts,  by = c("Species", "Country", "Region")) %>%
        left_join(lengths, by = c("Species", "Region")) %>%
        mutate(
            Nvar = replace_na(Nvar, 0),
            VKb  = Nvar / kb
        ) %>%
        filter(Region != "5S_rRNA")
}

density_all <- calc_density(var_region_all, coords)
density_val <- calc_density(var_region_val, coords)

stopifnot(all(is.finite(density_all$VKb)))
stopifnot(all(is.finite(density_val$VKb)))

# --- Plot functions ---
region_order <- c("16S_rRNA", "ITS", "23S_rRNA")

plot_box <- function(density_df, main_title) {
    ggplot(density_df, aes(x = factor(Region, levels = region_order), y = VKb)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(aes(color = Species),
                    width = 0.15, alpha = 0.5, size = 1.2, show.legend = FALSE) +
        facet_wrap(~ Country, nrow = 1) +
        labs(title = main_title, x = "Región funcional", y = "VKb (variantes por kb)") +
        theme_bw(base_size = 13) +
        theme(
            strip.background = element_rect(fill = "grey90", color = "black"),
            strip.text       = element_text(face = "bold"),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            panel.grid       = element_blank()
        )
}

plot_heat <- function(density_df, main_title) {
    density_df$Region <- factor(density_df$Region, levels = region_order)
    ggplot(density_df, aes(x = Region, y = Species, fill = VKb)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "blue") +
        facet_grid(. ~ Country) +
        labs(
            title    = main_title,
            subtitle = "Densidad de variantes (VKb), separado por país",
            x        = "Región del rrn", y = "Especie", fill = "VKb"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x      = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "grey93"),
            panel.grid       = element_blank(),
            axis.title.y     = element_blank(),
            axis.text.y      = element_text(size = 6),
            panel.spacing    = unit(0.3, "lines")
        )
}

# --- Generate and save plots ---
p_box_val  <- plot_box(density_val,  "Densidad VKb por región (solo variantes validadas)")
p_heat_val <- plot_heat(density_val, "Hotspots de variación (solo variantes validadas)")
p_box_all  <- plot_box(density_all,  "Densidad VKb por región (todas las variantes)")
p_heat_all <- plot_heat(density_all, "Hotspots de variación (todas las variantes)")

ggsave(file.path(out_dir, "p_box_VKb_validated.png"),  p_box_val,  width = 10, height = 4, dpi = 300)
ggsave(file.path(out_dir, "p_heat_VKb_validated.png"), p_heat_val, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "p_box_VKb_all.png"),        p_box_all,  width = 10, height = 4, dpi = 300)
ggsave(file.path(out_dir, "p_heat_VKb_all.png"),       p_heat_all, width = 12, height = 7, dpi = 300)

cat("Done. Plots saved in:", out_dir, "\n")
