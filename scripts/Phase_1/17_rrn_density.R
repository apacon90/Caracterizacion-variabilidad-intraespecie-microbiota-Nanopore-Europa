# 17_rrn_density.R — Variant density (VKb) per ribosomal region, species and country
# Input:  RESULTS_DIR/species_pool3_AF_val.csv  (validated variants, 15_export_variants.sh + manual curation)
#         RESULTS_DIR/master_rrn_coords_completo.tsv  (from 16_rrn_coords_full.sh)
# Output: RESULTS_DIR/figures/phase1/p_box_VKb_validated.png
#         RESULTS_DIR/figures/phase1/p_box_VKb_all.png
#         RESULTS_DIR/figures/phase1/p_heat_VKb_validated.png
#         RESULTS_DIR/figures/phase1/p_heat_VKb_all.png
#         RESULTS_DIR/tables/vkb_descriptive_by_region.tsv
#         RESULTS_DIR/tables/kruskal_vkb_region.txt
#         RESULTS_DIR/tables/kruskal_vkb_country_by_region.txt
#
# Usage: Rscript 17_rrn_density.R config.sh

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

# --- Build derived paths ---
TABLES_DIR  <- file.path(RESULTS_DIR, "tables")
FIGURES_P1  <- file.path(RESULTS_DIR, "figures", "phase1")
FIGURES_P2  <- file.path(RESULTS_DIR, "figures", "phase2")

# --- Paths from config ---
variants_file <- file.path(RESULTS_DIR, "species_pool3_AF_val.csv")
coords_file   <- file.path(RESULTS_DIR, "master_rrn_coords_completo.tsv")
out_dir       <- FIGURES_P1

dir.create(out_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load packages ---
for (pkg in c("dplyr", "tidyr", "stringr", "ggplot2", "viridis",
              "tibble", "rstatix")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(tibble)
library(rstatix)   # wilcox_test pairwise post-hoc

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

# --- Translate country names to Spanish ---
var_long <- var_long %>%
    mutate(Country = recode(Country,
        "Norway"   = "Noruega",
        "Portugal" = "Portugal",
        "Spain"    = "España"
    ))

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

# =============================================================================
# --- ANALYSIS 1: Descriptive summary — median and IQR by region ---
# =============================================================================

cat("\nTotal VKb observations (zeros included):", nrow(density_val), "\n")

desc_region <- density_val %>%
    group_by(Region) %>%
    summarise(
        n        = n(),
        mediana  = round(median(VKb), 3),
        Q1       = round(quantile(VKb, 0.25), 3),
        Q3       = round(quantile(VKb, 0.75), 3),
        IQR      = round(IQR(VKb), 3),
        .groups  = "drop"
    ) %>%
    mutate(mediana_IQR = paste0(mediana, " (", Q1, "–", Q3, ")"))

cat("\n--- Descriptive VKb by region (median and IQR) ---\n")
print(desc_region)

write.table(
    desc_region,
    file      = file.path(TABLES_DIR, "vkb_descriptive_by_region.tsv"),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
)

# Also by region and country
desc_region_country <- density_val %>%
    group_by(Region, Country) %>%
    summarise(
        n       = n(),
        mediana = round(median(VKb), 3),
        Q1      = round(quantile(VKb, 0.25), 3),
        Q3      = round(quantile(VKb, 0.75), 3),
        IQR     = round(IQR(VKb), 3),
        .groups = "drop"
    ) %>%
    mutate(mediana_IQR = paste0(mediana, " (", Q1, "–", Q3, ")"))

cat("\n--- Descriptive VKb by region and country (median and IQR) ---\n")
print(desc_region_country)

write.table(
    desc_region_country,
    file      = file.path(TABLES_DIR, "vkb_descriptive_by_region_country.tsv"),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
)

# =============================================================================
# --- ANALYSIS 2: Kruskal-Wallis + pairwise Wilcoxon (Bonferroni) ---
#
# Zeros are kept: they reflect real absence of variants and match the boxplot.
# =============================================================================

# ── 2a. Does VKb differ between 16S, ITS and 23S? ────────────────────────────

cat("\n--- Kruskal-Wallis: VKb by region (zeros included) ---\n")
kw_region <- kruskal.test(VKb ~ Region, data = density_val)
print(kw_region)

cat("\n--- Pairwise Wilcoxon (Bonferroni): VKb by region ---\n")
pw_region <- density_val %>%
    wilcox_test(VKb ~ Region, p.adjust.method = "bonferroni")
print(pw_region)

capture.output(
    {
        cat("=== Descriptive VKb by region (median, IQR) ===\n")
        print(as.data.frame(desc_region))
        cat("\n=== Kruskal-Wallis: VKb by region (zeros included) ===\n")
        print(kw_region)
        cat("\n=== Pairwise Wilcoxon (Bonferroni) ===\n")
        print(as.data.frame(pw_region))
    },
    file = file.path(TABLES_DIR, "kruskal_vkb_region.txt")
)

# ── 2b. Does VKb differ between countries within each region? ────────────────

cat("\n--- Kruskal-Wallis: VKb by country within each region (zeros included) ---\n")

kw_country_list <- list()
pw_country_list <- list()

for (reg in c("16S_rRNA", "ITS", "23S_rRNA")) {
    cat("\nRegion:", reg, "\n")

    df_reg <- density_val %>% filter(Region == reg)

    kw <- kruskal.test(VKb ~ Country, data = df_reg)
    print(kw)
    kw_country_list[[reg]] <- kw

    if (kw$p.value < 0.05) {
        cat("  Pairwise Wilcoxon (Bonferroni):\n")
        pw <- df_reg %>%
            wilcox_test(VKb ~ Country, p.adjust.method = "bonferroni")
        print(pw)
        pw_country_list[[reg]] <- pw
    } else {
        cat("  Not significant — no post-hoc performed\n")
    }
}

capture.output(
    {
        for (r in names(kw_country_list)) {
            cat("\n=== Region:", r, "===\n")
            cat("--- Kruskal-Wallis ---\n")
            print(kw_country_list[[r]])
            if (!is.null(pw_country_list[[r]])) {
                cat("--- Pairwise Wilcoxon (Bonferroni) ---\n")
                print(as.data.frame(pw_country_list[[r]]))
            }
        }
    },
    file = file.path(TABLES_DIR, "kruskal_vkb_country_by_region.txt")
)

# --- Translate country names in density tables ---
translate_country <- function(df) {
    df %>% mutate(Country = recode(Country,
        "Norway"   = "Noruega",
        "Portugal" = "Portugal",
        "Spain"    = "España"
    ))
}
density_all <- translate_country(density_all)
density_val <- translate_country(density_val)

# --- Plot functions ---
region_order  <- c("16S_rRNA", "ITS", "23S_rRNA")
country_order <- c("Noruega", "Portugal", "España")

plot_box <- function(density_df, main_title) {
    density_df$Country <- factor(density_df$Country, levels = country_order)
    ggplot(density_df, aes(x = factor(Region, levels = region_order), y = VKb)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "grey90") +
        geom_jitter(width = 0.15, alpha = 0.6, size = 1.5, color = "#0072B2") +
        facet_wrap(~ Country, nrow = 1) +
        labs(
            title = main_title,
            x     = "Región funcional",
            y     = "VKb (variantes por kb)"
        ) +
        theme_bw(base_size = 13) +
        theme(
            strip.background = element_rect(fill = "grey90", color = "black"),
            strip.text       = element_text(face = "bold"),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            panel.grid       = element_blank()
        )
}

plot_heat <- function(density_df, main_title) {
    density_df$Region  <- factor(density_df$Region,  levels = region_order)
    density_df$Country <- factor(density_df$Country, levels = country_order)
    ggplot(density_df, aes(x = Region, y = Species, fill = VKb)) +
        geom_tile() +
        scale_fill_viridis_c(option = "viridis", name = "VKb") +
        facet_grid(. ~ Country) +
        labs(
            title    = main_title,
            subtitle = "Densidad de variantes (VKb) por región funcional y país",
            x        = "Región del rrn",
            y        = "Especie"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x      = element_text(angle = 45, hjust = 1),
            axis.title.y     = element_text(face = "bold"),
            axis.text.y      = element_text(size = 8, face = "italic"),
            strip.background = element_rect(fill = "grey93", color = "black"),
            strip.text       = element_text(size = 14, face = "bold"),
            plot.title       = element_text(size = 18, face = "bold"),
            plot.subtitle    = element_text(size = 14),
            panel.grid       = element_blank(),
            panel.spacing    = unit(0.3, "lines")
        )
}

cat("\n--- Tipo de variante en validadas (n=418) ---\n")
variants %>%
    filter(VALIDADA == TRUE) %>%
    count(TYPE) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    print()

# --- Generate and save plots ---
p_box_val  <- plot_box(density_val,  "Densidad VKb por región (variantes validadas)")
p_heat_val <- plot_heat(density_val, "Puntos calientes de variación (variantes validadas)")
p_box_all  <- plot_box(density_all,  "Densidad VKb por región (todas las variantes candidatas)")
p_heat_all <- plot_heat(density_all, "Puntos calientes de variación (todas las variantes candidatas)")

ggsave(file.path(out_dir, "p_box_VKb_validated.png"),  p_box_val,  width = 10, height = 4, dpi = 300)
ggsave(file.path(out_dir, "p_heat_VKb_validated.png"), p_heat_val, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "p_box_VKb_all.png"),        p_box_all,  width = 10, height = 4, dpi = 300)
ggsave(file.path(out_dir, "p_heat_VKb_all.png"),       p_heat_all, width = 12, height = 7, dpi = 300)

cat("Done. Plots saved in:", out_dir, "\n")
