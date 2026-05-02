# 05_pca_af_exploratorio.R — PCA exploratorio sobre frecuencias alélicas
# Input:  RDS de 01_prepare_data.R y 03_variant_level_analysis.R
#         phase2_data/variants_tsv/*.FULL.tsv
# Output: figures/phase2/PCA_AF_scores.png
#         figures/phase2/PCA_AF_biplot_top_variantes.png
#         tables/PCA_AF_*.csv

library(tidyverse)
library(ggrepel)

############################################################
## 0. RUTAS
############################################################
base_dir <- "/media/apacon/SSD_APACON/1_SCRIPTS_FINALES_TFM_run_final/GROND_nr_Rep"

dir_full <- file.path(base_dir, "phase2_data", "variants_tsv")
out_dir  <- file.path(base_dir, "metadata")
fig_dir  <- file.path(out_dir, "figures", "phase2")
tab_dir  <- file.path(out_dir, "tables")
rds_dir  <- file.path(out_dir, "rds")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. CARGAR RDS
############################################################
meta_model  <- readRDS(file.path(rds_dir, "meta_model.rds"))
Xa_main     <- readRDS(file.path(rds_dir, "Xa_main.rds"))

# Colores / formas (paleta Wong, dalton-safe)
COL <- c("Norway" = "#E69F00", "Portugal" = "#CC79A7", "Spain" = "#0072B2")
SHP <- c("Man" = 16, "Woman" = 17)

COUNTRY_ES <- c("Norway" = "Noruega", "Portugal" = "Portugal", "Spain" = "España")
SEX_ES     <- c("Man" = "Hombre", "Woman" = "Mujer")

############################################################
## 2. PARÁMETROS
############################################################
min_dp          <- 5     # profundidad mínima para calcular AF
min_prop_global <- 0.60  # proporción mínima de individuos con AF observada
min_n_country   <- 3     # observaciones mínimas por país
top_n_loadings  <- 15    # top variantes a etiquetar en biplot

############################################################
## 3. FUNCIÓN: extraer AF desde FULL.tsv
############################################################
extract_af <- function(f, ids, min_dp = 5) {
  df <- read_tsv(
    f,
    col_types = cols(
      CHROM = col_character(), POS = col_integer(),
      REF   = col_character(), ALT = col_character(),
      .default = col_character()
    ),
    show_col_types = FALSE
  )
  
  sp  <- sub("\\.QUAL_.*$", "", basename(f))
  smp <- names(df)[6:ncol(df)]
  
  parse_af <- function(x) {
    sapply(x, function(v) {
      if (is.na(v) || v == "." || v == "") return(NA_real_)
      fields <- str_split(v, ":")[[1]]
      if (length(fields) < 6 || is.na(fields[6]) || fields[6] == ".") return(NA_real_)
      counts <- suppressWarnings(as.numeric(str_split(fields[6], ",")[[1]]))
      if (length(counts) < 2 || any(is.na(counts)) || sum(counts) < min_dp) return(NA_real_)
      counts[2] / sum(counts)
    }, USE.NAMES = FALSE)
  }
  
  X <- df %>%
    select(all_of(smp)) %>%
    mutate(across(everything(), parse_af))
  
  for (id in setdiff(ids, smp)) X[[id]] <- NA_real_
  X <- X[, ids, drop = FALSE]
  
  bind_cols(
    df %>% transmute(
      Species   = sp,
      CHROM     = CHROM, POS = POS, REF = REF, ALT = ALT,
      VariantID = paste(sp, CHROM, POS, REF, ALT, sep = "|")
    ),
    X
  )
}

############################################################
## 4. LEER AF Y FILTRAR AL UNIVERSO DE Xa_main
############################################################
full_files <- list.files(dir_full, pattern = "\\.FULL\\.tsv$", full.names = TRUE)
cat("Ficheros FULL encontrados:", length(full_files), "\n")

mega_af <- bind_rows(lapply(full_files, extract_af,
                            ids = meta_model$sample, min_dp = min_dp)) %>%
  filter(!grepl("^Ralstonia", Species)) %>%
  filter(VariantID %in% rownames(Xa_main))

cat("Variantes candidatas para PCA (antes de filtros extra):", nrow(mega_af), "\n")

############################################################
## 5. FORMATO LONG Y COBERTURA POR VARIANTE
############################################################
af_long <- mega_af %>%
  select(Species, VariantID, all_of(meta_model$sample)) %>%
  pivot_longer(
    cols      = all_of(meta_model$sample),
    names_to  = "sample",
    values_to = "AF"
  ) %>%
  left_join(meta_model %>% select(sample, Country, Sex, Age), by = "sample")

n_total_samples <- nrow(meta_model)

variant_cov <- af_long %>%
  group_by(VariantID, Species) %>%
  summarise(
    n_obs_global = sum(!is.na(AF)),
    prop_global  = n_obs_global / n_total_samples,
    n_Norway     = sum(!is.na(AF) & Country == "Norway"),
    n_Portugal   = sum(!is.na(AF) & Country == "Portugal"),
    n_Spain      = sum(!is.na(AF) & Country == "Spain"),
    mean_AF      = mean(AF, na.rm = TRUE),
    sd_AF        = sd(AF,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(sd_AF = replace_na(sd_AF, 0))

write_csv(variant_cov, file.path(tab_dir, "PCA_AF_variant_coverage_summary.csv"))

############################################################
## 6. FILTRO ESTRICTO PARA PCA
############################################################
keep_variants <- variant_cov %>%
  filter(
    prop_global >= min_prop_global,
    n_Norway    >= min_n_country,
    n_Portugal  >= min_n_country,
    n_Spain     >= min_n_country,
    sd_AF > 0
  )

cat("Variantes retenidas para PCA:", nrow(keep_variants), "\n")

if (nrow(keep_variants) < 3)
  stop("Muy pocas variantes tras el filtrado. Relaja min_prop_global o min_n_country.")

############################################################
## 7. MATRIZ MUESTRA × VARIANTE
############################################################
af_wide <- af_long %>%
  filter(VariantID %in% keep_variants$VariantID) %>%
  select(sample, VariantID, AF) %>%
  pivot_wider(names_from = VariantID, values_from = AF) %>%
  right_join(meta_model %>% select(sample, Country, Sex, Age), by = "sample") %>%
  relocate(sample, Country, Sex, Age)

mat <- af_wide %>%
  select(-sample, -Country, -Sex, -Age) %>%
  as.data.frame()
rownames(mat) <- af_wide$sample

############################################################
## 8. MISSINGNESS
############################################################
na_by_var <- colMeans(is.na(mat))
na_by_ind <- rowMeans(is.na(mat))

write_csv(
  tibble(VariantID = names(na_by_var), prop_na = as.numeric(na_by_var)),
  file.path(tab_dir, "PCA_AF_missingness_por_variante.csv")
)
write_csv(
  tibble(sample = rownames(mat), prop_na = na_by_ind) %>%
    left_join(meta_model %>% select(sample, Country), by = "sample"),
  file.path(tab_dir, "PCA_AF_missingness_por_individuo.csv")
)

cat("Missingness medio por variante:", round(mean(na_by_var), 3), "\n")
cat("Missingness medio por individuo:", round(mean(na_by_ind), 3), "\n")

############################################################
## 9. IMPUTACIÓN SIMPLE (media por variante)
############################################################
mat_imp <- mat

for (j in seq_len(ncol(mat_imp))) {
  v <- mat_imp[[j]]
  m <- mean(v, na.rm = TRUE)
  if (is.nan(m)) m <- 0
  v[is.na(v)] <- m
  mat_imp[[j]] <- v
}

# Eliminar columnas con varianza cero tras imputación
sds     <- sapply(mat_imp, sd, na.rm = TRUE)
mat_imp <- mat_imp[, sds > 0, drop = FALSE]

cat("Dimensión final matriz PCA:", nrow(mat_imp), "muestras ×", ncol(mat_imp), "variantes\n")

if (ncol(mat_imp) < 2)
  stop("Muy pocas variantes con varianza > 0 para PCA.")

############################################################
## 10. PCA
############################################################
pca     <- prcomp(mat_imp, center = TRUE, scale. = TRUE)
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

tabla_var_exp <- tibble(
  PC                 = paste0("PC", seq_along(var_exp)),
  variance_explained = var_exp,
  percent            = round(var_exp * 100, 2),
  cumulative_percent = round(cumsum(var_exp) * 100, 2)
)
write_csv(tabla_var_exp, file.path(tab_dir, "PCA_AF_varianza_explicada.csv"))

############################################################
## 11. SCORES DE MUESTRAS
############################################################
scores <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(meta_model, by = "sample") %>%
  mutate(
    Country_es = recode(Country, !!!COUNTRY_ES),
    Sex_es     = recode(Sex,     !!!SEX_ES)
  )

pct1 <- round(var_exp[1] * 100, 1)
pct2 <- round(var_exp[2] * 100, 1)
cat("PC1:", pct1, "% | PC2:", pct2, "%\n")

############################################################
## 12. LOADINGS
############################################################
loadings_df <- as.data.frame(pca$rotation[, 1:2]) %>%
  rownames_to_column("VariantID") %>%
  left_join(keep_variants %>% select(VariantID, Species, prop_global,
                                     n_Norway, n_Portugal, n_Spain),
            by = "VariantID") %>%
  mutate(
    contrib = sqrt(PC1^2 + PC2^2),
    label   = paste0(
      word(Species, 1), "_",
      sapply(str_split(VariantID, "\\|"), `[`, 3)
    )
  ) %>%
  arrange(desc(contrib))

write_csv(loadings_df, file.path(tab_dir, "PCA_AF_loadings_PC1_PC2.csv"))
top_load <- loadings_df %>% slice_head(n = top_n_loadings)

############################################################
## 13. FIGURAS
############################################################
n_country      <- scores %>% count(Country_es) %>% deframe()
country_labels <- setNames(paste0(names(n_country), " (n=", n_country, ")"),
                           names(n_country))

COL_ES <- setNames(COL, COUNTRY_ES[names(COL)])
SHP_ES <- setNames(SHP, SEX_ES[names(SHP)])

# --- Scores ---
p_scores <- ggplot(scores, aes(PC1, PC2, color = Country_es, shape = Sex_es)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Country_es), level = 0.95, linetype = 2) +
  scale_color_manual(values = COL_ES, labels = country_labels) +
  scale_shape_manual(values = SHP_ES) +
  theme_bw() +
  labs(
    title    = "PCA exploratorio con frecuencia alélica",
    subtitle = paste0("Variantes filtradas: n=", ncol(mat_imp),
                      " | Imputación por media de variante"),
    x        = paste0("PC1 (", pct1, "%)"),
    y        = paste0("PC2 (", pct2, "%)"),
    color    = "País", shape = "Sexo"
  )

ggsave(file.path(fig_dir, "PCA_AF_scores.png"),
       p_scores, width = 8, height = 6, dpi = 300)
cat("PCA scores guardado\n")

# --- Biplot ---
arrow_mult <- 5

p_biplot <- ggplot() +
  geom_point(data = scores,
             aes(PC1, PC2, color = Country_es, shape = Sex_es),
             size = 2.5, alpha = 0.7) +
  stat_ellipse(data = scores,
               aes(PC1, PC2, color = Country_es, group = Country_es),
               level = 0.95, linetype = 2) +
  geom_segment(data = top_load,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_mult, yend = PC2 * arrow_mult),
               arrow = arrow(length = unit(0.22, "cm")),
               color = "black", linewidth = 0.7) +
  geom_text_repel(data = top_load,
                  aes(x = PC1 * arrow_mult, y = PC2 * arrow_mult, label = label),
                  size = 3.2, max.overlaps = 30) +
  scale_color_manual(values = COL_ES, labels = country_labels) +
  scale_shape_manual(values = SHP_ES) +
  theme_bw() +
  labs(
    title    = "Biplot PCA exploratorio con frecuencia alélica",
    subtitle = "Top variantes por magnitud de loading en PC1-PC2",
    x        = paste0("PC1 (", pct1, "%)"),
    y        = paste0("PC2 (", pct2, "%)"),
    color    = "País", shape = "Sexo"
  )

ggsave(file.path(fig_dir, "PCA_AF_biplot_top_variantes.png"),
       p_biplot, width = 10, height = 7, dpi = 300)
cat("Biplot guardado\n")

############################################################
## 14. TABLAS TOP CONTRIBUYENTES
############################################################
write_csv(loadings_df %>% arrange(desc(PC1)) %>% slice_head(n = 15),
          file.path(tab_dir, "PCA_AF_top_PC1_positivos.csv"))
write_csv(loadings_df %>% arrange(PC1)        %>% slice_head(n = 15),
          file.path(tab_dir, "PCA_AF_top_PC1_negativos.csv"))
write_csv(loadings_df %>% arrange(desc(PC2)) %>% slice_head(n = 15),
          file.path(tab_dir, "PCA_AF_top_PC2_positivos.csv"))
write_csv(loadings_df %>% arrange(PC2)        %>% slice_head(n = 15),
          file.path(tab_dir, "PCA_AF_top_PC2_negativos.csv"))

############################################################
## 15. GUARDAR RDS
############################################################
saveRDS(mega_af,       file.path(rds_dir, "pca_af_mega_af_raw.rds"))
saveRDS(keep_variants, file.path(rds_dir, "pca_af_keep_variants.rds"))
saveRDS(mat_imp,       file.path(rds_dir, "pca_af_matrix_imputed.rds"))
saveRDS(pca,           file.path(rds_dir, "pca_af_prcomp.rds"))
saveRDS(scores,        file.path(rds_dir, "pca_af_scores.rds"))
saveRDS(loadings_df,   file.path(rds_dir, "pca_af_loadings.rds"))

cat("Done. PCA AF exploratorio finalizado.\n")

