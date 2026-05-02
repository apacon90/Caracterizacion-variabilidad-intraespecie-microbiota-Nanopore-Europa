# 03_variant_level_analysis.R — Análisis de variantes por país
# Input:  RDS de 01_prepare_data.R
#         phase2_data/variants_tsv/*.FULL.tsv
# Output: tables/tabla_final_variantes.csv
#         figures/phase2/heatmap_AF_top20.png

set.seed(42)
library(tidyverse)
library(pheatmap)
library(viridis)

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
meta_qc   <- readRDS(file.path(rds_dir, "meta_qc.rds"))
mega_main <- readRDS(file.path(rds_dir, "mega_main.rds"))
Xa_main   <- readRDS(file.path(rds_dir, "Xa_main.rds"))

############################################################
## 2. EXTRAER AF DESDE CAMPO AD
############################################################
extract_af <- function(f, ids, min_dp = 5) {
  df  <- read_tsv(f,
                  col_types = cols(CHROM = col_character(), POS = col_integer(),
                                   REF = col_character(), ALT = col_character(),
                                   .default = col_character()),
                  show_col_types = FALSE)
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
  
  X <- df %>% select(all_of(smp)) %>% mutate(across(everything(), parse_af))
  for (id in setdiff(ids, smp)) X[[id]] <- NA_real_
  X <- X[, ids, drop = FALSE]
  bind_cols(
    df %>% transmute(Species = sp, CHROM, POS, REF, ALT,
                     VariantID = paste(sp, CHROM, POS, REF, ALT, sep = "|")), X
  )
}

full_files <- list.files(dir_full, pattern = "FULL\\.tsv$", full.names = TRUE)
cat("Ficheros FULL encontrados:", length(full_files), "\n")

mega_full_main <- bind_rows(lapply(full_files, extract_af, ids = meta_qc$sample)) %>%
  filter(!grepl("^Ralstonia", Species)) %>%
  filter(VariantID %in% rownames(Xa_main))

cat("Variantes en mega_full_main:", nrow(mega_full_main), "\n")

############################################################
## 3. AF MEDIA POR PAÍS
############################################################
af_real <- mega_full_main %>%
  select(Species, VariantID, all_of(meta_qc$sample)) %>%
  pivot_longer(cols = all_of(meta_qc$sample), names_to = "sample", values_to = "AF_ind") %>%
  filter(!is.na(AF_ind)) %>%
  left_join(meta_qc %>% select(sample, Country), by = "sample") %>%
  group_by(VariantID, Species, Country) %>%
  summarise(AF = mean(AF_ind), n_obs = n(), .groups = "drop") %>%
  pivot_wider(names_from = Country, values_from = c(AF, n_obs)) %>%
  mutate(
    max_diff = pmax(abs(AF_Norway - AF_Portugal), abs(AF_Norway - AF_Spain),
                    abs(AF_Portugal - AF_Spain), na.rm = TRUE),
    max_diff = ifelse(is.finite(max_diff), max_diff, NA_real_)
  ) %>%
  arrange(desc(max_diff))

############################################################
## 4. FISHER POR VARIANTE
############################################################
fisher_res <- mega_main %>%
  select(Species, VariantID, all_of(meta_qc$sample)) %>%
  pivot_longer(cols = all_of(meta_qc$sample), names_to = "sample", values_to = "GT") %>%
  filter(!is.na(GT)) %>%
  left_join(meta_qc %>% select(sample, Country), by = "sample") %>%
  group_by(VariantID, Species) %>%
  summarise(
    n_Norway   = sum(Country == "Norway"),
    n_Portugal = sum(Country == "Portugal"),
    n_Spain    = sum(Country == "Spain"),
    p_fisher   = tryCatch(
      fisher.test(table(GT, Country), simulate.p.value = TRUE, B = 10000)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p_fisher, method = "fdr")) %>%
  arrange(fdr)

cat("Variantes FDR<0.05:", sum(fisher_res$fdr < 0.05, na.rm = TRUE),
    "| FDR<0.10:", sum(fisher_res$fdr < 0.10, na.rm = TRUE), "\n")

############################################################
## 5. TABLA FINAL
############################################################
tabla_final <- af_real %>%
  left_join(fisher_res %>% select(VariantID, n_Norway, n_Portugal, n_Spain, p_fisher, fdr),
            by = "VariantID") %>%
  arrange(fdr) %>%
  mutate(
    enriched_in = case_when(
      is.na(fdr) | fdr >= 0.05                                    ~ "ns",
      is.na(AF_Norway) | is.na(AF_Portugal) | is.na(AF_Spain)     ~ "nd",
      AF_Norway   >= AF_Portugal & AF_Norway   >= AF_Spain         ~ "Norway",
      AF_Portugal >= AF_Norway   & AF_Portugal >= AF_Spain         ~ "Portugal",
      TRUE                                                         ~ "Spain"
    )
  ) %>%
  select(Species, VariantID, AF_Norway, AF_Portugal, AF_Spain,
         n_obs_Norway, n_obs_Portugal, n_obs_Spain, max_diff, p_fisher, fdr, enriched_in)

write_csv(tabla_final, file.path(tab_dir, "tabla_final_variantes.csv"))

cat("Variantes enriquecidas por país (FDR<0.05):\n")
tabla_final %>% filter(fdr < 0.05) %>%
  count(enriched_in) %>% print()

############################################################
## 6. HEATMAP AF — TOP 20
############################################################
top20 <- af_real %>%
  filter(!is.na(AF_Norway), !is.na(AF_Portugal), !is.na(AF_Spain)) %>%
  slice_head(n = 20) %>%
  mutate(label = paste0(word(Species, 1), "_",
                        str_extract(VariantID, "(?<=\\|)\\d+(?=\\|)")))

mat_heat <- top20 %>%
  select(label, AF_Norway, AF_Portugal, AF_Spain) %>%
  column_to_rownames("label") %>%
  as.matrix()

colnames(mat_heat) <- c("Noruega", "Portugal", "España")

ann_row <- top20 %>% select(label, Species) %>% column_to_rownames("label")

row_labs <- as.expression(
  lapply(rownames(mat_heat), function(x) bquote(italic(.(x))))
)

pheatmap(
  mat_heat,
  labels_row       = row_labs,
  cluster_cols     = FALSE,
  cluster_rows     = TRUE,
  annotation_row   = ann_row,
  annotation_legend     = FALSE,
  annotation_names_row  = FALSE,
  color            = viridis(101),
  breaks           = seq(0, 1, length.out = 101),
  main             = "Frecuencia alélica media por país\nTop 20 variantes",
  fontsize_row     = 8,
  angle_col        = 45,
  border_color     = NA,
  cellwidth        = 35,
  cellheight       = 18,
  filename         = file.path(fig_dir, "heatmap_AF_top20.png"),
  width            = 6.5,
  height           = 7
)

cat("Heatmap saved\n")

############################################################
## 7. GUARDAR RDS
############################################################
saveRDS(af_real,     file.path(rds_dir, "af_real.rds"))
saveRDS(fisher_res,  file.path(rds_dir, "fisher_res.rds"))
saveRDS(tabla_final, file.path(rds_dir, "tabla_final.rds"))

cat("Done. Outputs guardados en:", out_dir, "\n")

