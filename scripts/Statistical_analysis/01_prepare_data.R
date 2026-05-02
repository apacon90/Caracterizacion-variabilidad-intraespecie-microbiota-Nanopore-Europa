# 01_prepare_data.R — Carga y preparación de datos
# Input:  dir_gt (ficheros .GT.tsv), sample_metadata.csv,
#         reads_by_species_matrix_unique.tsv,
#         validated_variants_blast.tsv (variantes validadas por BLASTn)
# Output: RDS en rds_dir, CSVs en tab_dir

############################################################
## 0. RUTAS
############################################################
base_dir   <- "/media/apacon/SSD_APACON/1_SCRIPTS_FINALES_TFM_run_final/GROND_nr_Rep"
setwd(base_dir)

dir_gt     <- file.path(base_dir, "phase2_data", "variants_tsv")
meta_csv   <- file.path(base_dir, "metadata", "sample_metadata.csv")
reads_file <- file.path(base_dir, "metadata", "tables", "reads_by_species_matrix_unique.tsv")
val_file   <- file.path(base_dir, "metadata", "tables", "validated_variants_blast.tsv")

out_dir <- file.path(base_dir, "metadata")
rds_dir <- file.path(out_dir, "rds")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures", "phase2")

dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)

############################################################
## 1. METADATA
############################################################
meta <- read_csv(meta_csv, show_col_types = FALSE) %>%
  rename(sample = Sample) %>%
  mutate(
    Sex     = na_if(Sex, "N/A") |> na_if("NA") |> na_if("") |> factor(),
    Age     = na_if(Age, "N/A") |> na_if("NA") |> na_if("") |> as.numeric(),
    Country = factor(Country)
  )

tabla_desc <- meta %>%
  group_by(Country) %>%
  summarise(
    n             = n(),
    n_hombres     = sum(Sex == "Man", na.rm = TRUE),
    n_mujeres     = sum(Sex == "Woman", na.rm = TRUE),
    n_sex_na      = sum(is.na(Sex)),
    n_age_na      = sum(is.na(Age)),
    n_any_missing = sum(is.na(Sex) | is.na(Age)),
    edad_media    = round(mean(Age, na.rm = TRUE), 1),
    edad_sd       = round(sd(Age, na.rm = TRUE), 1),
    edad_min      = min(Age, na.rm = TRUE),
    edad_max      = max(Age, na.rm = TRUE),
    .groups = "drop"
  )
print(tabla_desc)
write_csv(tabla_desc, file.path(tab_dir, "tabla_descriptiva_paises.csv"))

############################################################
## 2. CARGAR VARIANTES VALIDADAS (BLASTn core_nt o WGS)
############################################################
validated_raw <- read.table(val_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, check.names = FALSE)

# Eliminar columnas con nombre vacío o NA
validated_raw <- validated_raw[, nzchar(names(validated_raw)) & !is.na(names(validated_raw))]

# Ver nombres de columnas para confirmar
cat("Columnas:", paste(names(validated_raw), collapse = ", "), "\n")

validated_ids <- validated_raw %>%
  filter(VALIDACION_CORE == "SI" | `2_VALIDACION_WGS` == "SI") %>%
  mutate(VariantID = paste(Species, CHROM, POS, REF, ALT, sep = "|")) %>%
  pull(VariantID) %>%
  unique()

############################################################
## 3. LEER FICHEROS GT -> MEGA MATRIZ
############################################################
read_gt <- function(f, ids) {
  df <- read_tsv(f,
                 col_types = cols(CHROM = col_character(), POS = col_integer(),
                                  REF   = col_character(), ALT = col_character(),
                                  .default = col_character()),
                 show_col_types = FALSE)
  species     <- sub("\\.QUAL_.*$", "", basename(f))
  sample_cols <- names(df)[5:ncol(df)]
  X <- df %>%
    select(all_of(sample_cols)) %>%
    mutate(across(everything(), ~ as.numeric(na_if(.x, "."))))
  for (id in setdiff(ids, sample_cols)) X[[id]] <- NA_real_
  X <- X[, ids, drop = FALSE]
  bind_cols(
    df %>% transmute(
      Species   = species,
      CHROM     = CHROM, POS = POS, REF = REF, ALT = ALT,
      VariantID = paste(species, CHROM, POS, REF, ALT, sep = "|")
    ), X
  )
}

gt_files <- list.files(dir_gt, pattern = "\\.GT\\.tsv$", full.names = TRUE)
mega_raw <- bind_rows(lapply(gt_files, read_gt, ids = meta$sample))

cat("Variantes totales (antes de filtro validación):", nrow(mega_raw), "\n")

############################################################
## 4. FILTRAR A VARIANTES VALIDADAS
############################################################
mega <- mega_raw %>% filter(VariantID %in% validated_ids)
cat("Variantes tras filtro validación:", nrow(mega), "\n")

Xmega    <- as.matrix(mega[, -(1:6)])
rownames(Xmega) <- mega$VariantID

cat("Variantes:", nrow(Xmega), "| Individuos:", ncol(Xmega),
    "| %NA:", round(mean(is.na(Xmega)) * 100, 1), "%\n")

############################################################
## 5. COBERTURA POR ESPECIE Y PAÍS
############################################################
cobertura <- mega %>%
  select(Species, VariantID, all_of(meta$sample)) %>%
  pivot_longer(cols = all_of(meta$sample), names_to = "sample", values_to = "GT") %>%
  left_join(meta %>% select(sample, Country), by = "sample") %>%
  group_by(Species, Country) %>%
  summarise(n_individuos = n_distinct(sample[!is.na(GT)]), .groups = "drop") %>%
  pivot_wider(names_from = Country, values_from = n_individuos)
print(cobertura)
write_csv(cobertura, file.path(tab_dir, "cobertura_bacteria_pais.csv"))

############################################################
## 6. ORIGEN DE LOS NAs
############################################################
reads_mat <- read_tsv(reads_file, show_col_types = FALSE) %>%
  column_to_rownames("Species")

is_ral_tmp   <- grepl("^Ralstonia", mega$Species)
species_main <- intersect(unique(mega$Species[!is_ral_tmp]), rownames(reads_mat))
reads_main   <- reads_mat[species_main, , drop = FALSE]

n_tot <- prod(dim(reads_main))
n_0   <- sum(reads_main == 0)
n_low <- sum(reads_main > 0 & reads_main < 5)
cat("Origen de NAs — ausencia real (reads=0):", sprintf("%.1f%%", 100*n_0/n_tot),
    "| baja cobertura (1-4 reads):", sprintf("%.1f%%\n", 100*n_low/n_tot))

prev_especie <- data.frame(
  Species    = species_main,
  n_covered  = rowSums(reads_main >= 5),
  n_absent   = rowSums(reads_main == 0),
  n_low      = rowSums(reads_main > 0 & reads_main < 5),
  prevalence = round(100 * rowSums(reads_main >= 5) / ncol(reads_main), 1)
) %>% arrange(desc(prevalence))
print(prev_especie)
write_csv(prev_especie, file.path(tab_dir, "prevalencia_real_por_especie.csv"))

############################################################
## 7. FILTROS QC
##   - Variantes con >60% NA eliminadas
##   - Individuos con <10 variantes observadas eliminados
############################################################
prop_na_var <- rowMeans(is.na(Xmega))
Xa          <- Xmega[prop_na_var < 0.60, ]
mega_a      <- mega[prop_na_var < 0.60, ]
n_obs_ind   <- colSums(!is.na(Xa))
Xa          <- Xa[, n_obs_ind >= 10, drop = FALSE]
cat("Tras QC — Variantes:", nrow(Xa), "| Individuos:", ncol(Xa),
    "| %NA:", round(mean(is.na(Xa)) * 100, 1), "%\n")

############################################################
## 8. SEPARAR CON Y SIN RALSTONIA
############################################################
is_ral    <- grepl("^Ralstonia", rownames(Xa))
Xa_main   <- Xa[!is_ral, ]
mega_main <- mega_a[!is_ral, ]
cat("Sin Ralstonia — Variantes:", nrow(Xa_main), "| Individuos:", ncol(Xa_main), "\n")

############################################################
## 9. METADATA QC Y MODELO
############################################################
meta_qc    <- meta %>% filter(sample %in% colnames(Xa_main))
meta_model <- meta_qc %>% filter(!is.na(Age), !is.na(Sex))
cat("meta_qc:", nrow(meta_qc), "| meta_model (sin NA Age/Sex):", nrow(meta_model), "\n")

############################################################
## 10. GUARDAR RDS
############################################################
# Xmega_raw  → matriz completa (117 variantes × 337 individuos, con Ralstonia)
# Xa_qc      → matriz tras QC (117 variantes × 337 individuos, con Ralstonia)
# Xa_main    → matriz principal (105 variantes × 337 individuos, sin Ralstonia)
saveRDS(meta,       file.path(rds_dir, "meta_raw.rds"))
saveRDS(meta_qc,    file.path(rds_dir, "meta_qc.rds"))
saveRDS(meta_model, file.path(rds_dir, "meta_model.rds"))
saveRDS(mega,       file.path(rds_dir, "mega_raw.rds"))
saveRDS(Xmega,      file.path(rds_dir, "Xmega_raw.rds"))
saveRDS(Xa,         file.path(rds_dir, "Xa_qc.rds"))
saveRDS(Xa_main,    file.path(rds_dir, "Xa_main.rds"))
saveRDS(mega_main,  file.path(rds_dir, "mega_main.rds"))

cat("Done. RDS guardados en:", rds_dir, "\n")
