# 02_global_structure.R - Estructura global de la comunidad
# Input:  RDS de 01_prepare_data.R
# Output: PCoA_sin_Ralstonia.pdf, dispersion_por_pais.pdf,
#         PERMANOVA_principal.csv, PERMANOVA_chemistry.csv,
#         PERMANOVA_iberico.csv, tabla_sensibilidad_chemistry.csv
#         PERMANOVA_pairwise.csv, tabla_sensibilidad_Ralstonia.csv

set.seed(42)
library(tidyverse); library(vegan); library(ape); library(ggrepel)

############################################################
## RUTAS
############################################################
base_dir <- "/media/apacon/SSD_APACON/1_SCRIPTS_FINALES_TFM_run_final/GROND_nr_Rep"

out_dir <- file.path(base_dir, "metadata")
fig_dir <- file.path(out_dir, "figures", "phase2")
tab_dir <- file.path(out_dir, "tables")
rds_dir <- file.path(out_dir, "rds")
val_file <- file.path(tab_dir, "validated_variants_blast.tsv")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

meta_model <- readRDS(file.path(rds_dir, "meta_model.rds"))
Xa         <- readRDS(file.path(rds_dir, "Xa_qc.rds"))
Xa_main    <- readRDS(file.path(rds_dir, "Xa_main.rds"))

# Paleta 
COL <- c("Norway" = "#E69F00", "Portugal" = "#CC79A7", "Spain" = "#0072B2")
SHP <- c("Man" = 16, "Woman" = 17)

COUNTRY_ES <- c("Norway" = "Noruega", "Portugal" = "Portugal", "Spain" = "España")
SEX_ES     <- c("Man" = "Hombre", "Woman" = "Mujer")

############################################################
## JACCARD NA-aware
############################################################
calc_jaccard_na <- function(Xmat) {
  n    <- ncol(Xmat)
  nams <- colnames(Xmat)
  D     <- matrix(0, n, n, dimnames = list(nams, nams))
  Nused <- matrix(0, n, n, dimnames = list(nams, nams))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      a  <- Xmat[, i]; b <- Xmat[, j]
      ok <- !is.na(a) & !is.na(b)
      Nused[i, j] <- Nused[j, i] <- sum(ok)
      a <- a[ok]; b <- b[ok]
      inter   <- sum(a == 1 & b == 1)
      uni     <- sum(a == 1 | b == 1)
      D[i, j] <- D[j, i] <- ifelse(uni == 0, 0, 1 - inter / uni)
    }
  }
  list(dist = as.dist(D), Nused = Nused)
}

res_main <- calc_jaccard_na(Xa_main)
used     <- res_main$Nused[upper.tri(res_main$Nused)]
cat("Jaccard — % pares con <10 variantes comparables:",
    round(mean(used < 10) * 100, 1), "%  | mediana:", median(used), "\n")

############################################################
## PERMANOVA PRINCIPAL (sin Chemistry)
############################################################
D2    <- as.matrix(res_main$dist)[meta_model$sample, meta_model$sample]
dist2 <- as.dist(D2)

perm_main <- adonis2(
  dist2 ~ Country + Sex + Age,
  data = meta_model, permutations = 999, by = "margin"
)
print(perm_main)
write_csv(as.data.frame(perm_main) %>% rownames_to_column("Term"),
          file.path(tab_dir, "PERMANOVA_principal.csv"))

############################################################
## BETADISPER
############################################################
disp       <- betadisper(dist2, meta_model$Country)
disp_perm  <- permutest(disp)
disp_tukey <- TukeyHSD(disp)
print(disp_perm); print(disp_tukey)
capture.output(disp_perm,  file = file.path(tab_dir, "betadisper_permutest.txt"))
capture.output(disp_tukey, file = file.path(tab_dir, "betadisper_tukey.txt"))

dist_centroid <- tibble(
  sample   = meta_model$sample,
  Country  = meta_model$Country,
  distance = disp$distances
)
write_csv(dist_centroid, file.path(tab_dir, "distancia_centroide_por_muestra.csv"))
cat("Media distancia al centroide por país:\n")
print(round(tapply(disp$distances, meta_model$Country, mean), 3))

############################################################
## PCoA
############################################################
p_main <- pcoa(dist2, correction = "cailliez")
coords <- as.data.frame(p_main$vectors[, 1:2]) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("sample") %>%
  left_join(meta_model, by = "sample") %>%
  mutate(
    Country_es = recode(Country, !!!COUNTRY_ES),
    Sex_es     = recode(Sex,     !!!SEX_ES)
  )

pct1 <- round(p_main$values$Rel_corr_eig[1] * 100, 1)
pct2 <- round(p_main$values$Rel_corr_eig[2] * 100, 1)

n_country      <- coords %>% count(Country_es) %>% deframe()
country_labels <- setNames(
  paste0(names(n_country), " (n=", n_country, ")"),
  names(n_country)
)

# Paleta y formas con nombres 
COL_ES <- setNames(COL, COUNTRY_ES[names(COL)])
SHP_ES <- setNames(SHP, SEX_ES[names(SHP)])

ggplot(coords, aes(PCoA1, PCoA2, color = Country_es, shape = Sex_es)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Country_es), level = 0.95, linetype = 2) +
  scale_color_manual(values = COL_ES, labels = country_labels) +
  scale_shape_manual(values = SHP_ES) +
  theme_bw() +
  labs(
    title = expression(
      bold("PCoA - distancia Jaccard NA-aware") ~ italic("(sin Ralstonia insidiosa)")
    ),
    x = paste0("PCoA1 (", pct1, "%)"),
    y = paste0("PCoA2 (", pct2, "%)"),
    color = "País",
    shape = "Sexo"
  )
ggsave(file.path(fig_dir, "PCoA_sin_Ralstonia.png"), width = 8, height = 6, dpi = 300)

# Boxplot betadisper
dist_centroid_es <- dist_centroid %>%
  mutate(Country_es = recode(Country, !!!COUNTRY_ES))

n_country_box      <- dist_centroid_es %>% count(Country_es) %>% deframe()
country_labels_box <- setNames(
  paste0(names(n_country_box), " (n=", n_country_box, ")"),
  names(n_country_box)
)

ggplot(dist_centroid_es, aes(Country_es, distance, fill = Country_es)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Country_es), width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values  = COL_ES, labels = country_labels_box) +
  scale_color_manual(values = COL_ES, labels = country_labels_box) +
  scale_x_discrete(labels = country_labels_box) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Dispersión beta-diversidad por país (betadisper)",
    y     = "Distancia al centroide (Jaccard)",
    x     = ""
  )
ggsave(file.path(fig_dir, "dispersion_por_pais.png"), width = 6, height = 5, dpi = 300)

############################################################
## SENSIBILIDAD: con vs sin Ralstonia
############################################################
res_full  <- calc_jaccard_na(Xa)
Df        <- as.matrix(res_full$dist)[meta_model$sample, meta_model$sample]
perm_full <- adonis2(as.dist(Df) ~ Country + Sex + Age,
                     data = meta_model, permutations = 999, by = "margin")

tabla_sens <- tibble(
  Analisis    = c("Con Ralstonia (14 spp)", "Sin Ralstonia (13 spp)"),
  N_variantes = c(nrow(Xa), nrow(Xa_main)),
  Country_R2  = c(round(perm_full["Country", "R2"], 3),  round(perm_main["Country", "R2"], 3)),
  Country_p   = c(perm_full["Country", "Pr(>F)"],        perm_main["Country", "Pr(>F)"]),
  Age_R2      = c(round(perm_full["Age", "R2"], 3),      round(perm_main["Age", "R2"], 3)),
  Age_p       = c(perm_full["Age", "Pr(>F)"],            perm_main["Age", "Pr(>F)"]),
  Sex_p       = c(perm_full["Sex", "Pr(>F)"],            perm_main["Sex", "Pr(>F)"])
)
print(tabla_sens)
write_csv(tabla_sens, file.path(tab_dir, "tabla_sensibilidad_Ralstonia.csv"))

############################################################
## SENSIBILIDAD: efecto química de secuenciación
############################################################

# 1. PERMANOVA con Chemistry (confirmación colinealidad)
perm_chemistry <- adonis2(
  dist2 ~ Country + Sex + Age + Chemistry,
  data = meta_model, permutations = 999, by = "margin"
)
print(perm_chemistry)
write_csv(as.data.frame(perm_chemistry) %>% rownames_to_column("Term"),
          file.path(tab_dir, "PERMANOVA_chemistry.csv"))

# 2. PERMANOVA solo ibéricos (misma química R9.4.1)
ids_iberic  <- meta_model$sample[meta_model$Country %in% c("Spain", "Portugal")]
meta_iberic <- meta_model %>% filter(sample %in% ids_iberic) %>% droplevels()
D_iberic    <- as.matrix(res_main$dist)[ids_iberic, ids_iberic]

perm_iberic <- adonis2(
  as.dist(D_iberic) ~ Country + Sex + Age,
  data = meta_iberic, permutations = 999, by = "margin"
)
print(perm_iberic)
write_csv(as.data.frame(perm_iberic) %>% rownames_to_column("Term"),
          file.path(tab_dir, "PERMANOVA_iberico.csv"))

# 3. Tabla comparativa
tabla_chemistry <- tibble(
  Analisis   = c("3 países sin Chemistry",
                 "3 países con Chemistry",
                 "Solo ibéricos (R9.4.1)"),
  N          = c(332, 332, 239),
  Country_R2 = c(round(perm_main["Country",     "R2"], 3),   # 0.076
                 round(perm_chemistry["Country", "R2"], 3),   # 0.021
                 round(perm_iberic["Country",    "R2"], 3)),  # 0.029
  Country_F  = c(round(perm_main["Country",     "F"], 3),
                 round(perm_chemistry["Country", "F"], 3),
                 round(perm_iberic["Country",    "F"], 3)),
  Country_p  = c(perm_main["Country",     "Pr(>F)"],
                 perm_chemistry["Country","Pr(>F)"],
                 perm_iberic["Country",   "Pr(>F)"])
)
print(tabla_chemistry)
write_csv(tabla_chemistry, file.path(tab_dir, "tabla_sensibilidad_chemistry.csv"))

############################################################
## PERMANOVA POR PARES
############################################################
pares <- list(c("Norway", "Portugal"), c("Norway", "Spain"), c("Portugal", "Spain"))

tabla_pares <- bind_rows(lapply(pares, function(par) {
  ids      <- meta_model$sample[meta_model$Country %in% par]
  meta_par <- meta_model %>% filter(sample %in% ids) %>% droplevels()
  D_par    <- as.matrix(res_main$dist)[ids, ids]
  perm_par <- adonis2(as.dist(D_par) ~ Country + Sex + Age,
                      data = meta_par, permutations = 999, by = "margin")
  tibble(
    Comparacion = paste(COUNTRY_ES[par], collapse = " vs "),
    R2          = round(perm_par["Country", "R2"], 3),
    F           = round(perm_par["Country", "F"],  3),
    p_valor     = perm_par["Country", "Pr(>F)"]
  )
})) %>% mutate(p_bonferroni = p.adjust(p_valor, method = "bonferroni"))

print(tabla_pares)
write_csv(tabla_pares, file.path(tab_dir, "PERMANOVA_pairwise.csv"))

############################################################
## GUARDAR RDS
############################################################
saveRDS(res_main,        file.path(rds_dir, "jaccard_main.rds"))
saveRDS(coords,          file.path(rds_dir, "coords_pcoa.rds"))
saveRDS(perm_main,       file.path(rds_dir, "perm_main.rds"))
saveRDS(perm_chemistry,  file.path(rds_dir, "perm_chemistry.rds"))
saveRDS(perm_iberic,     file.path(rds_dir, "perm_iberic.rds"))
saveRDS(disp,            file.path(rds_dir, "betadisper_main.rds"))
saveRDS(country_labels,  file.path(rds_dir, "country_labels.rds"))
saveRDS(list(pct1 = pct1, pct2 = pct2), file.path(rds_dir, "pcoa_axis_pct.rds"))

