# 02_global_structure.R — Estructura global de la comunidad
# Input:  RDS de 01_prepare_data.R
# Output: PCoA_sin_Ralstonia.pdf, dispersion_por_pais.pdf,
#         PERMANOVA_principal.csv, PERMANOVA_pairwise.csv,
#         tabla_sensibilidad_Ralstonia.csv

set.seed(42)
library(tidyverse); library(vegan); library(ape); library(ggrepel)

out_dir <- "/media/apacon/SSD_APACON/1_SCRIPTS_FINALES_TFM/figures_phase2"
rds_dir <- file.path(out_dir, "rds")

meta_model <- readRDS(file.path(rds_dir, "meta_model.rds"))
Xa         <- readRDS(file.path(rds_dir, "Xa_qc.rds"))
Xa_main    <- readRDS(file.path(rds_dir, "Xa_main.rds"))

COL <- c("Norway" = "#E69F00", "Portugal" = "#CC79A7", "Spain" = "#0072B2")
SHP <- c("Man" = 16, "Woman" = 17)

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
## PERMANOVA
############################################################
D2    <- as.matrix(res_main$dist)[meta_model$sample, meta_model$sample]
dist2 <- as.dist(D2)

perm_main <- adonis2(
  dist2 ~ Country + Sex + Age,
  data = meta_model, permutations = 999, by = "margin"
)
print(perm_main)
write_csv(as.data.frame(perm_main) %>% rownames_to_column("Term"),
          file.path(out_dir, "PERMANOVA_principal.csv"))

############################################################
## BETADISPER
############################################################
disp       <- betadisper(dist2, meta_model$Country)
disp_perm  <- permutest(disp)
disp_tukey <- TukeyHSD(disp)
print(disp_perm); print(disp_tukey)
capture.output(disp_perm,  file = file.path(out_dir, "betadisper_permutest.txt"))
capture.output(disp_tukey, file = file.path(out_dir, "betadisper_tukey.txt"))

dist_centroid <- tibble(
  sample   = meta_model$sample,
  Country  = meta_model$Country,
  distance = disp$distances
)
write_csv(dist_centroid, file.path(out_dir, "distancia_centroide_por_muestra.csv"))
cat("Media distancia al centroide por país:\n")
print(round(tapply(disp$distances, meta_model$Country, mean), 3))

############################################################
## PCoA
############################################################
p_main <- pcoa(dist2, correction = "cailliez")
coords <- as.data.frame(p_main$vectors[, 1:2]) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("sample") %>%
  left_join(meta_model, by = "sample")

pct1 <- round(p_main$values$Rel_corr_eig[1] * 100, 1)
pct2 <- round(p_main$values$Rel_corr_eig[2] * 100, 1)

n_country      <- coords %>% count(Country) %>% deframe()
country_labels <- setNames(paste0(names(n_country), " (n=", n_country, ")"), names(n_country))

ggplot(coords, aes(PCoA1, PCoA2, color = Country, shape = Sex)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Country), level = 0.95, linetype = 2) +
  scale_color_manual(values = COL, labels = country_labels) +
  scale_shape_manual(values = SHP, labels = c("Man" = "Hombre", "Woman" = "Mujer")) +
  theme_bw() +
  labs(title = "PCoA — distancia Jaccard NA-aware (sin Ralstonia)",
       x = paste0("PCoA1 (", pct1, "%)"), y = paste0("PCoA2 (", pct2, "%)"),
       color = "País", shape = "Sexo")
ggsave(file.path(out_dir, "PCoA_sin_Ralstonia.pdf"), width = 8, height = 6)

ggplot(dist_centroid, aes(Country, distance, fill = Country)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Country), width = 0.15, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = COL) +
  scale_color_manual(values = COL) +
  scale_x_discrete(labels = country_labels) +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "Dispersión beta-diversidad por país (betadisper)",
       y = "Distancia al centroide (Jaccard)", x = "")
ggsave(file.path(out_dir, "dispersion_por_pais.pdf"), width = 6, height = 5)

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
  Country_R2  = c(round(perm_full["Country", "R2"], 3), round(perm_main["Country", "R2"], 3)),
  Country_p   = c(perm_full["Country", "Pr(>F)"],       perm_main["Country", "Pr(>F)"]),
  Age_R2      = c(round(perm_full["Age", "R2"], 3),     round(perm_main["Age", "R2"], 3)),
  Age_p       = c(perm_full["Age", "Pr(>F)"],           perm_main["Age", "Pr(>F)"]),
  Sex_p       = c(perm_full["Sex", "Pr(>F)"],           perm_main["Sex", "Pr(>F)"])
)
print(tabla_sens)
write_csv(tabla_sens, file.path(out_dir, "tabla_sensibilidad_Ralstonia.csv"))

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
    Comparacion = paste(par, collapse = " vs "),
    R2          = round(perm_par["Country", "R2"], 3),
    F           = round(perm_par["Country", "F"],  3),
    p_valor     = perm_par["Country", "Pr(>F)"]
  )
})) %>% mutate(p_bonferroni = p.adjust(p_valor, method = "bonferroni"))

print(tabla_pares)
write_csv(tabla_pares, file.path(out_dir, "PERMANOVA_pairwise.csv"))

############################################################
## GUARDAR RDS
############################################################
saveRDS(res_main,       file.path(rds_dir, "jaccard_main.rds"))
saveRDS(coords,         file.path(rds_dir, "coords_pcoa.rds"))
saveRDS(perm_main,      file.path(rds_dir, "perm_main.rds"))
saveRDS(disp,           file.path(rds_dir, "betadisper_main.rds"))
saveRDS(country_labels, file.path(rds_dir, "country_labels.rds"))
saveRDS(list(pct1 = pct1, pct2 = pct2), file.path(rds_dir, "pcoa_axis_pct.rds"))

