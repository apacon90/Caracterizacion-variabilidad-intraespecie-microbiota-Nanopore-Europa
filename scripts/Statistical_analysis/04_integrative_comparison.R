# 04_integrative_comparison.R — Comparación integrativa
# Input:  RDS de 02 y 03
# Output: biplot_principal_sin_Ralstonia.pdf, scatter_envfit_vs_AF.pdf,
#         tabla_concordancia_variantes.csv

library(tidyverse); library(vegan); library(ggrepel)

out_dir <- "/media/apacon/SSD_APACON/1_SCRIPTS_FINALES_TFM/figures_phase2"
rds_dir <- file.path(out_dir, "rds")

meta_model  <- readRDS(file.path(rds_dir, "meta_model.rds"))
Xa_main     <- readRDS(file.path(rds_dir, "Xa_main.rds"))
coords      <- readRDS(file.path(rds_dir, "coords_pcoa.rds"))
tabla_final <- readRDS(file.path(rds_dir, "tabla_final.rds"))
pct         <- readRDS(file.path(rds_dir, "pcoa_axis_pct.rds"))

COL <- c("Norway" = "#E69F00", "Portugal" = "#CC79A7", "Spain" = "#0072B2")
SHP <- c("Man" = 16, "Woman" = 17)

n_country      <- meta_model %>% count(Country) %>% deframe()
country_labels <- setNames(paste0(names(n_country), " (n=", n_country, ")"), names(n_country))

############################################################
## ENVFIT
## NAs imputados con la media por variante (requerimiento técnico de envfit)
############################################################
Xv <- Xa_main[, meta_model$sample]
for (k in seq_len(nrow(Xv))) {
  g <- Xv[k, ]; m <- mean(g, na.rm = TRUE)
  Xv[k, is.na(g)] <- ifelse(is.nan(m), 0, m)
}

fit <- envfit(coords[, c("PCoA1", "PCoA2")], t(Xv), permutations = 999)

envfit_df <- as.data.frame(scores(fit, display = "vectors")) %>%
  rownames_to_column("VariantID") %>%
  mutate(envfit_r2  = fit$vectors$r,
         envfit_p   = fit$vectors$pvals,
         envfit_fdr = p.adjust(envfit_p, "fdr"))

cat("Variantes significativas envfit FDR<0.05:",
    sum(envfit_df$envfit_fdr < 0.05, na.rm = TRUE), "\n")

############################################################
## BIPLOT — TOP 10
############################################################
fdf <- envfit_df %>%
  filter(envfit_fdr < 0.05) %>%
  arrange(desc(envfit_r2)) %>%
  slice_head(n = 10) %>%
  mutate(label = paste0(word(sapply(str_split(VariantID, "\\|"), `[`, 1), 1), "_",
                        sapply(str_split(VariantID, "\\|"), `[`, 3)))

ggplot() +
  geom_point(data = coords, aes(PCoA1, PCoA2, color = Country, shape = Sex),
             size = 2.5, alpha = 0.7) +
  stat_ellipse(data = coords, aes(PCoA1, PCoA2, color = Country, group = Country),
               level = 0.95, linetype = 2) +
  geom_segment(data = fdf, aes(x = 0, y = 0, xend = PCoA1 * 0.5, yend = PCoA2 * 0.5),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.8) +
  geom_text_repel(data = fdf, aes(x = PCoA1 * 0.5, y = PCoA2 * 0.5, label = label),
                  size = 3.5, fontface = "bold", max.overlaps = 30) +
  scale_color_manual(values = COL, labels = country_labels) +
  scale_shape_manual(values = SHP, labels = c("Man" = "Hombre", "Woman" = "Mujer")) +
  theme_bw() +
  labs(title = "Biplot PCoA — top 10 variantes asociadas a la separación entre países",
       x = paste0("PCoA1 (", pct$pct1, "%)"), y = paste0("PCoA2 (", pct$pct2, "%)"),
       color = "País", shape = "Sexo")
ggsave(file.path(out_dir, "biplot_principal_sin_Ralstonia.pdf"), width = 10, height = 7)

############################################################
## TABLA DE CONCORDANCIA
############################################################
tabla_concordancia <- envfit_df %>%
  select(VariantID, envfit_r2, envfit_p, envfit_fdr) %>%
  left_join(tabla_final %>% select(VariantID, Species, max_diff, p_fisher, fdr,
                                   enriched_in, AF_Norway, AF_Portugal, AF_Spain),
            by = "VariantID") %>%
  mutate(sig_envfit  = envfit_fdr < 0.05,
         sig_fisher  = fdr < 0.05,
         concordante = sig_envfit & sig_fisher) %>%
  arrange(desc(envfit_r2))

cat("Significativas envfit:", sum(tabla_concordancia$sig_envfit, na.rm = TRUE),
    "| Fisher:", sum(tabla_concordancia$sig_fisher, na.rm = TRUE),
    "| Concordantes:", sum(tabla_concordancia$concordante, na.rm = TRUE), "\n")

write_csv(tabla_concordancia, file.path(out_dir, "tabla_concordancia_variantes.csv"))

############################################################
## SCATTER envfit R² vs max_diff AF
############################################################
scatter_df <- tabla_concordancia %>%
  filter(!is.na(max_diff)) %>%
  mutate(
    categoria = case_when(
      concordante ~ "Ambas (envfit + Fisher)",
      sig_envfit  ~ "Solo envfit",
      sig_fisher  ~ "Solo Fisher",
      TRUE        ~ "Ninguna"
    ),
    categoria = factor(categoria, levels = c("Ambas (envfit + Fisher)",
                                             "Solo envfit", "Solo Fisher", "Ninguna"))
  )

ggplot(scatter_df, aes(x = max_diff, y = envfit_r2, color = categoria)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_text_repel(
    data = scatter_df %>% filter(concordante),
    aes(label = paste0(word(Species, 1), "_",
                       str_extract(VariantID, "(?<=\\|)\\d+(?=\\|)"))),
    size = 3, max.overlaps = 20, fontface = "italic"
  ) +
  scale_color_manual(values = c("Ambas (envfit + Fisher)" = "#D55E00",
                                "Solo envfit"             = "#E69F00",
                                "Solo Fisher"             = "#0072B2",
                                "Ninguna"                 = "grey70")) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(title = "Concordancia entre estructura global y diferencias por variante",
       x = "Diferencia máxima de AF entre países (max_diff)",
       y = "Correlación con el PCoA (envfit R²)", color = "")
ggsave(file.path(out_dir, "scatter_envfit_vs_AF.pdf"), width = 9, height = 7)

