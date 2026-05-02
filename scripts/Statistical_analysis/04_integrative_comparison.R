# 04_integrative_comparison.R - Comparación integrativa
# Input:  RDS de 02 y 03
# Output: biplot_principal_sin_Ralstonia.png, scatter_envfit_vs_AF.png,
#         tabla_concordancia_variantes.csv

library(tidyverse); library(vegan); library(ggrepel)

############################################################
## RUTAS
############################################################
out_dir <- "/media/apacon/SSD_APACON/1_SCRIPTS_FINALES_TFM_run_final/GROND_nr_Rep/metadata"
fig_dir <- file.path(out_dir, "figures", "phase2")
tab_dir <- file.path(out_dir, "tables")
rds_dir <- file.path(out_dir, "rds")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

meta_model  <- readRDS(file.path(rds_dir, "meta_model.rds"))
Xa_main     <- readRDS(file.path(rds_dir, "Xa_main.rds"))
coords      <- readRDS(file.path(rds_dir, "coords_pcoa.rds"))
tabla_final <- readRDS(file.path(rds_dir, "tabla_final.rds"))
pct         <- readRDS(file.path(rds_dir, "pcoa_axis_pct.rds"))

# Paleta 
COL <- c("Noruega" = "#E69F00", "Portugal" = "#CC79A7", "España" = "#0072B2")
SHP <- c("Hombre" = 16, "Mujer" = 17)

COUNTRY_ES <- c("Norway" = "Noruega", "Portugal" = "Portugal", "Spain" = "España")
SEX_ES     <- c("Man" = "Hombre", "Woman" = "Mujer")

# Metadatos en coords
coords <- coords %>%
  mutate(
    Country_es = recode(Country, !!!COUNTRY_ES),
    Sex_es     = recode(Sex,     !!!SEX_ES)
  )

n_country      <- coords %>% count(Country_es) %>% deframe()
country_labels <- setNames(
  paste0(names(n_country), " (n=", n_country, ")"),
  names(n_country)
)

############################################################
## ENVFIT
############################################################
Xv <- Xa_main[, meta_model$sample]
for (k in seq_len(nrow(Xv))) {
  g <- Xv[k, ]; m <- mean(g, na.rm = TRUE)
  Xv[k, is.na(g)] <- ifelse(is.nan(m), 0, m)
}

fit <- envfit(coords[, c("PCoA1", "PCoA2")], t(Xv), permutations = 999)

envfit_df <- as.data.frame(scores(fit, display = "vectors")) %>%
  rownames_to_column("VariantID") %>%
  mutate(
    envfit_r2  = fit$vectors$r,
    envfit_p   = fit$vectors$pvals,
    envfit_fdr = p.adjust(envfit_p, "fdr")
  )

cat("Variantes significativas envfit FDR<0.05:",
    sum(envfit_df$envfit_fdr < 0.05, na.rm = TRUE), "\n")

############################################################
## BIPLOT - TOP 10
############################################################
fdf <- envfit_df %>%
  filter(envfit_fdr < 0.05) %>%
  arrange(desc(envfit_r2)) %>%
  slice_head(n = 10) %>%
  mutate(label = paste0(word(sapply(str_split(VariantID, "\\|"), `[`, 1), 1), "_",
                        sapply(str_split(VariantID, "\\|"), `[`, 3)))

ggplot() +
  geom_point(data = coords,
             aes(PCoA1, PCoA2, color = Country_es, shape = Sex_es),
             size = 2.5, alpha = 0.7) +
  stat_ellipse(data = coords,
               aes(PCoA1, PCoA2, color = Country_es, group = Country_es),
               level = 0.95, linetype = 2) +
  geom_segment(data = fdf,
               aes(x = 0, y = 0, xend = PCoA1 * 0.5, yend = PCoA2 * 0.5),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.8) +
  geom_text_repel(data = fdf,
                  aes(x = PCoA1 * 0.5, y = PCoA2 * 0.5, label = label),
                  size = 3.5, fontface = "bold", max.overlaps = 30) +
  scale_color_manual(values = COL, labels = country_labels) +
  scale_shape_manual(values = SHP) +
  theme_bw() +
  labs(
    title  = "Biplot PCoA — top 10 variantes asociadas a la separación entre países",
    x      = paste0("PCoA1 (", pct$pct1, "%)"),
    y      = paste0("PCoA2 (", pct$pct2, "%)"),
    color  = "País",
    shape  = "Sexo"
  )
ggsave(file.path(fig_dir, "biplot_principal_sin_Ralstonia.png"),width = 10, height = 7, dpi = 300)

############################################################
## TABLA DE CONCORDANCIA
############################################################
tabla_concordancia <- envfit_df %>%
  select(VariantID, envfit_r2, envfit_p, envfit_fdr) %>%
  left_join(tabla_final %>% select(VariantID, Species, max_diff, p_fisher, fdr,
                                   enriched_in, AF_Norway, AF_Portugal, AF_Spain),
            by = "VariantID") %>%
  mutate(
    sig_envfit  = envfit_fdr < 0.05,
    sig_fisher  = fdr < 0.05,
    concordante = sig_envfit & sig_fisher
  ) %>%
  arrange(desc(envfit_r2))

cat("Significativas envfit:", sum(tabla_concordancia$sig_envfit, na.rm = TRUE),
    "| Fisher:", sum(tabla_concordancia$sig_fisher, na.rm = TRUE),
    "| Concordantes:", sum(tabla_concordancia$concordante, na.rm = TRUE), "\n")

write_csv(tabla_concordancia, file.path(tab_dir, "tabla_concordancia_variantes.csv"))

############################################################
## SCATTER envfit R² vs max_diff AF 
############################################################

# Paleta dalton-safe para categorías
CAT_COL <- c(
  "Ambas (envfit + Fisher)" = "#D55E00",   # naranja-rojo
  "Solo envfit"             = "#E69F00",   # naranja
  "Solo Fisher"             = "#0072B2",   # azul
  "Ninguna"                 = "grey70"
)

scatter_df <- tabla_concordancia %>%
  filter(!is.na(max_diff)) %>%
  mutate(
    categoria = case_when(
      concordante ~ "Ambas (envfit + Fisher)",
      sig_envfit  ~ "Solo envfit",
      sig_fisher  ~ "Solo Fisher",
      TRUE        ~ "Ninguna"
    ),
    categoria = factor(categoria, levels = names(CAT_COL))
  )

ggplot(scatter_df, aes(x = max_diff, y = envfit_r2, color = categoria)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_text_repel(
    data = scatter_df %>% filter(concordante),
    aes(label = paste0(word(Species, 1), "_",
                       str_extract(VariantID, "(?<=\\|)\\d+(?=\\|)"))),
    size = 3, max.overlaps = 20, fontface = "italic"
  ) +
  scale_color_manual(values = CAT_COL) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    title = "Concordancia entre estructura global y diferencias por variante",
    x     = "Diferencia máxima de FA entre países (max_diff)",
    y     = expression("Correlación con el PCoA (envfit R"^2*")"),
    color = ""
  )
ggsave(file.path(fig_dir, "scatter_envfit_vs_AF.png"), width = 10, height = 7, dpi = 300)

