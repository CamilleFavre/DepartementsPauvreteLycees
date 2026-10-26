# =====================================================================
# CODE 0b - CARTES DES COLLÈGES (à lancer UNE FOI en amont)
# Produit : carte_colleges_prives.png et carte_colleges_publics.png
# Le .qmd se contente ensuite d'inclure les images -> rendu instantané.
# =====================================================================
library(arrow); library(dplyr); library(tidyr); library(sf)
library(ggplot2); library(scales); library(stringr)

os_name <- Sys.info()[["sysname"]]
if (os_name == "Windows") {
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
  DOSSIER_QMD <- DOSSIER_OUT   # <- adapter : dossier du .qmd
} else if (os_name == "Darwin") {
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
  DOSSIER_QMD <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"  # <- adapter : dossier du .qmd
} else stop("Système d'exploitation non pris en charge.")

# ---------------------------------------------------------------------
# 1. CONTOURS (caches light du Code 0) et DONNÉES
# ---------------------------------------------------------------------
communes <- readRDS(file.path(DOSSIER_OUT, "cache_communes.rds")) %>%
  mutate(code_insee = str_pad(as.character(code_insee), 5, "left"))
regions  <- readRDS(file.path(DOSSIER_OUT, "cache_regions.rds"))

map_plm <- function(code) {
  dplyr::case_when(
    grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
    grepl("^6938[1-9]$", code)             ~ "69123",
    grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
    TRUE ~ code
  )
}

d_cartes <- read_parquet(file.path(DOSSIER_OUT, "df_analyse_colleges.parquet")) %>%
  mutate(code_insee = map_plm(str_pad(as.character(code_insee), 5, "left"))) %>%
  group_by(code_insee) %>%
  summarise(
    nb_colleges_prives  = sum(coalesce(nb_colleges_prives, 0L)),
    nb_colleges_publics = sum(coalesce(nb_colleges_publics, 0L)),
    .groups = "drop"
  )

# ---------------------------------------------------------------------
# 2. FONCTION DE CARTE (identique à la version Quarto)
# ---------------------------------------------------------------------
carte_colleges <- function(var, titre_legende) {
  carte <- communes %>%
    left_join(d_cartes, by = "code_insee") %>%
    mutate(
      nb = replace_na(.data[[var]], 0),
      tranche = case_when(
        nb == 0  ~ "0", nb == 1  ~ "1", nb == 2  ~ "2",
        nb == 3  ~ "3", nb == 4  ~ "4",
        nb >= 5 & nb <= 10 ~ "5-10", nb >= 11 ~ "11+"
      ),
      tranche = factor(tranche, levels = c("0", "1", "2", "3", "4", "5-10", "11+"))
    )
  
  ggplot() +
    geom_sf(data = carte, aes(fill = tranche), color = NA, linewidth = 0) +
    scale_fill_manual(values = c("0" = "white", "1" = "#fde725", "2" = "#aadc32",
                                 "3" = "#5ec962", "4" = "#21918c", "5-10" = "#3b528b",
                                 "11+" = "#440154"),
                      name = titre_legende, na.translate = FALSE) +
    geom_sf(data = regions, fill = NA, color = "black", linewidth = 0.3, alpha = 0.7) +
    coord_sf(xlim = c(80000, 1100000), ylim = c(6040000, 7120000),
             expand = FALSE, datum = NA) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1,
                               override.aes = list(color = "grey80", linewidth = 0.2))) +
    theme_void(base_family = "sans") +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 10, face = "bold", colour = "grey10"),
          legend.text = element_text(size = 8, colour = "grey20"),
          legend.margin = margin(t = 10),
          plot.margin = margin(10, 10, 10, 10),
          plot.background = element_rect(fill = "white", colour = NA))
}

# ---------------------------------------------------------------------
# 3. EXPORT PNG (300 dpi, ~25x22 cm : qualité publication)
# ---------------------------------------------------------------------
p1 <- carte_colleges("nb_colleges_prives",  "Number of private collèges")
p2 <- carte_colleges("nb_colleges_publics", "Number of public collèges")

ggsave(file.path(DOSSIER_OUT, "carte_colleges_prives.png"),
       p1, width = 9, height = 8, dpi = 300, bg = "white")
ggsave(file.path(DOSSIER_OUT, "carte_colleges_publics.png"),
       p2, width = 9, height = 8, dpi = 300, bg = "white")

cat("PNG écrits dans :", DOSSIER_QMD, "\n")