library(arrow)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(scales)
library(stringr)

DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
CHEMIN_COMMUNES <- "C:/Users/camil/Downloads/ADMIN-EXPRESS_4-0__GPKG_LAMB93_FXX_2026-10-15/ADMIN-EXPRESS_4-0__GPKG_LAMB93_FXX_2026-10-15/ADMIN-EXPRESS/1_DONNEES_LIVRAISON_2026-10-00164/ADE_4-0_GPKG_LAMB93_FXX-ED2026-10-15/ADE_4-0_GPKG_LAMB93_FXX-ED2026-10-15.gpkg"

# ---------------------------------------------------------------------
# 1. DONNÉES : agréger Paris/Lyon/Marseille vers la commune entière
# ---------------------------------------------------------------------
map_plm <- function(code) {
  dplyr::case_when(
    grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
    grepl("^6938[1-9]$", code)             ~ "69123",
    grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
    TRUE ~ code
  )
}

d <- read_parquet(file.path(DOSSIER_OUT, "df_analyse_all.parquet")) %>%
  mutate(code_insee = map_plm(str_pad(as.character(code_insee), 5, "left"))) %>%
  group_by(code_insee) %>%
  summarise(nb_lycees_prives = sum(coalesce(nb_lycees_prives, 0L)), .groups = "drop")

# ---------------------------------------------------------------------
# 2. CONTOURS : Charger TOUTES les communes et les régions
# ---------------------------------------------------------------------
communes <- st_read(CHEMIN_COMMUNES, layer = "commune", quiet = TRUE) %>% 
  st_transform(2154)

regions <- st_read(CHEMIN_COMMUNES, layer = "region", quiet = TRUE) %>% 
  st_transform(2154)

col_code <- intersect(names(communes), c("INSEE_COM", "insee_com", "code_insee"))[1]
communes <- communes %>%
  mutate(code_insee = str_pad(as.character(.data[[col_code]]), 5, "left"))

# Jointure et découpage en classes discrètes selon vos tranches
carte <- communes %>% 
  left_join(d, by = "code_insee") %>%
  mutate(
    nb_lycees_prives = tidyr::replace_na(nb_lycees_prives, 0),
    tranche_lycees = case_when(
      nb_lycees_prives == 0  ~ "0",
      nb_lycees_prives == 1  ~ "1",
      nb_lycees_prives == 2  ~ "2",
      nb_lycees_prives == 3  ~ "3",
      nb_lycees_prives == 4  ~ "4",
      nb_lycees_prives >= 5 & nb_lycees_prives <= 10 ~ "5-10",
      nb_lycees_prives >= 11 ~ "11 et +"
    ),
    # Ordonner proprement les facteurs pour la légende
    tranche_lycees = factor(tranche_lycees, levels = c("0", "1", "2", "3", "4", "5-10", "11 et +"))
  )

# ---------------------------------------------------------------------
# 3. CARTE AVEC CLASSES DISCRÈTES ET FOND BLANC POUR 0
# ---------------------------------------------------------------------
# Définition des couleurs pour chaque classe (0 = blanc, puis dégradé visible du rose au violet foncé)
couleurs_classes <- c(
  "0"      = "white",
  "1"      = "#fde725", # Jaune vif uniquement pour les communes à 1 (ou changez si vous ne voulez pas de jaune du tout)
  "2"      = "#aadc32",
  "3"      = "#5ec962",
  "4"      = "#21918c",
  "5-10"   = "#3b528b",
  "11 et +"= "#440154"  # Violet foncé pour les très fortes concentrations
)

p <- ggplot() +
  # A. Rendu de toutes les communes par classe
  geom_sf(
    data = carte, 
    aes(fill = tranche_lycees), 
    color = NA, 
    linewidth = 0
  ) +
  
  # B. Application de la palette discrète
  scale_fill_manual(
    values = couleurs_classes,
    name   = "Nombre de lycées privés"
  ) +
  
  # C. Frontières régionales fines en NOIR
  geom_sf(data = regions, fill = NA, color = "black", linewidth = 0.3, alpha = 0.7) +
  
  coord_sf(xlim = c(80000, 1100000), ylim = c(6040000, 7120000),
           expand = FALSE, datum = NA) +
  
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5,
                             nrow = 1, override.aes = list(color = "grey80", linewidth = 0.2))) +
  theme_void(base_family = "sans") +
  theme(
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 10, face = "bold", colour = "grey10"),
    legend.text          = element_text(size = 8, colour = "grey20"),
    legend.margin        = margin(t = 10),
    plot.margin          = margin(10, 10, 10, 10),
    plot.background      = element_rect(fill = "white", colour = NA)
  )

print(p)











library(arrow)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(scales)
library(stringr)

DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
CHEMIN_COMMUNES <- "C:/Users/camil/Downloads/ADMIN-EXPRESS_4-0__GPKG_LAMB93_FXX_2026-10-15/ADMIN-EXPRESS_4-0__GPKG_LAMB93_FXX_2026-10-15/ADMIN-EXPRESS/1_DONNEES_LIVRAISON_2026-10-00164/ADE_4-0_GPKG_LAMB93_FXX-ED2026-10-15/ADE_4-0_GPKG_LAMB93_FXX-ED2026-10-15.gpkg"

# ---------------------------------------------------------------------
# 1. DONNÉES : agréger Paris/Lyon/Marseille vers la commune entière
# ---------------------------------------------------------------------
map_plm <- function(code) {
  dplyr::case_when(
    grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
    grepl("^6938[1-9]$", code)             ~ "69123",
    grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
    TRUE ~ code
  )
}

d <- read_parquet(file.path(DOSSIER_OUT, "df_analyse_all.parquet")) %>%
  mutate(code_insee = map_plm(str_pad(as.character(code_insee), 5, "left"))) %>%
  group_by(code_insee) %>%
  summarise(nb_lycees_publics = sum(coalesce(nb_lycees_publics, 0L)), .groups = "drop")

# ---------------------------------------------------------------------
# 2. CONTOURS : Charger TOUTES les communes et les régions
# ---------------------------------------------------------------------
communes <- st_read(CHEMIN_COMMUNES, layer = "commune", quiet = TRUE) %>% 
  st_transform(2154)

regions <- st_read(CHEMIN_COMMUNES, layer = "region", quiet = TRUE) %>% 
  st_transform(2154)

col_code <- intersect(names(communes), c("INSEE_COM", "insee_com", "code_insee"))[1]
communes <- communes %>%
  mutate(code_insee = str_pad(as.character(.data[[col_code]]), 5, "left"))

# Jointure et découpage en classes discrètes selon vos tranches
carte <- communes %>% 
  left_join(d, by = "code_insee") %>%
  mutate(
    nb_lycees_publics = tidyr::replace_na(nb_lycees_publics, 0),
    tranche_lycees = case_when(
      nb_lycees_publics == 0  ~ "0",
      nb_lycees_publics == 1  ~ "1",
      nb_lycees_publics == 2  ~ "2",
      nb_lycees_publics == 3  ~ "3",
      nb_lycees_publics == 4  ~ "4",
      nb_lycees_publics >= 5 & nb_lycees_publics <= 10 ~ "5-10",
      nb_lycees_publics >= 11 ~ "11 et +"
    ),
    # Ordonner proprement les facteurs pour la légende
    tranche_lycees = factor(tranche_lycees, levels = c("0", "1", "2", "3", "4", "5-10", "11 et +"))
  )

# ---------------------------------------------------------------------
# 3. CARTE AVEC CLASSES DISCRÈTES ET FOND BLANC POUR 0
# ---------------------------------------------------------------------
# Définition exacte des mêmes couleurs que pour les lycées privés
couleurs_classes <- c(
  "0"      = "white",
  "1"      = "#fde725",
  "2"      = "#aadc32",
  "3"      = "#5ec962",
  "4"      = "#21918c",
  "5-10"   = "#3b528b",
  "11 et +"= "#440154"
)

p <- ggplot() +
  # A. Rendu de toutes les communes par classe
  geom_sf(
    data = carte, 
    aes(fill = tranche_lycees), 
    color = NA, 
    linewidth = 0
  ) +
  
  # B. Application de la palette discrète
  scale_fill_manual(
    values = couleurs_classes,
    name   = "Nombre de lycées publics"
  ) +
  
  # C. Frontières régionales fines en NOIR
  geom_sf(data = regions, fill = NA, color = "black", linewidth = 0.3, alpha = 0.7) +
  
  coord_sf(xlim = c(80000, 1100000), ylim = c(6040000, 7120000),
           expand = FALSE, datum = NA) +
  
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5,
                             nrow = 1, override.aes = list(color = "grey80", linewidth = 0.2))) +
  theme_void(base_family = "sans") +
  theme(
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.title         = element_text(size = 10, face = "bold", colour = "grey10"),
    legend.text          = element_text(size = 8, colour = "grey20"),
    legend.margin        = margin(t = 10),
    plot.margin          = margin(10, 10, 10, 10),
    plot.background      = element_rect(fill = "white", colour = NA)
  )

print(p)