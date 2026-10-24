library(sf)
library(dplyr)
library(tidyr)
library(leaflet)
library(arrow)
library(data.table)

 
# 3. Définition du dossier de travail (sans slash final)
DOSSIER <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"

# 4. Lecture du premier fichier Parquet (Chemin absolu)
# mobsco23 <- read_parquet("C:/Users/camil/Downloads/RP2023_mobsco.parquet")
# 
# 5. Affichage des colonnes pour vérification
# colnames(mobsco23)

# 6. Lecture du second fichier Parquet (Chemin construit proprement)
annuairelycees <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))
donnees_carreaux <- fread("G:/Mon Drive/Education/Master APE 25-27/Mémoire/france_metro/carreaux_1km_met.csv")



# =====================================================================
# CODE 1bis - Distances aux lycées (remplace le Code 1 et les blocs
# de distance du Code 3)
#
# Variables créées (en km, depuis le centroïde de population de la commune) :
#   dist_lycee_general_km        lycée proposant la voie générale (public ou privé)
#   dist_lycee_public_km         lycée public (toute voie)
#   dist_lycee_prive_km          lycée privé (toute voie), y compris dans la commune
#   dist_lycee_prive_autre_km    lycée privé le plus proche HORS de la commune
#   dist_lycee_public_general_km lycée public proposant la voie générale
#   dist_lycee_prive_general_km  lycée privé proposant la voie générale
#   dist_lycee_public_gt_km      lycée public général ou techno (contrôle du Code 4)
#
# Corrections par rapport au script d'origine :
#  1. Seuls les LYCÉES OUVERTS sont retenus. Le filtre précédent
#     (statut == "Public") incluait écoles, collèges, services administratifs...
#  2. Centroïdes calculés à partir des carreaux Filosofi (pondérés par la
#     population) : les arrondissements de Paris/Lyon/Marseille ne sont plus perdus
#     (le geojson des centroïdes n'a que des communes entières).
#  3. Distance au privé le plus proche calculée aussi HORS commune : voir note.
#
# NOTE IMPORTANTE (endogénéité) : si la commune contient un lycée privé, la
# distance au privé le plus proche est ~0 : c'est quasi mécaniquement une
# fonction de ta variable dépendante (nb_lycees_prives > 0). Ne l'utilise PAS
# comme régresseur. Pour les modèles, prends dist_lycee_prive_autre_km (ou le
# nombre de privés dans un rayon, déjà dans le Code 3).
# =====================================================================

library(dplyr)
library(stringr)
library(arrow)
library(sf)

# ---------------------------------------------------------------------
# 0. PARAMÈTRES (AGREGER_PLM : identique aux Codes 2 et 3)
# ---------------------------------------------------------------------
DOSSIER      <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
DOSSIER_ECO  <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/data_eco"
DOSSIER_OUT  <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
AGREGER_PLM  <- FALSE
CENTROIDES_DEPUIS_CARREAUX <- TRUE    # FALSE = utiliser le geojson (perd les arrondissements)
CHEMIN_CENTROIDES <- "C:/Users/camil/Downloads/centroides-communes-population.geojson"
DEMI_CARREAU <- 100   # m ; ajouté à x_l93/y_l93 si ces coordonnées sont le coin sud-ouest du carreau
# (mettre 0 si ce sont déjà les centres)

map_plm <- function(code) {
  dplyr::case_when(
    grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
    grepl("^6938[1-9]$", code)             ~ "69123",
    grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
    TRUE ~ code)
}
harmoniser_code <- function(x) {
  code <- str_pad(as.character(x), width = 5, pad = "0")
  if (AGREGER_PLM) map_plm(code) else code
}
flag <- function(x) as.character(x) %in% c("1", "TRUE", "true", "Oui", "OUI", "oui")

# ---------------------------------------------------------------------
# 1. LYCÉES (ouverts, avec coordonnées)
# ---------------------------------------------------------------------
annuairelycees <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))

lycees <- annuairelycees %>%
  filter(grepl("^Lyc", type_etablissement),          # uniquement les lycées
         toupper(etat) == "OUVERT",
         !is.na(latitude), !is.na(longitude)) %>%
  mutate(
    code_insee = harmoniser_code(code_commune),
    secteur = case_when(statut_public_prive == "Public"     ~ "public",
                        grepl("^Priv", statut_public_prive) ~ "prive",
                        TRUE ~ NA_character_),
    general = flag(voie_generale),
    gt      = flag(voie_generale) | flag(voie_technologique)
  ) %>%
  filter(!is.na(secteur)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(2154)

cat("Lycées ouverts géolocalisés :", nrow(lycees),
    "| publics :", sum(lycees$secteur == "public"),
    "| privés :", sum(lycees$secteur == "prive"),
    "| avec voie générale :", sum(lycees$general), "\n")

# ---------------------------------------------------------------------
# 2. CENTROÏDES DE POPULATION DES COMMUNES DE L'ÉCHANTILLON
# ---------------------------------------------------------------------
df_analyse <- read_parquet(file.path(DOSSIER_OUT, "df_analyse.parquet"))
codes_cibles <- unique(df_analyse$code_insee)

if (CENTROIDES_DEPUIS_CARREAUX) {
  files <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*\\.parquet", full.names = TRUE)
  cent_df <- bind_rows(lapply(files, function(f) {
    read_parquet(f, col_select = c("lcog_geo", "x_l93", "y_l93", "ind")) %>%
      mutate(code_insee = harmoniser_code(substr(as.character(lcog_geo), 1, 5))) %>%
      filter(code_insee %in% codes_cibles, is.finite(ind), ind > 0) %>%
      group_by(code_insee) %>%
      summarise(x = weighted.mean(x_l93 + DEMI_CARREAU, ind),
                y = weighted.mean(y_l93 + DEMI_CARREAU, ind), .groups = "drop")
  }))
  centroides <- st_as_sf(cent_df, coords = c("x", "y"), crs = 2154)
} else {
  centroides <- st_read(CHEMIN_CENTROIDES, quiet = TRUE) %>%
    rename(code_insee = INSEE_COM) %>%
    mutate(code_insee = harmoniser_code(code_insee)) %>%
    filter(code_insee %in% codes_cibles) %>%
    st_transform(2154)
}
cat("Centroïdes calculés :", nrow(centroides), "sur", length(codes_cibles), "communes\n")

# ---------------------------------------------------------------------
# 3. DISTANCES
#    Matrice communes x lycées (quelques centaines x quelques milliers : OK).
#    Pour toutes les communes de France (35 000), utiliser st_nearest_feature.
# ---------------------------------------------------------------------
dist_min_km <- function(cent, cible, exclure_meme_commune = FALSE) {
  D <- st_distance(cent, cible)
  D <- matrix(as.numeric(D), nrow = nrow(cent))
  if (exclure_meme_commune)
    D[outer(cent$code_insee, cible$code_insee, "==")] <- Inf
  apply(D, 1, min) / 1000
}

sous <- function(...) lycees[with(lycees, ...), ]

cible_general        <- sous(general)
cible_public         <- sous(secteur == "public")
cible_prive          <- sous(secteur == "prive")
cible_public_general <- sous(secteur == "public" & general)
cible_prive_general  <- sous(secteur == "prive" & general)
cible_public_gt      <- sous(secteur == "public" & gt)

distances <- centroides %>%
  mutate(
    dist_lycee_general_km        = dist_min_km(., cible_general),
    dist_lycee_public_km         = dist_min_km(., cible_public),
    dist_lycee_prive_km          = dist_min_km(., cible_prive),
    dist_lycee_prive_autre_km    = dist_min_km(., cible_prive, exclure_meme_commune = TRUE),
    dist_lycee_public_general_km = dist_min_km(., cible_public_general),
    dist_lycee_prive_general_km  = dist_min_km(., cible_prive_general),
    dist_lycee_public_gt_km      = dist_min_km(., cible_public_gt)
  ) %>%
  st_drop_geometry()

# ---------------------------------------------------------------------
# 4. FUSION (remplace les anciennes colonnes de distance de df_analyse)
# ---------------------------------------------------------------------
anciennes <- grep("^dist_lycee_", names(df_analyse), value = TRUE)
df_analyse <- df_analyse %>%
  select(-all_of(anciennes)) %>%
  left_join(distances, by = "code_insee")

write_parquet(df_analyse, file.path(DOSSIER_OUT, "df_analyse.parquet"))

# ---------------------------------------------------------------------
# 5. CONTRÔLES
# ---------------------------------------------------------------------
cat("\n=== CONTRÔLES ===\n")
cat("Communes sans distance (NA) :", sum(is.na(df_analyse$dist_lycee_public_km)), "\n")
cat("Résumé des distances (km) :\n")
print(summary(df_analyse[, grep("^dist_lycee_", names(df_analyse))]))

# Cohérence : une commune avec un lycée privé doit avoir une distance au privé proche de 0
cat("\nDistance médiane au privé le plus proche (km) :\n")
cat("  communes avec >= 1 lycée privé :",
    round(median(df_analyse$dist_lycee_prive_km[df_analyse$nb_lycees_prives > 0], na.rm = TRUE), 2), "\n")
cat("  communes sans lycée privé      :",
    round(median(df_analyse$dist_lycee_prive_km[df_analyse$nb_lycees_prives == 0], na.rm = TRUE), 2), "\n")
cat("-> la première valeur doit être très faible : dist_lycee_prive_km est ENDOGÈNE à y.\n")

cat("\nCorrélation de la distance au privé (hors commune) avec nb_lycees_prives :",
    round(cor(df_analyse$dist_lycee_prive_autre_km, df_analyse$nb_lycees_prives,
              use = "complete.obs"), 3), "\n")
cat("Communes ayant un lycée public mais distance publique > 5 km (erreur de géocodage ?) :",
    sum(df_analyse$nb_lycees_publics > 0 & df_analyse$dist_lycee_public_km > 5, na.rm = TRUE), "\n")