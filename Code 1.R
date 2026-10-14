library(sf)
library(dplyr)
library(tidyr)
library(leaflet)
library(arrow)
library(data.table)

rm(list = ls())

# 3. Définition du dossier de travail (sans slash final)
DOSSIER <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"

# 4. Lecture du premier fichier Parquet (Chemin absolu)
mobsco23 <- read_parquet("C:/Users/camil/Downloads/RP2023_mobsco.parquet")

# 5. Affichage des colonnes pour vérification
colnames(mobsco23)

# 6. Lecture du second fichier Parquet (Chemin construit proprement)
annuairelycees <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))
donnees_carreaux <- fread("G:/Mon Drive/Education/Master APE 25-27/Mémoire/france_metro/carreaux_1km_met.csv")



# ---------------------------------------------------------------------
# 1. PRÉPARATION DES LYCÉES PUBLICS
# ---------------------------------------------------------------------
# (On suppose que ton objet 'annuairelycees' est déjà chargé dans ton environnement)

lycees_clean <- annuairelycees %>%
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

lycees_publics <- lycees_clean %>% 
  filter(statut_public_prive == "Public")


# ---------------------------------------------------------------------
# 2. CHARGEMENT ET HARMONISATION DES CENTROÏDES DE COMMUNE
# ---------------------------------------------------------------------
chemin_centroides <- "C:/Users/camil/Downloads/centroides-communes-population.geojson"
centroides <- st_read(chemin_centroides)

# Renommer la colonne pour harmoniser avec tes autres bases
centroides <- centroides %>%
  rename(code_commune = INSEE_COM)


# ---------------------------------------------------------------------
# 3. PROJECTION ET CALCUL DE LA DISTANCE AU LYCÉE PUBLIC LE PLUS PROCHE
# ---------------------------------------------------------------------
# On projette en Lambert 93 (EPSG:2154) pour obtenir des distances précises en mètres
centroides_l93 <- st_transform(centroides, 2154)
lycees_publics_l93 <- st_transform(lycees_publics, 2154)

# Recherche de l'index du lycée public le plus proche
nearest_idx <- st_nearest_feature(centroides_l93, lycees_publics_l93)

# Calcul effectif de la distance
centroides_l93$dist_lycee_public_m <- st_distance(
  centroides_l93, 
  lycees_publics_l93[nearest_idx, ], 
  by_element = TRUE
)

# Conversion de la distance en kilomètres pour plus de lisibilité en économétrie
centroides_l93 <- centroides_l93 %>%
  mutate(dist_lycee_public_km = as.numeric(dist_lycee_public_m) / 1000)

