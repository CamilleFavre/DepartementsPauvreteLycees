# =====================================================================
# CODE 3 (adapté - TOUTES COMMUNES) - Variable dépendante (lycées) et contrôles
#
# Entrée  : df_final.parquet     (sortie de l'étape 5 du Code 2, toutes communes)
# Sortie  : df_analyse_complet.parquet (une ligne par commune, ~35 000 communes)
#
# Contenu :
#  A. Comptages de lycées par commune (publics / privés, voie GT / pro,
#     sous contrat / hors contrat) à partir de l'annuaire
#  B. Contrôles géographiques :
#     distance au lycée public le plus proche, nb de lycées dans un
#     rayon autour du centroïde (hors commune elle-même)
#  C. Variables de contexte scolaire (IPS, effectifs, réussite) :
#     ATTENTION, probablement endogènes -> à ne pas mettre dans le modèle principal
# =====================================================================

library(dplyr)
library(stringr)
library(arrow)
library(sf)

# ---------------------------------------------------------------------
# 0. PARAMÈTRES (AGREGER_PLM doit être IDENTIQUE à celui du Code 2)
# ---------------------------------------------------------------------
# Détection automatique du système d'exploitation
os_name <- Sys.info()[["sysname"]]

# Attribution des chemins selon l'OS
if (os_name == "Windows") {
  DOSSIER           <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT       <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
  CHEMIN_CENTROIDES <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/centroides-communes-population.geojson"
  
} else if (os_name == "Darwin") { # macOS
  DOSSIER           <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT       <- "/Users/camillefavre/DepartementsPauvreteLycees"
  CHEMIN_CENTROIDES <- "/Users/camillefavre/DepartementsPauvreteLycees/centroides-communes-population.geojson" # Ajustez si besoin
  
} else {
  stop("Système d'exploitation non pris en charge.")
}

AGREGER_PLM       <- FALSE
RAYON_KM          <- 10

map_plm <- function(code) {
  dplyr::case_when(
    grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
    grepl("^6938[1-9]$", code)             ~ "69123",
    grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
    TRUE ~ code
  )
}

harmoniser_code <- function(x) {
  code <- str_pad(as.character(x), width = 5, pad = "0")
  if (AGREGER_PLM) map_plm(code) else code
}

flag <- function(x) as.character(x) %in% c("1", "TRUE", "true", "Oui", "OUI", "oui")

# ---------------------------------------------------------------------
# 1. LECTURE (df_final au lieu de df_poplarge)
# ---------------------------------------------------------------------
message("Lecture du référentiel complet des communes (df_final.parquet)...")
df_toutes_communes <- read_parquet(file.path(DOSSIER_OUT, "df_final.parquet"))

annuairelycees <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))
reussitelycee  <- read_parquet(file.path(DOSSIER, "fr-en-indicateurs-de-resultat-des-lycees-denseignement-general-et-technologique.parquet"))
IPSlycee       <- read_parquet(file.path(DOSSIER, "fr-en-ips-lycees-ap2023.parquet"))

# Vérifications de modalités
cat("type_etablissement :\n");     print(table(annuairelycees$type_etablissement, useNA = "ifany"))
cat("statut_public_prive :\n");    print(table(annuairelycees$statut_public_prive, useNA = "ifany"))
cat("etat :\n");                   print(table(annuairelycees$etat, useNA = "ifany"))
cat("type_contrat_prive :\n");     print(table(annuairelycees$type_contrat_prive, useNA = "ifany"))

# ---------------------------------------------------------------------
# 2. A. LYCÉES : NETTOYAGE ET CATÉGORIES (mutuellement exclusives)
# ---------------------------------------------------------------------
lycees <- annuairelycees %>%
  filter(
    grepl("^Lyc", type_etablissement),
    toupper(coalesce(as.character(etat), "")) == "OUVERT"
  ) %>%
  mutate(
    code_insee = harmoniser_code(code_commune),
    secteur = case_when(
      statut_public_prive == "Public"        ~ "public",
      grepl("^Priv", statut_public_prive)    ~ "prive",
      TRUE                                   ~ NA_character_
    ),
    gt       = flag(voie_generale) | flag(voie_technologique),
    pro_seul = flag(voie_professionnelle) & !gt,
    
    sous_contrat = secteur == "prive" & (
      grepl("CONTRAT", type_contrat_prive) & !grepl("HORS CONTRAT", type_contrat_prive)
    ),
    hors_contrat = secteur == "prive" & (
      grepl("HORS CONTRAT", type_contrat_prive)
    )
  ) %>%
  filter(!is.na(secteur))

stopifnot(nrow(lycees) > 0)
cat("\nNb de lycées ouverts retenus :", nrow(lycees),
    "| publics :", sum(lycees$secteur == "public"),
    "| privés :",  sum(lycees$secteur == "prive"), "\n")

lycees_agg <- lycees %>%
  group_by(code_insee) %>%
  summarise(
    nb_lycees_publics = sum(secteur == "public"),
    nb_lycees_prives  = sum(secteur == "prive"),
    
    nb_publics_gt       = sum(secteur == "public" & gt),
    nb_publics_pro_seul = sum(secteur == "public" & pro_seul),
    nb_prives_gt        = sum(secteur == "prive" & gt),
    nb_prives_pro_seul  = sum(secteur == "prive" & pro_seul),
    
    nb_prives_sous_contrat = sum(sous_contrat),
    nb_prives_hors_contrat = sum(hors_contrat),
    .groups = "drop"
  )

# ---------------------------------------------------------------------
# 3. B. CONTRÔLES GÉOGRAPHIQUES (pour toutes les communes)
# ---------------------------------------------------------------------
lycees_sf <- lycees %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(2154)

publics    <- lycees_sf %>% filter(secteur == "public")
publics_gt <- lycees_sf %>% filter(secteur == "public", gt)
prives     <- lycees_sf %>% filter(secteur == "prive")

# Centroïdes de population calculés depuis les carreaux Filosofi
DOSSIER_ECO  <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/data_eco"
DEMI_CARREAU <- 100   # m

files_eco <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*\\.parquet", full.names = TRUE)

message("Calcul des centroïdes pour toutes les communes...")
centroides <- bind_rows(lapply(files_eco, function(f) {
  read_parquet(f, col_select = c("lcog_geo", "x_l93", "y_l93", "ind")) %>%
    mutate(code_insee = harmoniser_code(substr(as.character(lcog_geo), 1, 5))) %>%
    filter(code_insee %in% df_toutes_communes$code_insee, is.finite(ind), ind > 0) %>%
    group_by(code_insee) %>%
    summarise(x = weighted.mean(x_l93 + DEMI_CARREAU, ind),
              y = weighted.mean(y_l93 + DEMI_CARREAU, ind), .groups = "drop")
})) %>%
  st_as_sf(coords = c("x", "y"), crs = 2154)

cat("Communes sans centroïde dans les carreaux :", 
    length(setdiff(df_toutes_communes$code_insee, centroides$code_insee)), "\n")

dist_plus_proche_km <- function(cent, cible) {
  idx <- st_nearest_feature(cent, cible)
  as.numeric(st_distance(cent, cible[idx, ], by_element = TRUE)) / 1000
}

nb_dans_rayon_autres <- function(cent, cible) {
  idx <- st_is_within_distance(cent, cible, dist = RAYON_KM * 1000)
  mapply(function(ids, code) sum(cible$code_insee[ids] != code),
         idx, cent$code_insee)
}

message("Calcul des distances et de l'environnement concurrentiel...")
controles_geo <- centroides %>%
  mutate(
    dist_lycee_public_km    = dist_plus_proche_km(., publics),
    dist_lycee_public_gt_km = dist_plus_proche_km(., publics_gt),
    nb_prives_10km_autres   = nb_dans_rayon_autres(., prives),
    nb_publics_10km_autres  = nb_dans_rayon_autres(., publics)
  ) %>%
  st_drop_geometry() %>%
  select(code_insee, dist_lycee_public_km, dist_lycee_public_gt_km,
         nb_prives_10km_autres, nb_publics_10km_autres) %>%
  distinct(code_insee, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 4. C. CONTEXTE SCOLAIRE (IPS, effectifs, réussite)
# ---------------------------------------------------------------------

# 4a. IPS
ips_com <- IPSlycee %>%
  filter(rentree_scolaire == max(rentree_scolaire, na.rm = TRUE)) %>%
  mutate(
    code_insee = harmoniser_code(code_insee_de_la_commune),
    ips = suppressWarnings(as.numeric(ips_etab)),
    sect = ifelse(grepl("^pub", tolower(secteur)), "public", "prive")
  ) %>%
  filter(is.finite(ips)) %>%
  group_by(code_insee) %>%
  summarise(
    ips_pub_moy   = mean(ips[sect == "public"], na.rm = TRUE),
    ips_pub_ecart = ifelse(sum(sect == "public") >= 2,
                           max(ips[sect == "public"]) - min(ips[sect == "public"]), NA_real_),
    ips_prive_moy = mean(ips[sect == "prive"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

# 4b. Effectifs GT et réussite
code_par_uai <- annuairelycees %>%
  transmute(uai = identifiant_de_l_etablissement,
            code_insee = harmoniser_code(code_commune)) %>%
  distinct(uai, .keep_all = TRUE)

reussite_com <- reussitelycee %>%
  filter(annee == max(annee, na.rm = TRUE)) %>%
  mutate(
    uai  = code_etablissement,
    sect = ifelse(grepl("^pu", tolower(secteur_public_pu_prive_pr)), "public", "prive"),
    eff  = rowSums(cbind(suppressWarnings(as.numeric(effectif_de_seconde)),
                         suppressWarnings(as.numeric(effectif_de_premiere)),
                         suppressWarnings(as.numeric(effectif_de_terminale))), na.rm = TRUE),
    va   = suppressWarnings(as.numeric(va_reu_total))
  ) %>%
  left_join(code_par_uai, by = "uai") %>%
  filter(!is.na(code_insee)) %>%
  group_by(code_insee) %>%
  summarise(
    eff_gt_public = sum(eff[sect == "public"], na.rm = TRUE),
    eff_gt_prive  = sum(eff[sect == "prive"],  na.rm = TRUE),
    va_reu_pub_moy = mean(va[sect == "public"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    part_eff_prive_gt = ifelse(eff_gt_public + eff_gt_prive > 0,
                               eff_gt_prive / (eff_gt_public + eff_gt_prive), NA_real_),
    va_reu_pub_moy = ifelse(is.nan(va_reu_pub_moy), NA_real_, va_reu_pub_moy)
  )

# ---------------------------------------------------------------------
# 5. ASSEMBLAGE
# ---------------------------------------------------------------------
message("Assemblage final du data frame...")
df_analyse_all <- df_toutes_communes %>%
  left_join(lycees_agg,    by = "code_insee") %>%
  left_join(controles_geo, by = "code_insee") %>%
  left_join(ips_com,       by = "code_insee") %>%
  left_join(reussite_com,  by = "code_insee") %>%
  # Les NA de comptage pour les communes sans lycée de tel type -> 0
  mutate(across(starts_with("nb_lycees_") | starts_with("nb_publics_gt") |
                  starts_with("nb_publics_pro") | starts_with("nb_prives_gt") |
                  starts_with("nb_prives_pro") | starts_with("nb_prives_sous") |
                  starts_with("nb_prives_hors"),
                ~ coalesce(., 0L))) %>%
  mutate(dept = ifelse(startsWith(code_insee, "97"),
                       substr(code_insee, 1, 3), substr(code_insee, 1, 2)))

# Sauvegarde
write_parquet(df_analyse_all, file.path(DOSSIER_OUT, "df_analyse_all.parquet"))
message("Sauvegarde effectuée dans : ", file.path(DOSSIER_OUT, "df_analyse.parquet"))

# ---------------------------------------------------------------------
# 6. CONTRÔLES DE QUALITÉ
# ---------------------------------------------------------------------
cat("\n=== CONTRÔLES QUALITÉ ===\n")
cat("Communes dans df_analyse :", nrow(df_analyse), "\n")

hors <- lycees %>% filter(!code_insee %in% df_toutes_communes$code_insee)
cat("Lycées situés hors du référentiel des communes :", nrow(hors), "sur", nrow(lycees), "\n")

cat("Communes sans centroïde (distance NA) :",
    sum(is.na(df_analyse$dist_lycee_public_km)), "\n")

cat("\nDistribution de nb_lycees_prives sur TOUTES les communes :\n")
print(table(df_analyse$nb_lycees_prives))
cat(sprintf("Moyenne = %.2f | Variance = %.2f | Part de zéros = %.1f%%\n",
            mean(df_analyse$nb_lycees_prives), var(df_analyse$nb_lycees_prives),
            100 * mean(df_analyse$nb_lycees_prives == 0)))