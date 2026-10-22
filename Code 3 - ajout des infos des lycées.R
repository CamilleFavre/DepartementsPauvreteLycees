# =====================================================================
# CODE 3 (adapté) - Variable dépendante (lycées) et contrôles
#
# Entrée  : df_poplarge.parquet (sortie du Code 2 corrigé)
# Sortie  : df_analyse.parquet  (une ligne par commune)
#
# Contenu :
#  A. Comptages de lycées par commune (publics / privés, voie GT / pro,
#     sous contrat / hors contrat) à partir de l'annuaire
#  B. Contrôles géographiques (Code 1 fusionné ici) :
#     distance au lycée public le plus proche, nb de lycées dans un
#     rayon autour du centroïde (hors commune elle-même)
#  C. Variables de contexte scolaire (IPS, effectifs, réussite) :
#     ATTENTION, probablement endogènes (cf. commentaires) -> à ne pas
#     mettre dans le modèle principal
# =====================================================================

library(dplyr)
library(stringr)
library(arrow)
library(sf)

# ---------------------------------------------------------------------
# 0. PARAMÈTRES  (AGREGER_PLM doit être IDENTIQUE à celui du Code 2)
# ---------------------------------------------------------------------
DOSSIER           <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
DOSSIER_OUT       <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
CHEMIN_CENTROIDES <- "C:/Users/camil/Downloads/centroides-communes-population.geojson"
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
# Les indicatrices de l'annuaire peuvent être "1"/"0", 1/0 ou TRUE/FALSE
flag <- function(x) as.character(x) %in% c("1", "TRUE", "true", "Oui", "OUI", "oui")

# ---------------------------------------------------------------------
# 1. LECTURE
# ---------------------------------------------------------------------
df_poplarge   <- read_parquet(file.path(DOSSIER_OUT, "df_poplarge.parquet"))
annuairelycees <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))
reussitelycee  <- read_parquet(file.path(DOSSIER,
                                         "fr-en-indicateurs-de-resultat-des-lycees-denseignement-general-et-technologique.parquet"))
IPSlycee       <- read_parquet(file.path(DOSSIER, "fr-en-ips-lycees-ap2023.parquet"))

# Vérifications de modalités (à lire : adapter les filtres si besoin)
cat("type_etablissement :\n");     print(table(annuairelycees$type_etablissement, useNA = "ifany"))
cat("statut_public_prive :\n");    print(table(annuairelycees$statut_public_prive, useNA = "ifany"))
cat("etat :\n");                   print(table(annuairelycees$etat, useNA = "ifany"))
cat("type_contrat_prive :\n");     print(table(annuairelycees$type_contrat_prive, useNA = "ifany"))

# ---------------------------------------------------------------------
# 2. A. LYCÉES : NETTOYAGE ET CATÉGORIES (mutuellement exclusives)
# ---------------------------------------------------------------------
lycees <- annuairelycees %>%
  filter(
    grepl("^Lyc", type_etablissement),  # Robuste à l'encodage (Lycée / Lycee)
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
    
    # Correction de la détection du contrat privé :
    sous_contrat = secteur == "prive" & (
      grepl("CONTRAT", type_contrat_prive) & !grepl("HORS CONTRAT", type_contrat_prive)
    ),
    hors_contrat = secteur == "prive" & (
      grepl("HORS CONTRAT", type_contrat_prive)
    )
  ) %>%
  filter(!is.na(secteur))

table(lycees$hors_contrat)
table(annuairelycees$statut_public_prive)
table(annuairelycees$type_contrat_prive)

stopifnot(nrow(lycees) > 0)
cat("\nNb de lycées ouverts retenus :", nrow(lycees),
    "| publics :", sum(lycees$secteur == "public"),
    "| privés :",  sum(lycees$secteur == "prive"), "\n")

lycees_agg <- lycees %>%
  group_by(code_insee) %>%
  summarise(
    nb_lycees_publics = sum(secteur == "public"),
    nb_lycees_prives  = sum(secteur == "prive"),
    
    # Les catégories GT / pro-seul s'additionnent au total (pas de double compte)
    nb_publics_gt       = sum(secteur == "public" & gt),
    nb_publics_pro_seul = sum(secteur == "public" & pro_seul),
    nb_prives_gt        = sum(secteur == "prive" & gt),
    nb_prives_pro_seul  = sum(secteur == "prive" & pro_seul),
    
    # Sous contrat (financés par l'État, programmes nationaux) vs hors contrat
    nb_prives_sous_contrat = sum(sous_contrat),
    nb_prives_hors_contrat = sum(hors_contrat),
    .groups = "drop"
  )

# ---------------------------------------------------------------------
# 3. B. CONTRÔLES GÉOGRAPHIQUES (ex-Code 1)
# ---------------------------------------------------------------------
lycees_sf <- lycees %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(2154)

publics    <- lycees_sf %>% filter(secteur == "public")
publics_gt <- lycees_sf %>% filter(secteur == "public", gt)
prives     <- lycees_sf %>% filter(secteur == "prive")

# Centroïdes de POPULATION calculés à partir des carreaux Filosofi
# (fonctionne pour les arrondissements, contrairement au geojson)
DOSSIER_ECO  <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/data_eco"
DEMI_CARREAU <- 100   # m ; mettre 0 si x_l93/y_l93 sont déjà les centres des carreaux

files_eco <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*\\.parquet", full.names = TRUE)

centroides <- bind_rows(lapply(files_eco, function(f) {
  read_parquet(f, col_select = c("lcog_geo", "x_l93", "y_l93", "ind")) %>%
    mutate(code_insee = harmoniser_code(substr(as.character(lcog_geo), 1, 5))) %>%
    filter(code_insee %in% df_poplarge$code_insee, is.finite(ind), ind > 0) %>%
    group_by(code_insee) %>%
    summarise(x = weighted.mean(x_l93 + DEMI_CARREAU, ind),
              y = weighted.mean(y_l93 + DEMI_CARREAU, ind), .groups = "drop")
})) %>%
  st_as_sf(coords = c("x", "y"), crs = 2154)

# Contrôle : doit afficher 0
cat("Communes sans centroïde :", length(setdiff(df_poplarge$code_insee, centroides$code_insee)), "\n")


dist_plus_proche_km <- function(cent, cible) {
  idx <- st_nearest_feature(cent, cible)
  as.numeric(st_distance(cent, cible[idx, ], by_element = TRUE)) / 1000
}

# Nombre de lycées dans le rayon, HORS lycées de la commune elle-même
# (la commune est déjà couverte par nb_lycees_*) : contrôle de concurrence / débordements
nb_dans_rayon_autres <- function(cent, cible) {
  idx <- st_is_within_distance(cent, cible, dist = RAYON_KM * 1000)
  mapply(function(ids, code) sum(cible$code_insee[ids] != code),
         idx, cent$code_insee)
}

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
#    Ces variables sont des RÉSULTATS ENDOGÈNES de la présence du privé
#    (écrémage social des publics) : à utiliser en analyse complémentaire,
#    pas comme contrôles du modèle principal.
# ---------------------------------------------------------------------

# 4a. IPS (indice de position sociale) : dernière rentrée disponible
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

# 4b. Effectifs GT et réussite (dernière année), rattachés à la commune via l'UAI
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
df_analyse <- df_poplarge %>%
  left_join(lycees_agg,    by = "code_insee") %>%
  left_join(controles_geo, by = "code_insee") %>%
  left_join(ips_com,       by = "code_insee") %>%
  left_join(reussite_com,  by = "code_insee") %>%
  # Les NA de comptage = aucun lycée de ce type dans la commune -> 0
  mutate(across(starts_with("nb_lycees_") | starts_with("nb_publics_gt") |
                  starts_with("nb_publics_pro") | starts_with("nb_prives_gt") |
                  starts_with("nb_prives_pro") | starts_with("nb_prives_sous") |
                  starts_with("nb_prives_hors"),
                ~ coalesce(., 0L))) %>%
  mutate(dept = ifelse(startsWith(code_insee, "97"),
                       substr(code_insee, 1, 3), substr(code_insee, 1, 2)))

write_parquet(df_analyse, file.path(DOSSIER_OUT, "df_analyse.parquet"))

# ---------------------------------------------------------------------
# 6. CONTRÔLES DE QUALITÉ DE LA JOINTURE
# ---------------------------------------------------------------------
cat("\n=== CONTRÔLES ===\n")
cat("Communes dans df_analyse :", nrow(df_analyse), "\n")

# Lycées de l'annuaire dont la commune n'est pas dans df_poplarge (normal : petites communes)
# -> mais vérifier que Paris / Lyon / Marseille ne sont pas perdues par un problème de codes
hors <- lycees %>% filter(!code_insee %in% df_poplarge$code_insee)
cat("Lycées situés hors des communes retenues :", nrow(hors), "sur", nrow(lycees), "\n")
cat("Lycées non appariés à Paris/Lyon/Marseille (codes d'arrondissement ou de commune) :",
    sum(grepl("^(751[0-2][0-9]|6938[1-9]|132[0-1][0-9]|75056|69123|13055)$", hors$code_insee)), "\n")
cat("Communes de df_poplarge sans centroïde (distance NA) :",
    sum(is.na(df_analyse$dist_lycee_public_km)), "\n")
print(df_analyse %>% filter(is.na(dist_lycee_public_km)) %>% select(code_insee, nom_commune))

cat("\nDistribution de nb_lycees_prives :\n")
print(table(df_analyse$nb_lycees_prives))
cat(sprintf("Moyenne = %.2f | Variance = %.2f | Part de zéros = %.1f%%\n",
            mean(df_analyse$nb_lycees_prives), var(df_analyse$nb_lycees_prives),
            100 * mean(df_analyse$nb_lycees_prives == 0)))
cat("Communes sans lycée public :", sum(df_analyse$nb_lycees_publics == 0), "\n")