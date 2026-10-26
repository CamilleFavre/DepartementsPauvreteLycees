# =====================================================================
# CODE 3b - COLLÈGES : variable dépendante et contrôles
#
# Entrée  : df_poplarge.parquet (Code 2, inchangé)
# Sortie  : df_analyse_colleges.parquet
#
# Différences avec la version lycées :
#  - filtre "Collège" dans l'annuaire, PAS de notion de voie (GT/pro)
#  - rayon concurrentiel réduit à 5 km (zone de recrutement collège ~3-4 km)
#  - effectifs/IPS via les fichiers collèges (IPS collèges AP2023,
#    résultats au brevet) -> à adapter aux noms de vos fichiers
# =====================================================================

library(dplyr); library(stringr); library(arrow); library(sf)

# Détection automatique du système d'exploitation
os_name <- Sys.info()[["sysname"]]

# Attribution des chemins de base selon l'OS
if (os_name == "Windows") {
  DOSSIER     <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
  
} else if (os_name == "Darwin") { # macOS
  DOSSIER     <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
  
} else {
  stop("Système d'exploitation non pris en charge.")
}

# Dérivation automatique du dossier ECO (gère les barres de manière propre selon l'OS)
DOSSIER_ECO <- file.path(DOSSIER_OUT, "data_eco")

AGREGER_PLM  <- FALSE
RAYON_KM     <- 5          # <- collège : marché local
DEMI_CARREAU <- 100

map_plm <- function(code) dplyr::case_when(
  grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
  grepl("^6938[1-9]$", code)             ~ "69123",
  grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
  TRUE ~ code)
harmoniser_code <- function(x) {
  code <- str_pad(as.character(x), width = 5, pad = "0")
  if (AGREGER_PLM) map_plm(code) else code
}
flag <- function(x) as.character(x) %in% c("1","TRUE","true","Oui","OUI","oui")

# ---------------------------------------------------------------------
# 1. LECTURE
# ---------------------------------------------------------------------
df_poplarge    <- read_parquet(file.path(DOSSIER_OUT, "df_poplarge.parquet"))
annuaire       <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))
IPScoll        <- read_parquet(file.path(DOSSIER, "fr-en-ips-colleges-ap2023.parquet"))
# Résultats au brevet (effectifs 3e) : ADAPTER au nom de votre fichier, ex. :
# brevet <- read_parquet(file.path(DOSSIER, "fr-en-indicateurs-de-resultat-des-colleges...parquet"))
brevet <- NULL   # <- mettre NULL tant que le fichier n'est pas identifié

# ---------------------------------------------------------------------
# 2. COLLÈGES : nettoyage
# ---------------------------------------------------------------------
colleges <- annuaire %>%
  filter(
    grepl("^Coll", type_etablissement),        # <- le seul vrai changement de filtre
    toupper(coalesce(as.character(etat), "")) == "OUVERT"
  ) %>%
  mutate(
    code_insee   = harmoniser_code(code_commune),
    secteur      = case_when(
      statut_public_prive == "Public"     ~ "public",
      grepl("^Priv", statut_public_prive) ~ "prive",
      TRUE                                ~ NA_character_),
    sous_contrat = grepl("CONTRAT", type_contrat_prive) &
      !grepl("HORS CONTRAT", type_contrat_prive),
    hors_contrat = grepl("HORS CONTRAT", type_contrat_prive)
  ) %>%
  filter(!is.na(secteur))

stopifnot(nrow(colleges) > 0)
cat("Collèges ouverts :", nrow(colleges),
    "| publics :", sum(colleges$secteur == "public"),
    "| privés :",  sum(colleges$secteur == "prive"), "\n")

colleges_agg <- colleges %>%
  group_by(code_insee) %>%
  summarise(
    nb_colleges_publics = sum(secteur == "public"),
    nb_colleges_prives  = sum(secteur == "prive"),
    nb_prives_sous_contrat = sum(sous_contrat),
    nb_prives_hors_contrat = sum(hors_contrat),
    .groups = "drop")

# ---------------------------------------------------------------------
# 3. CONTRÔLES GÉOGRAPHIQUES (rayon 5 km)
# ---------------------------------------------------------------------
colleges_sf <- colleges %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(2154)

publics_c <- colleges_sf %>% filter(secteur == "public")
prives_c  <- colleges_sf %>% filter(secteur == "prive")

files_eco <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*\\.parquet", full.names = TRUE)
centroides <- bind_rows(lapply(files_eco, function(f) {
  read_parquet(f, col_select = c("lcog_geo","x_l93","y_l93","ind")) %>%
    mutate(code_insee = harmoniser_code(substr(as.character(lcog_geo), 1, 5))) %>%
    filter(code_insee %in% df_poplarge$code_insee, is.finite(ind), ind > 0) %>%
    group_by(code_insee) %>%
    summarise(x = weighted.mean(x_l93 + DEMI_CARREAU, ind),
              y = weighted.mean(y_l93 + DEMI_CARREAU, ind), .groups = "drop")
})) %>% st_as_sf(coords = c("x","y"), crs = 2154)

dist_plus_proche_km <- function(cent, cible) {
  idx <- st_nearest_feature(cent, cible)
  as.numeric(st_distance(cent, cible[idx, ], by_element = TRUE)) / 1000
}
nb_dans_rayon_autres <- function(cent, cible) {
  idx <- st_is_within_distance(cent, cible, dist = RAYON_KM * 1000)
  mapply(function(ids, code) sum(cible$code_insee[ids] != code),
         idx, cent$code_insee)
}

controles_geo <- centroides %>%
  mutate(
    dist_college_public_km   = dist_plus_proche_km(., publics_c),
    nb_prives_5km_autres     = nb_dans_rayon_autres(., prives_c),
    nb_publics_5km_autres    = nb_dans_rayon_autres(., publics_c)
  ) %>% st_drop_geometry() %>%
  select(code_insee, dist_college_public_km, nb_prives_5km_autres, nb_publics_5km_autres) %>%
  distinct(code_insee, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 4. IPS COLLÈGES (contexte, endogène -> analyse complémentaire)
# ---------------------------------------------------------------------
ips_com <- IPScoll %>%
  filter(rentree_scolaire == max(rentree_scolaire, na.rm = TRUE)) %>%
  mutate(code_insee = harmoniser_code(code_insee_de_la_commune),
         ips = suppressWarnings(as.numeric(ips)),   # <- déjà numérique en fait
         sect = ifelse(grepl("^pub", tolower(secteur)), "public", "prive")) %>%
  filter(is.finite(ips)) %>%
  group_by(code_insee) %>%
  summarise(ips_pub_moy   = mean(ips[sect == "public"], na.rm = TRUE),
            ips_prive_moy = mean(ips[sect == "prive"],  na.rm = TRUE),
            ips_pub_ecart = ifelse(sum(sect == "public") >= 2,
                                   max(ips[sect=="public"]) - min(ips[sect=="public"]), NA_real_),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

# Effectifs 3e (brevet) : à activer quand le fichier est identifié
# eff_com <- brevet %>% ... (même logique que reussite_com côté lycées)

# ---------------------------------------------------------------------
# 5. ASSEMBLAGE
# ---------------------------------------------------------------------
df_analyse_col <- df_poplarge %>%
  left_join(colleges_agg,   by = "code_insee") %>%
  left_join(controles_geo, by = "code_insee") %>%
  left_join(ips_com,        by = "code_insee") %>%
  mutate(across(starts_with("nb_colleges_") | starts_with("nb_prives_sous") |
                  starts_with("nb_prives_hors"),
                ~ coalesce(., 0L))) %>%
  mutate(dept = ifelse(startsWith(code_insee, "97"),
                       substr(code_insee, 1, 3), substr(code_insee, 1, 2)))

write_parquet(df_analyse_col, file.path(DOSSIER_OUT, "df_analyse_colleges.parquet"))

# ---------------------------------------------------------------------
# 6. CONTRÔLES
# ---------------------------------------------------------------------
cat("\n=== CONTRÔLES ===\n")
cat("N =", nrow(df_analyse_col), "\n")
cat("Distribution de nb_colleges_prives :\n"); print(table(df_analyse_col$nb_colleges_prives))
cat(sprintf("Moyenne = %.2f | Variance = %.2f | Zéros = %.1f%%\n",
            mean(df_analyse_col$nb_colleges_prives),
            var(df_analyse_col$nb_colleges_prives),
            100 * mean(df_analyse_col$nb_colleges_prives == 0)))
cat("Sans collège public :", sum(df_analyse_col$nb_colleges_publics == 0), "\n")