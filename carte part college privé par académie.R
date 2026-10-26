# =====================================================================
# CARTE — PART DES COLLÉGIENS DANS UN COLLÈGE PRIVÉ PAR ACADÉMIE
# Source : annuaire éducation (statut) + IVP collèges (effectifs)
# =====================================================================
library(arrow); library(dplyr); library(sf); library(ggplot2); library(stringi)
library(stringr); library(scales)

# ---------------------------------------------------------------------
# 0. PARAMÈTRES — détection du système d'exploitation
# ---------------------------------------------------------------------
os_name <- Sys.info()[["sysname"]]

if (os_name == "Windows") {
  DOSSIER           <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT       <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
  CHEMIN_GEO        <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/fr-en-contour-academies-2020.geojson"
  
} else if (os_name == "Darwin") { # macOS
  DOSSIER           <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT       <- "/Users/camillefavre/DepartementsPauvreteLycees"
  CHEMIN_GEO        <- "/Users/camillefavre/DepartementsPauvreteLycees/fr-en-contour-academies-2020.geojson" # Ajustez si besoin
  
} else {
  stop("Système d'exploitation non pris en charge.")
}

FICHIER_EFF_COL <- "fr-en-indicateurs-valeur-ajoutee-colleges.parquet"

# Normalisation des noms d'académies (accents, casse, préfixe)
norm_acad <- function(x) {
  x <- tolower(stri_trans_general(as.character(x), "Latin-ASCII"))
  x <- trimws(gsub("^academie (de |d')?", "", x))
  gsub("[^a-z0-9]", "", x)
}

# ---------------------------------------------------------------------
# 1. STATUT DES COLLÈGES
# ---------------------------------------------------------------------
annuaire <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))

colleges_statut <- annuaire %>%
  filter(grepl("^Coll", type_etablissement),
         toupper(coalesce(as.character(etat), "")) == "OUVERT",
         !str_detect(code_postal, "^973")) %>%
  transmute(
    uai       = identifiant_de_l_etablissement,
    academie  = norm_acad(libelle_academie),
    secteur   = case_when(
      statut_public_prive == "Public"     ~ "public",
      grepl("^Priv", statut_public_prive) ~ "prive",
      TRUE                                ~ NA_character_),
    sous_contrat = grepl("CONTRAT", type_contrat_prive) &
      !grepl("HORS CONTRAT", type_contrat_prive)
  ) %>%
  filter(!is.na(secteur))

cat("Collèges ouverts :", nrow(colleges_statut),
    "| publics :", sum(colleges_statut$secteur == "public"),
    "| privés :",  sum(colleges_statut$secteur == "prive"), "\n")

# ---------------------------------------------------------------------
# 2. EFFECTIFS PAR COLLÈGE (proxy : candidats au brevet, dernière session)
# ---------------------------------------------------------------------
eff_brut <- read_parquet(file.path(DOSSIER, FICHIER_EFF_COL))

cat("Sessions disponibles :\n"); print(sort(unique(eff_brut$session)))

# session est une Date -> extraire l'année
eff_brut <- eff_brut %>%
  mutate(annee_session = as.numeric(format(as.Date(session), "%Y")))

derniere <- max(eff_brut$annee_session, na.rm = TRUE)
cat("Session retenue :", derniere, "\n")

eff <- eff_brut %>%
  filter(annee_session == derniere) %>%
  transmute(
    uai      = as.character(uai),
    academie = norm_acad(academie),
    secteur  = ifelse(toupper(as.character(secteur)) == "PR", "prive", "public"),
    eff      = suppressWarnings(as.numeric(nb_candidats_g)) +
      coalesce(suppressWarnings(as.numeric(nb_candidats_p)), 0)   # <<< ici
  ) %>%
  filter(is.finite(eff), !is.na(academie))

cat("Collèges avec candidats :", nrow(eff),
    "| somme des candidats :", sum(eff$eff), "\n")
cat("Répartition secteurs :\n"); print(table(eff$secteur))
# Contrôles attendus : ~750 000 candidats (une classe d'âge de 3e),
# part nationale ~20-22 % (imprimée en section 3).

# ---------------------------------------------------------------------
# 3. PART DU PRIVÉ PAR ACADÉMIE (base : candidats au brevet, proxy des 3e)
# ---------------------------------------------------------------------
part_academie <- eff %>%
  group_by(academie) %>%
  summarise(eff_total = sum(eff),
            eff_prive = sum(eff[secteur == "prive"]),
            .groups = "drop") %>%
  mutate(part_prive_contrat = 100 * eff_prive / eff_total)

cat("\nPart nationale (collégiens 3e, proxy candidats brevet) :",
    round(100 * sum(part_academie$eff_prive) / sum(part_academie$eff_total), 1), "%\n")
# Contrôle attendu : ~20-22 % (MEN).

print(part_academie %>% arrange(desc(part_prive_contrat)), n = 30)

# ---------------------------------------------------------------------
# 4. CONTOURS DES ACADÉMIES (fichier local, colonne détectée automatiquement)
# ---------------------------------------------------------------------
carte_acad <- st_read(CHEMIN_GEO, quiet = TRUE)
cand  <- names(carte_acad)[grepl("acad|lib|nom", tolower(names(carte_acad)))]
col_acad <- cand[1]
if (is.na(col_acad)) stop("Colonne du nom d'académie introuvable : voir names(carte_acad)")
cat("Colonne utilisée pour le nom d'académie :", col_acad, "\n")

carte <- carte_acad %>%
  mutate(academie = norm_acad(.data[[col_acad]])) %>%
  filter(!academie %in% c("guyane", "guadeloupe", "mayotte", "lareunion", "martinique")) %>%
  select(academie, geometry) %>%
  left_join(part_academie, by = "academie")

# Diagnostic de fusion
cat("Académies sans effectifs (verifier le mapping des noms) :",
    sum(is.na(carte$part_prive_contrat)), "sur", nrow(carte), "\n")
if (sum(is.na(carte$part_prive_contrat)) > 0) {
  print(anti_join(carte %>% st_drop_geometry(), part_academie, by = "academie") %>% select(academie))
  print(anti_join(part_academie, carte %>% st_drop_geometry(), by = "academie") %>% select(academie))
}

# ---------------------------------------------------------------------
# 5. CARTE (version publication)
# ---------------------------------------------------------------------
carte_pub <- st_transform(carte, 2154)   # Lambert-93

bornes <- c(floor(min(carte_pub$part_prive_contrat, na.rm = TRUE) / 5) * 5,
            ceiling(max(carte_pub$part_prive_contrat, na.rm = TRUE) / 5) * 5)

p <- ggplot(carte_pub) +
  geom_sf(aes(fill = part_prive_contrat), color = "white", linewidth = 0.2) +
  coord_sf(xlim = c(80000, 1100000), ylim = c(6040000, 7120000),
           expand = FALSE, datum = NA) +
  scale_fill_gradientn(
    colours  = c("#fbf3e4", "#f5cf9b", "#e8896a", "#c2465a", "#6e1d4f", "#2b0f3a"),
    limits   = bornes,
    breaks   = seq(bornes[1], bornes[2], by = 10),
    labels   = label_number(suffix = " %"),
    na.value = "grey90",
    name     = "Share of collège students in a private institution"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 1,
    ticks = FALSE, frame.colour = NA
  )) +
  theme_void(base_family = "sans") +
  theme(
    legend.position      = "top",
    legend.justification = "right",
    legend.direction     = "horizontal",
    legend.key.width     = unit(1, "cm"),
    legend.key.height    = unit(0.25, "cm"),
    legend.title         = element_text(size = 9, colour = "grey20"),
    legend.text          = element_text(size = 8, colour = "grey30"),
    legend.margin        = margin(b = 6),
    plot.margin          = margin(8, 8, 8, 8),
    plot.background      = element_rect(fill = "white", colour = NA)
  )

print(p)
ggsave(file.path(DOSSIER_OUT, "carte_part_prive_colleges.png"),
       p, width = 9, height = 7.5, dpi = 300)
cat("\nExport :", file.path(DOSSIER_OUT, "carte_part_prive_colleges.png"), "\n")