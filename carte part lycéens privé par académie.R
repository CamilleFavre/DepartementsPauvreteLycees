library(arrow); library(dplyr); library(sf); library(ggplot2); library(stringi)
library(stringr)

DOSSIER     <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
CHEMIN_GEO  <- file.path(DOSSIER_OUT, "fr-en-contour-academies-2020.geojson")
ANNEE       <- 2022

flag <- function(x) as.character(x) %in% c("1","TRUE","true","Oui","OUI","oui")

# Normalisation des noms d'académies (accents, casse, préfixe)
norm_acad <- function(x) {
  x <- tolower(stri_trans_general(as.character(x), "Latin-ASCII"))
  x <- trimws(gsub("^academie (de |d')?", "", x))
  gsub("[^a-z0-9]", "", x)   # suppress accents ok, on retire espaces et tirets
}

# ---------------------------------------------------------------------
# 1. STATUT DES LYCÉES
# ---------------------------------------------------------------------
annuaire <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))

lycees_statut <- annuaire %>%
  filter(grepl("^Lyc", type_etablissement),
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

# ---------------------------------------------------------------------
# 2. EFFECTIFS PAR LYCÉE (voie GT, rentrée 2022)
# ---------------------------------------------------------------------

reussite <- read_parquet(file.path(DOSSIER,"fr-en-indicateurs-de-resultat-des-lycees-denseignement-general-et-technologique.parquet"))

cat("Années disponibles :"); print(table(reussite$annee))

eff <- reussite %>%
  filter(annee == ANNEE) %>%
  transmute(
    uai = code_etablissement,
    eff = rowSums(cbind(suppressWarnings(as.numeric(effectif_de_seconde)),
                        suppressWarnings(as.numeric(effectif_de_premiere)),
                        suppressWarnings(as.numeric(effectif_de_terminale))),
                  na.rm = TRUE)) %>%
  group_by(uai) %>% summarise(eff = sum(eff), .groups = "drop")

# ---------------------------------------------------------------------
# 3. PART DU PRIVÉ SOUS CONTRAT PAR ACADÉMIE
# ---------------------------------------------------------------------
part_academie <- eff %>%
  left_join(select(lycees_statut, uai, academie, secteur, sous_contrat), by = "uai") %>%
  filter(!is.na(academie)) %>%
  group_by(academie) %>%
  summarise(eff_total = sum(eff),
            eff_prive = sum(eff[secteur == "prive" & sous_contrat]),
            .groups = "drop") %>%
  mutate(part_prive_contrat = 100 * eff_prive / eff_total)

cat("\nPart nationale (GT) :",
    round(100 * sum(part_academie$eff_prive) / sum(part_academie$eff_total), 1), "%\n")

# ---------------------------------------------------------------------
# 4. CONTOURS DES ACADÉMIES (fichier local, colonne détectée automatiquement)
# ---------------------------------------------------------------------
carte_acad <- st_read(CHEMIN_GEO, quiet = TRUE)
names(carte_acad)   # <- regardez la sortie pour identifier la colonne du nom

# Détection : colonne contenant "acad" dans son nom, sinon contenant "lib" ou "nom"
cand <- names(carte_acad)[grepl("acad|lib|nom", tolower(names(carte_acad)))]
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
library(scales)

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
    name     = "Share of high-scool students in a private institution"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 1,     # titre aligné à droite
    ticks = FALSE, frame.colour = NA
  )) +
  theme_void(base_family = "sans") +
  theme(
    # --- légende AU-DESSUS de la carte, alignée à droite ---
    legend.position      = "top",
    legend.justification = "right",
    legend.direction     = "horizontal",
    legend.key.width     = unit(1, "cm"),        # longueur de la barre
    legend.key.height    = unit(0.25, "cm"),     # PETITE barre fine
    legend.title         = element_text(size = 9, colour = "grey20"),
    legend.text          = element_text(size = 8, colour = "grey30"),
    legend.margin        = margin(b = 6),
    plot.margin          = margin(8, 8, 8, 8),
    plot.background      = element_rect(fill = "white", colour = NA)
  )

print(p)
