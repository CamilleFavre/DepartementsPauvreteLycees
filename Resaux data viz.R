# =====================================================================
# CODE 2c - EXPORT JSON : RESEAUX DE CARREAUX PAR VILLE
#
# Pour chaque ville cible, exporte :
#   - les noeuds : un carreau Filosofi 200 m (coordonnees, revenu,
#     taux de pauvrete, population, menages, menages pauvres)
#   - les liens : paires de carreaux voisins a < 285 m (voisinage
#     intra-communal, identique au Code 2 v2 ; recours kNN pour les
#     carreaux isoles afin que le reseau reste connexe)
#   - un resume communal (population, revenu moyen, taux de pauvrete)
#
# Le fichier reseaux_villes.json est PUISé dans l'application React
# "reseau-quartiers" (bouton Importer > coller le contenu du JSON).
# Les quadrants LISA (HH/LL/HL/LH), z-scores et contrastes sont
# calcules dans l'application elle-meme, a partir des liens.
# =====================================================================

library(arrow)
library(dplyr)
library(jsonlite)
library(spdep)

# ---------------------------------------------------------------------
# 0. PARAMETRES
# ---------------------------------------------------------------------
os_name <- Sys.info()[["sysname"]]
if (os_name == "Windows") {
  DOSSIER_ECO <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/data_eco"
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
} else if (os_name == "Darwin") {
  DOSSIER_ECO <- "/Users/camillefavre/DepartementsPauvreteLycees/data_eco"
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
} else {
  stop("Systeme d'exploitation non pris en charge.")
}

DEPTS_CIBLES <- c("35")  # NULL = tous les fichiers
VILLES        <- NULL   # codes INSEE precis, ex. c("69123","33063") ; NULL = top N
TOP_N         <- 8      # nb de villes (les plus peuplees) si VILLES = NULL
MIN_POP       <- 20000  # plancher de population pour les villes retenues
DIST_VOISIN   <- 285
K_KNN         <- 4
MIN_CARREAUX  <- 10

FICHIER_JSON <- file.path(DOSSIER_OUT, "reseaux_villes.json")

# ---------------------------------------------------------------------
# 1. EXTRACTION D'UNE VILLE -> liste   JSON (noeuds + liens + resume)
# ---------------------------------------------------------------------
extraire_ville <- function(d) {
  d <- d[is.finite(d$revenu_moyen) & is.finite(d$taux_pauvrete) &
           !is.na(d$x_l93) & !is.na(d$y_l93) & d$ind > 0, ]
  n <- nrow(d)
  if (n < MIN_CARREAUX) return(NULL)
  
  coords <- cbind(d$x_l93, d$y_l93)
  nb <- dnearneigh(coords, d1 = 0, d2 = DIST_VOISIN)
  card0 <- card(nb)
  if (any(card0 == 0)) {                      # recours kNN pour les isoles
    kk <- min(K_KNN, n - 1)
    knnb <- as.list(knn2nb(knearneigh(coords, k = kk)))
    nbl <- as.list(nb)
    for (j in which(card0 == 0)) nbl[[j]] <- knnb[[j]]
    nb <- structure(nbl, class = "nb",
                    region.id = attr(nb, "region.id"), sym = attr(nb, "sym"))
  }
  
  # Liens (paires uniques, indices 0-based pour le JS)
  liens <- lapply(seq_len(n), function(i) {
    v <- nb[[i]]; v <- v[v > i]
    if (length(v) == 0) NULL else cbind(i - 1L, v - 1L)
  })
  liens <- matrix(do.call(rbind, liens), ncol = 2)
  
  # Noeuds
  noeuds <- data.frame(
    i    = seq_len(n) - 1L,
    x    = round(d$x_l93), y = round(d$y_l93),
    rev  = round(d$revenu_moyen),
    pauv = round(d$taux_pauvrete, 1),
    ind  = as.integer(d$ind),
    men  = as.integer(d$men),
    menPauv = as.integer(d$men_pauv)
  )
  
  list(
    code     = as.character(d$code_insee[1]),
    nom      = as.character(d$nom_ville[1]),
    pop      = as.integer(sum(d$ind)),
    revMoyen = round(weighted.mean(d$revenu_moyen, d$ind, na.rm = TRUE)),
    tauxPauv = round(sum(d$men_pauv, na.rm = TRUE) / sum(d$men, na.rm = TRUE) * 100, 1),
    nCarreaux = n,
    nLiens    = nrow(liens),
    nodes = noeuds,
    edges = liens
  )
}

# ---------------------------------------------------------------------
# 2. BOUCLE : identifier les villes cibles et extraire
# ---------------------------------------------------------------------
files <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*[.]parquet$", full.names = TRUE)
if (!is.null(DEPTS_CIBLES))
  files <- files[substr(basename(files), 10, 11) %in% DEPTS_CIBLES]
cat("Departements traites :", length(files), "\n")

toutes_villes <- c()
if (!is.null(VILLES)) toutes_villes <- VILLES

resultats <- list()
for (f in files) {
  df <- read_parquet(f) %>%
    mutate(code_insee = substr(as.character(lcog_geo), 1, 5),
           taux_pauvrete = ifelse(!is.na(men) & men > 0, men_pauv / men * 100, NA_real_),
           nom_ville = code_insee)   # nom remplace ensuite via l'API si dispo
  
  # pop par commune du fichier
  pop_com <- df %>% group_by(code_insee) %>%
    summarise(pop = sum(ind, na.rm = TRUE), .groups = "drop") %>%
    filter(pop >= MIN_POP) %>%
    arrange(desc(pop))
  
  # si pas de liste manuelle : top N des dept traites
  if (is.null(VILLES)) {
    candidats <- pop_com$code_insee
    if (length(toutes_villes) + length(candidats) > TOP_N)
      candidats <- candidats[seq_len(max(TOP_N - length(toutes_villes), 0))]
    toutes_villes <- c(toutes_villes, candidats)
  }
  
  cibles <- intersect(pop_com$code_insee, toutes_villes)
  if (length(cibles) == 0) next
  
  for (cv in cibles) {
    res <- extraire_ville(df[df$code_insee == cv, ])
    if (!is.null(res)) {
      resultats[[res$code]] <- res
      cat(sprintf("  %s (%s) : %d carreaux, %d liens\n",
                  res$code, basename(f), res$nCarreaux, res$nLiens))
    }
  }
}

if (length(resultats) == 0) stop("Aucune ville extraite : verifier VILLES / TOP_N / MIN_POP.")

# ---------------------------------------------------------------------
# 3. NOMS DE COMMUNES (API geo.api.gouv.fr) + EXPORT
# ---------------------------------------------------------------------
ref <- tryCatch(
  fromJSON("https://geo.api.gouv.fr/communes?fields=nom,code&format=json"),
  error = function(e) NULL)
if (!is.null(ref)) {
  for (k in names(resultats)) {
    m <- match(k, ref$code)
    if (!is.na(m)) resultats[[k]]$nom <- ref$nom[m]
  }
}

json <- toJSON(resultats, auto_unbox = TRUE, dataframe = "rows", digits = NA)
writeLines(json, FICHIER_JSON)
cat("\nExport :", FICHIER_JSON, "\n")
cat("Taille :", round(file.size(FICHIER_JSON) / 1024), "Ko |",
    length(resultats), "villes\n")
cat(">>> Ouvrir reseaux_villes.json, copier tout le contenu, et le coller\n",
    "    dans l'application (bouton Importer) pour remplacer la ville demo.\n")












library(jsonlite)

DOSSIER <- "/Users/camillefavre/DepartementsPauvreteLycees"
app <- file.path(DOSSIER, "réseau-des-quartiers-version-sombre.html")

# 1. JSON de Rennes (code INSEE 35238) : soit tout le fichier, soit rennes.json
j <- fromJSON(file.path(DOSSIER, "reseaux_villes.json"), simplifyVector = FALSE)
rennes <- Filter(function(v) v$code == "35238", j)[[1]]
stopifnot(!is.null(rennes$nodes), !is.null(rennes$edges))
json_rennes <- toJSON(rennes, auto_unbox = TRUE)

# 2. Patcher l'app
html <- paste(readLines(app, warn = FALSE), collapse = "\n")

# 2a. Données embarquées à la place de DATA vide
html <- sub('var DATA={};',
            paste0('var DATA={"35238":', json_rennes, '};'),
            html, fixed = TRUE)

# 2b. Démarrage direct sur Rennes (au lieu de l'écran d'import)
html <- sub('loadCity(null);', 'loadCity("35238");', html, fixed = TRUE)

# 2c. Masquer le bouton d'import (inutile pour les lecteurs du mémoire)
html <- sub('<button class="btn" id="importbtn">Importer les données…</button>',
            '<button class="btn" id="importbtn" style="display:none">Importer les données…</button>',
            html, fixed = TRUE)

out <- file.path(DOSSIER, "reseau-rennes.html")
writeLines(html, out)
cat("OK :", out, "-", round(file.size(out) / 1024), "Ko\n")