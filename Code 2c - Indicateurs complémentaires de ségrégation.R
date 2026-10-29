# =====================================================================
# CODE 2c - INDICATEURS COMPLÉMENTAIRES DE SÉGRÉGATION (VERSION CORRIGÉE)
#
# Entrée  : fichiers eco_dept_*.parquet (carreaux Filosofi 200m)
# Sortie  : df_segregation_complementaire.parquet (une ligne par commune)
#
# Indicateurs calculés :
#   - Exposition spatiale aux riches (Pxx) - version Massey & Denton
#   - Entropie relative pondérée (diversité)
#   - Diversité binaire pauvre/non-pauvre (indice de Simpson)
#   - Indice de polarisation (Wolfson)
#   - Part population dans top/bottom déciles/quintiles
#   - Distance moyenne pondérée au carreau le plus riche/pauvre
#   - Getis-Ord Gi* (concentration spatiale de la richesse)
#   - Indice de ségrégation multidimensionnel (revenu + pauvreté)
#
# À placer : APRÈS le Code 2, AVANT le Code 3b
# =====================================================================

library(arrow)
library(dplyr)
library(spdep)

# ---------------------------------------------------------------------
# 0. PARAMÈTRES (identiques au Code 2)
# ---------------------------------------------------------------------
os_name <- Sys.info()[["sysname"]]

if (os_name == "Windows") {
  DOSSIER_ECO <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/data_eco"
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
} else if (os_name == "Darwin") {
  DOSSIER_ECO <- "/Users/camillefavre/DepartementsPauvreteLycees/data_eco"
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
} else {
  stop("Système d'exploitation non pris en charge.")
}

SEUIL_POP_SPATIAL <- 10000   # communes avec assez de carreaux
DIST_VOISIN       <- 285     # m : 8 voisins d'un carreau 200m
MIN_CARREAUX      <- 30      # minimum de carreaux pour les indicateurs spatiaux

set.seed(2025)

# ---------------------------------------------------------------------
# 1. FONCTIONS AUXILIAIRES
# ---------------------------------------------------------------------

# Médiane pondérée (fonction autonome)
weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) return(NA_real_)
  
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

# Gini pondéré (repris du Code 2 pour autonomie)
weighted_gini <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2 || sum(x * w) <= 0) return(NA_real_)
  o <- order(x); x <- x[o]; w <- w[o]
  p <- c(0, cumsum(w) / sum(w))
  L <- c(0, cumsum(x * w) / sum(x * w))
  g <- 1 - sum(diff(p) * (L[-1] + L[-length(L)]))
  max(0, min(1, g))
}

# Indice de dissimilarité (repris du Code 2)
dissimilarity <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 2 || sum(a) <= 0 || sum(b) <= 0) return(NA_real_)
  0.5 * sum(abs(a / sum(a) - b / sum(b)))
}

# ---------------------------------------------------------------------
# 2. INDICATEURS COMPLÉMENTAIRES
# ---------------------------------------------------------------------

# Exposition spatiale aux riches (Pxx) - Massey & Denton (1988)
# Probabilité qu'un riche rencontre un riche dans son voisinage immédiat
exposition_riches_spatial <- function(x, y, revenu, ind, dist_voisin = 285) {
  ok <- is.finite(x) & is.finite(y) & is.finite(revenu) & is.finite(ind) & ind > 0
  if (sum(ok) < 5) return(NA_real_)
  
  x <- x[ok]; y <- y[ok]; rev <- revenu[ok]; pop <- ind[ok]
  
  # Définir riche (top 25% du revenu communal)
  seuil <- quantile(rev, 0.75, na.rm = TRUE)
  pop_riche <- ifelse(rev >= seuil, pop, 0)
  total_riche <- sum(pop_riche)
  
  if (total_riche == 0) return(NA_real_)
  
  # Voisinage spatial (distance)
  coords <- cbind(x, y)
  nb <- dnearneigh(coords, d1 = 0, d2 = dist_voisin)
  
  # Proportion de riches dans le voisinage fermé de chaque carreau
  prop_riche_voisin <- sapply(seq_along(nb), function(i) {
    voisins <- c(i, nb[[i]])  # voisinage fermé (inclut le carreau lui-même)
    sum(pop_riche[voisins]) / sum(pop[voisins])
  })
  
  # Exposition pondérée : moyenne de la proportion de riches dans le voisinage,
  # pondérée par la population riche de chaque carreau
  exposition <- sum((pop_riche / total_riche) * prop_riche_voisin, na.rm = TRUE)
  
  return(exposition)
}

# Entropie relative pondérée (indice de diversité)
# 0 = un seul groupe présent, 1 = équirépartition parfaite entre groupes
entropie_relative_ponderee <- function(groupes, poids) {
  ok <- is.finite(groupes) & is.finite(poids) & poids > 0
  if (sum(ok) < 2) return(NA_real_)
  
  groupes <- groupes[ok]
  poids <- poids[ok]
  
  # Nombre d'individus par groupe
  tab <- tapply(poids, groupes, sum)
  tab <- tab[!is.na(tab) & tab > 0]
  
  if (length(tab) < 2) return(NA_real_)
  
  # Parts de population
  p <- tab / sum(tab)
  
  # Entropie de Shannon
  H <- -sum(p * log(p))
  H_max <- log(length(tab))
  
  if (H_max == 0) return(NA_real_)
  
  return(H / H_max)
}

# Diversité binaire pauvre/non-pauvre (indice de Simpson)
# 0 = un seul groupe, 0.5 = équirépartition parfaite (50/50)
diversite_binaire <- function(men_pauv, men) {
  ok <- is.finite(men_pauv) & is.finite(men) & men > 0
  if (sum(ok) < 2) return(NA_real_)
  
  pauvres <- sum(men_pauv[ok], na.rm = TRUE)
  non_pauvres <- sum(pmax(men[ok] - men_pauv[ok], 0), na.rm = TRUE)
  total <- pauvres + non_pauvres
  
  if (total == 0) return(NA_real_)
  
  p_pauv <- pauvres / total
  p_non_pauv <- non_pauvres / total
  
  # Indice de Simpson (diversité)
  1 - (p_pauv^2 + p_non_pauv^2)
}

# Indice de polarisation (Wolfson 1994) - concentration aux extrêmes
indice_polarisation <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 3) return(NA_real_)
  
  med <- weighted_median(x, w)
  moy <- sum(x * w) / sum(w)
  var <- sum(w * (x - moy)^2) / sum(w)
  sd <- sqrt(var)
  
  if (sd == 0) return(NA_real_)
  abs(moy - med) / sd
}

# Part de population dans les déciles/quintiles extrêmes
part_pop_extremes <- function(revenu, ind, type = "top", p = 0.10) {
  ok <- is.finite(revenu) & is.finite(ind) & ind > 0
  if (sum(ok) == 0) return(NA_real_)
  
  rev <- revenu[ok]; pop <- ind[ok]
  
  seuil <- if (type == "top") {
    quantile(rev, 1 - p, na.rm = TRUE)
  } else {
    quantile(rev, p, na.rm = TRUE)
  }
  
  masque <- if (type == "top") rev >= seuil else rev <= seuil
  
  sum(pop[masque]) / sum(pop)
}

# Distance moyenne pondérée au carreau le plus riche (ou pauvre si rev négatif)
distance_moyenne_extreme <- function(x, y, revenu, ind) {
  ok <- is.finite(x) & is.finite(y) & is.finite(revenu) & is.finite(ind) & ind > 0
  if (sum(ok) < 2) return(NA_real_)
  
  x <- x[ok]; y <- y[ok]; rev <- revenu[ok]; pop <- ind[ok]
  
  # Carreau extrême (max si revenu positif, min si négatif pour "pauvre")
  i_ext <- which.max(rev)
  x_ext <- x[i_ext]; y_ext <- y[i_ext]
  
  # Distance euclidienne en km
  dist_km <- sqrt((x - x_ext)^2 + (y - y_ext)^2) / 1000
  
  # Distance moyenne pondérée par la population
  weighted.mean(dist_km, pop)
}

# Getis-Ord Gi* pour les hot-spots de richesse
getis_ord_riches <- function(x, y, revenu, ind, dist_max = 1000) {
  ok <- is.finite(x) & is.finite(y) & is.finite(revenu) & is.finite(ind) & ind > 0
  if (sum(ok) < 5) return(NA_real_)
  
  x <- x[ok]; y <- y[ok]; rev <- revenu[ok]
  
  coords <- cbind(x, y)
  
  # Voisinage par distance
  nb <- dnearneigh(coords, d1 = 0, d2 = dist_max)
  
  # Si pas de voisins (cas isolé), retourner NA
  if (all(card(nb) == 0)) return(NA_real_)
  
  # Calculer Gi* local
  gi <- tryCatch({
    localG(rev, nb2listw(nb, style = "W", zero.policy = TRUE))
  }, error = function(e) rep(NA_real_, length(rev)))
  
  # Agréger : moyenne des Gi* positifs significatifs (> 1.96)
  gi_sig <- gi[!is.na(gi) & gi > 1.96]
  if (length(gi_sig) == 0) return(NA_real_)
  
  return(mean(gi_sig))
}

# Indice de ségrégation multidimensionnel (revenu + pauvreté)
segregation_multidim <- function(revenu, ind, men_pauv, men) {
  gini_rev <- weighted_gini(revenu, ind)
  dissim <- dissimilarity(men_pauv, pmax(men - men_pauv, 0))
  
  if (is.na(gini_rev) || is.na(dissim)) return(NA_real_)
  
  # Moyenne des deux indicateurs normalisés [0,1]
  (gini_rev + dissim) / 2
}

# ---------------------------------------------------------------------
# 3. FONCTION PRINCIPALE PAR COMMUNE
# ---------------------------------------------------------------------
calculer_indicateurs_commune <- function(d) {
  
  # d = data.frame d'une commune avec colonnes :
  # lcog_geo, x_l93, y_l93, ind, men, men_pauv, revenu_moyen
  
  # Filtrer carreaux valides
  d <- d %>%
    filter(is.finite(revenu_moyen), is.finite(ind), ind > 0,
           !is.na(x_l93), !is.na(y_l93))
  
  n_car <- nrow(d)
  
  if (n_car < MIN_CARREAUX) {
    return(tibble(
      exposition_riches = NA_real_,
      entropie_diversite = NA_real_,
      diversite_binaire = NA_real_,
      polarisation = NA_real_,
      part_pop_top10 = NA_real_,
      part_pop_bottom10 = NA_real_,
      part_pop_top20 = NA_real_,
      part_pop_bottom20 = NA_real_,
      dist_moy_carreau_riche = NA_real_,
      dist_moy_carreau_pauvre = NA_real_,
      getis_ord_riches = NA_real_,
      segregation_multidim = NA_real_,
      n_carreaux_utilises = n_car
    ))
  }
  
  # Créer déciles de revenu pour l'entropie
  d <- d %>%
    mutate(decile_rev = ntile(revenu_moyen, 10))
  
  # Calculs
  res <- tibble(
    n_carreaux_utilises = n_car,
    
    # Exposition spatiale aux riches (Pxx)
    exposition_riches = exposition_riches_spatial(d$x_l93, d$y_l93, d$revenu_moyen, d$ind, DIST_VOISIN),
    
    # Entropie pondérée (mixité des déciles selon la population)
    entropie_diversite = entropie_relative_ponderee(d$decile_rev, d$ind),
    
    # Diversité binaire pauvre/non-pauvre
    diversite_binaire = diversite_binaire(d$men_pauv, d$men),
    
    # Polarisation
    polarisation = indice_polarisation(d$revenu_moyen, d$ind),
    
    # Parts population dans les extrêmes
    part_pop_top10 = part_pop_extremes(d$revenu_moyen, d$ind, "top", 0.10),
    part_pop_bottom10 = part_pop_extremes(d$revenu_moyen, d$ind, "bottom", 0.10),
    part_pop_top20 = part_pop_extremes(d$revenu_moyen, d$ind, "top", 0.20),
    part_pop_bottom20 = part_pop_extremes(d$revenu_moyen, d$ind, "bottom", 0.20),
    
    # Distances moyennes aux carreaux extrêmes
    dist_moy_carreau_riche = distance_moyenne_extreme(d$x_l93, d$y_l93, d$revenu_moyen, d$ind),
    dist_moy_carreau_pauvre = distance_moyenne_extreme(d$x_l93, d$y_l93, -d$revenu_moyen, d$ind),
    
    # Getis-Ord (hot-spots de richesse)
    getis_ord_riches = getis_ord_riches(d$x_l93, d$y_l93, d$revenu_moyen, d$ind),
    
    # Ségrégation multidimensionnelle
    segregation_multidim = segregation_multidim(d$revenu_moyen, d$ind, d$men_pauv, d$men)
  )
  
  return(res)
}

# ---------------------------------------------------------------------
# 4. BOUCLE SUR LES DÉPARTEMENTS
# ---------------------------------------------------------------------
files <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*[.]parquet$", full.names = TRUE)
n_files <- length(files)
results_list <- vector("list", n_files)

cat("Calcul des indicateurs complémentaires de ségrégation...\n")

for (i in seq_along(files)) {
  t0 <- Sys.time()
  
  df <- read_parquet(files[i])
  
  # Préparer les données
  df <- df %>%
    mutate(
      code_insee = substr(as.character(lcog_geo), 1, 5),
      taux_pauvrete = ifelse(!is.na(men) & men > 0, men_pauv / men * 100, NA_real_)
    )
  
  # Agréger par commune pour avoir la population totale
  pop_com <- df %>%
    group_by(code_insee) %>%
    summarise(pop_totale = sum(ind, na.rm = TRUE), .groups = "drop")
  
  # Filtrer les communes assez peuplées
  communes_cibles <- pop_com %>%
    filter(pop_totale >= SEUIL_POP_SPATIAL) %>%
    pull(code_insee)
  
  df_cible <- df %>% filter(code_insee %in% communes_cibles)
  
  if (nrow(df_cible) > 0) {
    # Calculer les indicateurs par commune
    liste_com <- split(df_cible, df_cible$code_insee)
    
    res_com <- bind_rows(lapply(liste_com, calculer_indicateurs_commune),
                         .id = "code_insee")
    
    results_list[[i]] <- res_com
  } else {
    results_list[[i]] <- tibble(code_insee = character())
  }
  
  el <- as.numeric(Sys.time() - t0, units = "secs")
  if (i %% 10 == 0 || i == n_files) {
    cat(sprintf("[%02d/%02d] %s | %d communes | %.1fs | reste ~%.0fs\n",
                i, n_files, basename(files[i]), length(communes_cibles),
                el, (n_files - i) * el / max(i, 1) * 10))
  }
}

# ---------------------------------------------------------------------
# 5. ASSEMBLAGE, STANDARDISATION ET SAUVEGARDE
# ---------------------------------------------------------------------
df_segregation <- bind_rows(results_list)

# Liste des indicateurs numériques à standardiser
indicateurs_num <- c("exposition_riches", "entropie_diversite", "diversite_binaire",
                     "polarisation", "part_pop_top10", "part_pop_bottom10",
                     "part_pop_top20", "part_pop_bottom20",
                     "dist_moy_carreau_riche", "dist_moy_carreau_pauvre",
                     "getis_ord_riches", "segregation_multidim")

# Standardiser (z-scores) pour comparabilité dans les régressions
df_segregation <- df_segregation %>%
  mutate(across(all_of(indicateurs_num),
                ~ as.numeric(scale(.)),  # z-score
                .names = "z_{.col}"))

# Sauvegarder
write_parquet(df_segregation, file.path(DOSSIER_OUT, "df_segregation_complementaire.parquet"))

# ---------------------------------------------------------------------
# 6. CONTRÔLES
# ---------------------------------------------------------------------
cat("\n=== CONTRÔLES ===\n")
cat("Communes traitées :", nrow(df_segregation), "\n\n")

cat("Résumé des indicateurs (valeurs brutes) :\n")
print(summary(df_segregation[, indicateurs_num]))

cat("\nCorrélations entre indicateurs (valeurs brutes) :\n")
cor_mat <- cor(df_segregation[, indicateurs_num], use = "pairwise.complete.obs")
print(round(cor_mat, 2))

cat("\nNombre de NA par indicateur :\n")
print(colSums(is.na(df_segregation[, indicateurs_num])))

cat("\nFichier sauvegardé : df_segregation_complementaire.parquet\n")
cat("Colonnes z-scores disponibles :", paste0("z_", indicateurs_num), "\n")