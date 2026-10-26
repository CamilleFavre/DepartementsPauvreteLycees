# =====================================================================
# CODE 2 (v2) - Indicateurs d'heterogeneite intra-communale
# a partir des carreaux Filosofi 200 m
#
# Principales corrections par rapport a la version precedente :
#  1. Pas d'imputation des NA a 0 (revenu non renseigne = exclu, pas nul)
#  2. Voisinage ET standardisation calcules DANS la commune
#  3. Gini pondere corrige ; NA (et non 0) quand l'indicateur n'est pas defini
#  4. Suppression des Gini/Theil non ponderes sur des comptages ;
#     remplaces par un indice de dissimilarite (pauvres / non pauvres)
#  5. Option pour agreger les arrondissements de Paris/Lyon/Marseille
#
# NOUVEAU en v2 :
#  6. PLUS DE FILTRE DE SIGNIFICATIVITE sur les quadrants HH/LL/HL/LH.
#     Typologie descriptive exhaustive : chaque carreau est classe
#     selon (z, wz), sans permutation locale ni FDR. Justification :
#     avec NSIM permutations, la p-value est plancher a 1/(NSIM+1) ;
#     le FDR multiplie par n tests -> pour une commune de n < ~300
#     carreaux, AUCUN carreau ne peut jamais passer le seuil, quelle
#     que soit la realite des clusters (verifie empiriquement : 100 %
#     de NS sur les communes rurales). Le test du Moran GLOBAL
#     (moran.mc) est conserve : il n'est pas concerne par le probleme.
#  7. PROBLEME DE FRONTIERE TRAITE :
#     - les carreaux isoles (zero voisin a < 285 m dans la commune,
#       typiquement en extremite de commune ou dans un ilot entoure de
#       carreaux sans menage) ne sont plus EXCLUS : ils recourent a un
#       voisinage kNN (k plus proches voisins dans la commune).
#     - nouvelle colonne part_carreaux_bord : part des carreaux ayant
#       un voisinage incomplet (< 8 voisins) -- effet de bord residual,
#       a utiliser comme variable de controle ou flag qualite.
#     - NB : les carreaux de bord gardent moins de voisins (3-7), c'est
#       incompressible si on reste intra-communal ; c'est maintenant
#       mesure et documente au lieu d'ignore.
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

SEUIL_POP_SPATIAL <- 10000
SEUIL_POP_LARGE   <- 10000
MIN_CARREAUX      <- 30
DIST_VOISIN       <- 285     # m : couvre les 8 voisins d'un carreau de 200 m
K_KNN             <- 4       # voisins pour le recours kNN des carreaux isoles
NSIM              <- 499     # permutations pour le Moran GLOBAL seulement
AGREGER_PLM       <- FALSE   # TRUE = Paris, Lyon, Marseille en 1 commune
# (Code 3 doit utiliser le MEME niveau pour les lycees)

set.seed(2025)

# ---------------------------------------------------------------------
# 1. FONCTIONS D'INDICATEURS
# ---------------------------------------------------------------------

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

weighted_theil <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & x > 0 & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(NA_real_)
  w <- w / sum(w)
  rel <- x / sum(w * x)
  sum(w * rel * log(rel))
}

weighted_cv <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(NA_real_)
  m <- sum(w * x) / sum(w)
  if (m <= 0) return(NA_real_)
  sqrt(sum(w * (x - m)^2) / sum(w)) / m
}

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}

dissimilarity <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 2 || sum(a) <= 0 || sum(b) <= 0) return(NA_real_)
  0.5 * sum(abs(a / sum(a) - b / sum(b)))
}

map_plm <- function(code) {
  dplyr::case_when(
    grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
    grepl("^6938[1-9]$", code)             ~ "69123",
    grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
    TRUE ~ code
  )
}

# ---------------------------------------------------------------------
# 2. INDICATEURS SPATIAUX (LISA) POUR UNE VARIABLE, DANS UNE COMMUNE
#    Quadrants SANS filtre de significativite (typologie exhaustive)
# ---------------------------------------------------------------------
lisa_variable <- function(x, nb, lw, suffixe) {
  
  noms <- paste0(c("moran_I_", "moran_p_",
                   "part_HH_", "part_LL_", "part_HL_", "part_LH_",
                   "contraste_abs_", "contraste_rel_"), suffixe)
  vide <- setNames(as.list(rep(NA_real_, length(noms))), noms)
  
  if (!all(is.finite(x)) || sd(x) == 0) return(vide)
  
  # Standardisation DANS la commune
  z  <- as.numeric(scale(x))
  wz <- lag.listw(lw, z, zero.policy = TRUE)
  
  # Moran global par permutation (indicateur continu, non filtre)
  mi <- moran.mc(x, lw, nsim = NSIM, zero.policy = TRUE)
  
  # Quadrants : chaque carreau est classe, sans test statistique.
  # Exhaustif : les cas d'egalite exacte (z=0 ou wz=0, rares) sont
  # rattaches au quadrant positif le plus proche.
  quad <- ifelse(wz >= 0,
                 ifelse(z >= 0, "HH", "LH"),
                 ifelse(z >= 0, "HL", "LL"))
  
  # Contraste local : ecart absolu a la moyenne des voisins
  voisins_moy <- lag.listw(lw, x, zero.policy = TRUE)
  contraste   <- abs(x - voisins_moy)
  
  out <- list(
    unname(mi$statistic), mi$p.value,
    mean(quad == "HH"), mean(quad == "LL"),
    mean(quad == "HL"), mean(quad == "LH"),
    mean(contraste),
    ifelse(mean(x) > 0, mean(contraste) / mean(x), NA_real_)
  )
  setNames(out, noms)
}

# ---------------------------------------------------------------------
# 3. VOISINAGE INTRA-COMMUNAL AVEC RECOURS kNN POUR LES ISOLES
#    - dnearneigh (285 m) : voisinage principal
#    - les carreaux sans AUCUN voisin a 285 m recourent aux K_KNN plus
#      proches carreaux de la commune (plus d'exclusion)
#    - part_carreaux_bord : part des carreaux a voisinage incomplet
# ---------------------------------------------------------------------
voisinage_commune <- function(coords) {
  nb <- dnearneigh(coords, d1 = 0, d2 = DIST_VOISIN)
  card0 <- card(nb)
  
  part_isoles <- mean(card0 == 0)   # recourent au kNN
  part_bord   <- mean(card0 < 8)    # voisinage incomplet (bord, trous)
  
  if (any(card0 == 0)) {
    kk   <- min(K_KNN, nrow(coords) - 1)
    knnb <- as.list(knn2nb(knearneigh(coords, k = kk)))
    nbl  <- as.list(nb)
    for (j in which(card0 == 0)) nbl[[j]] <- knnb[[j]]
    nb <- structure(nbl, class = "nb",
                    region.id = attr(nb, "region.id"),
                    sym = attr(nb, "sym"))
  }
  list(nb = nb, part_isoles = part_isoles, part_bord = part_bord)
}

# Indicateurs spatiaux pour UNE commune
spatial_commune <- function(d) {
  
  vide <- tibble(
    part_carreaux_isoles = NA_real_,
    part_carreaux_bord   = NA_real_,
    n_carreaux_spatial   = NA_integer_,
    !!!setNames(as.list(rep(NA_real_, 16)),
                c(paste0(c("moran_I_", "moran_p_",
                           "part_HH_", "part_LL_", "part_HL_", "part_LH_",
                           "contraste_abs_", "contraste_rel_"), "rev"),
                  paste0(c("moran_I_", "moran_p_",
                           "part_HH_", "part_LL_", "part_HL_", "part_LH_",
                           "contraste_abs_", "contraste_rel_"), "pauv")))
  )
  
  # Carreaux exploitables : revenu et taux de pauvrete observes (pas d'imputation)
  d <- d[is.finite(d$revenu_moyen) & is.finite(d$taux_pauvrete) &
           !is.na(d$x_l93) & !is.na(d$y_l93) & d$ind > 0, ]
  n0 <- nrow(d)
  if (n0 < MIN_CARREAUX) return(vide)
  
  coords <- cbind(d$x_l93, d$y_l93)
  v      <- voisinage_commune(coords)
  lw     <- nb2listw(v$nb, style = "W", zero.policy = TRUE)
  
  res_rev  <- lisa_variable(d$revenu_moyen,  v$nb, lw, "rev")
  res_pauv <- lisa_variable(d$taux_pauvrete, v$nb, lw, "pauv")
  
  tibble(
    part_carreaux_isoles = v$part_isoles,   # info : part ayant recours au kNN
    part_carreaux_bord   = v$part_bord,     # info : part a voisinage incomplet
    n_carreaux_spatial   = nrow(d),          # PLUS AUCUN carreau n'est exclu
    !!!res_rev, !!!res_pauv
  )
}

# ---------------------------------------------------------------------
# 4. REFERENTIEL DES COMMUNES (noms)
# ---------------------------------------------------------------------
message("Recuperation de la table des communes...")
communes_ref <- fromJSON("https://geo.api.gouv.fr/communes?fields=nom,code&format=json") %>%
  rename(code_insee = code, nom_commune = nom)

# ---------------------------------------------------------------------
# 5. BOUCLE SUR LES DEPARTEMENTS
# ---------------------------------------------------------------------
files <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*[.]parquet$", full.names = TRUE)
n_files <- length(files)
results_list <- vector("list", n_files)
start_total <- Sys.time()

colonnes_requises <- c("lcog_geo", "x_l93", "y_l93", "ind", "men", "men_pauv", "revenu_moyen")

for (i in seq_along(files)) {
  t0 <- Sys.time()
  df <- read_parquet(files[i])
  manquantes <- setdiff(colonnes_requises, names(df))
  if (length(manquantes) > 0)
    stop("Colonnes manquantes dans ", basename(files[i]), " : ",
         paste(manquantes, collapse = ", "))
  
  df <- df %>%
    mutate(
      code_brut      = substr(as.character(lcog_geo), 1, 5),
      code_insee     = if (AGREGER_PLM) map_plm(code_brut) else code_brut,
      taux_pauvrete  = ifelse(!is.na(men) & men > 0, men_pauv / men * 100, NA_real_),
      men_non_pauv   = pmax(men - men_pauv, 0)
    )
  
  # ---- 5a. Indicateurs non spatiaux (toutes communes) ----
  res_nonspatial <- df %>%
    group_by(code_insee) %>%
    summarise(
      nb_carreaux        = n(),
      pop_totale         = sum(ind, na.rm = TRUE),
      men_totale         = sum(men, na.rm = TRUE),
      part_carreaux_rev  = mean(is.finite(revenu_moyen) & ind > 0, na.rm = TRUE),
      taux_pauv_commune  = ifelse(sum(men, na.rm = TRUE) > 0,
                                  sum(men_pauv, na.rm = TRUE) / sum(men, na.rm = TRUE) * 100,
                                  NA_real_),
      revenu_moyen_com   = weighted.mean(revenu_moyen, ind, na.rm = TRUE),
      gini_revenu        = weighted_gini(revenu_moyen, ind),
      theil_revenu       = weighted_theil(revenu_moyen, ind),
      cv_revenu          = weighted_cv(revenu_moyen, ind),
      ratio_p90_p10_rev  = {
        q <- weighted_quantile(revenu_moyen, ind, c(0.1, 0.9))
        ifelse(is.finite(q[1]) && q[1] > 0, q[2] / q[1], NA_real_)
      },
      dissim_pauv        = dissimilarity(men_pauv, men_non_pauv),
      gini_densite       = ifelse(n() > 1, weighted_gini(ind, rep(1, n())), NA_real_),
      .groups = "drop"
    )
  
  # ---- 5b. Indicateurs spatiaux (communes assez peuplees seulement) ----
  communes_cibles <- res_nonspatial$code_insee[res_nonspatial$pop_totale >= SEUIL_POP_SPATIAL]
  df_cible <- df %>% filter(code_insee %in% communes_cibles)
  
  if (nrow(df_cible) > 0) {
    liste_com <- split(df_cible, df_cible$code_insee)
    res_spatial <- bind_rows(lapply(liste_com, spatial_commune), .id = "code_insee")
    res_nonspatial <- left_join(res_nonspatial, res_spatial, by = "code_insee")
  }
  
  results_list[[i]] <- res_nonspatial
  
  el <- as.numeric(Sys.time() - start_total, units = "secs")
  cat(sprintf("[%02d/%02d] %s | %d communes spatiales | %.1fs | reste ~%.0fs\n",
              i, n_files, basename(files[i]), length(communes_cibles),
              as.numeric(Sys.time() - t0, units = "secs"), (n_files - i) * el / i))
}
cat("Traitement termine !\n")

# ---------------------------------------------------------------------
# 6. ASSEMBLAGE ET EXPORT
# ---------------------------------------------------------------------
df_final <- bind_rows(results_list) %>%
  left_join(communes_ref, by = "code_insee") %>%
  mutate(nom_commune = coalesce(nom_commune, code_insee)) %>%
  select(code_insee, nom_commune, everything())

df_final <- df_final %>%
  mutate(
    part_mixte_rev    = part_HL_rev  + part_LH_rev,
    part_segreg_rev   = part_HH_rev  + part_LL_rev,
    part_mixte_pauv   = part_HL_pauv + part_LH_pauv,
    part_segreg_pauv  = part_HH_pauv + part_LL_pauv
  )

write_parquet(df_final, file.path(DOSSIER_OUT, "df_final.parquet"))

df_poplarge <- df_final %>%
  filter(pop_totale > SEUIL_POP_LARGE) %>%
  arrange(desc(pop_totale))

write_parquet(df_poplarge, file.path(DOSSIER_OUT, "df_poplarge.parquet"))

# ---------------------------------------------------------------------
# 7. CONTROLES RAPIDES (a regarder avant toute regression)
# ---------------------------------------------------------------------
cat("\nNb de communes (pop > ", SEUIL_POP_LARGE, ") : ", nrow(df_poplarge), "\n", sep = "")
cat("Couverture du revenu (part de carreaux renseignes), resume :\n")
print(summary(df_poplarge$part_carreaux_rev))
cat("Part de carreaux isoles (recours kNN), resume :\n")
print(summary(df_poplarge$part_carreaux_isoles))
cat("Part de carreaux de bord (voisinage incomplet), resume :\n")
print(summary(df_poplarge$part_carreaux_bord))
cat("Verification : parts HH+LL+HL+LH = 1 partout (quadrants exhaustifs) ?\n")
print(all(abs((df_poplarge$part_segreg_rev + df_poplarge$part_mixte_rev) - 1) < 1e-8,
          na.rm = TRUE))
cat("Correlations entre indicateurs principaux :\n")
print(round(cor(df_poplarge[, c("gini_revenu", "theil_revenu", "cv_revenu", "dissim_pauv",
                                "moran_I_rev", "contraste_rel_rev",
                                "part_segreg_rev", "part_mixte_rev")],
                use = "pairwise.complete.obs"), 2))
cat("Correlation avec la taille (log pop) :\n")
print(round(cor(log(df_poplarge$pop_totale),
                df_poplarge[, c("gini_revenu", "theil_revenu", "dissim_pauv",
                                "moran_I_rev", "part_segreg_rev", "part_carreaux_bord")],
                use = "pairwise.complete.obs"), 2))
cat("Biais de peripherie a surveiller : cor(part_carreaux_bord, taux_pauv_commune)\n")
print(cor(df_poplarge$part_carreaux_bord, df_poplarge$taux_pauv_commune,
          use = "pairwise.complete.obs"))