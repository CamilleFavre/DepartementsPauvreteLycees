# =====================================================================
# CODE 2 (corrigé) - Indicateurs d'hétérogénéité intra-communale
# à partir des carreaux Filosofi 200 m
#
# Principales corrections par rapport à la version précédente :
#  1. Pas d'imputation des NA à 0 (revenu non renseigné = exclu, pas nul)
#  2. Voisinage ET standardisation calculés DANS la commune
#     (plus de voisins à travers les limites communales, plus de
#      z-score départemental)
#  3. Carreaux sans voisin exclus (NA), au lieu de lag = 0
#  4. Vrai LISA : p-values par permutation + correction FDR ;
#     "NS" = non significatif, pas "quadrant inconnu"
#  5. Parts (et non effectifs) de carreaux HH/LL/HL/LH
#  6. Gini pondéré corrigé (le premier segment de la courbe de Lorenz
#     était omis) ; NA (et non 0) quand l'indicateur n'est pas défini
#  7. Suppression des Gini/Theil non pondérés sur des comptages
#     (men_pauv, log_soc) : ils mesuraient la concentration, pas
#     l'hétérogénéité, et dépendaient mécaniquement de la taille.
#     Remplacés par un indice de dissimilarité (pauvres / non pauvres).
#  8. Option pour agréger les arrondissements de Paris/Lyon/Marseille
# =====================================================================

library(arrow)
library(dplyr)
library(jsonlite)
library(spdep)

# ---------------------------------------------------------------------
# 0. PARAMÈTRES
# ---------------------------------------------------------------------
DOSSIER_ECO <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/data_eco"
DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"

SEUIL_POP_SPATIAL <- 10000   # calcul spatial (LISA) seulement au-dessus de ce seuil (pop. de la commune)
SEUIL_POP_LARGE   <- 10000   # définition de df_poplarge
MIN_CARREAUX      <- 30      # nb minimal de carreaux (avec voisin) pour un indicateur spatial
DIST_VOISIN       <- 285     # m : couvre les 8 voisins (diagonale = 283 m) d'un carreau de 200 m
NSIM              <- 499     # nb de permutations pour le Moran local / global
ALPHA             <- 0.05    # seuil de significativité (après correction FDR)
AGREGER_PLM       <- FALSE   # TRUE = Paris, Lyon, Marseille traitées comme 1 seule commune
# (FALSE = chaque arrondissement est une unité)
# Attention : Code 3 doit utiliser le MÊME niveau pour les lycées.

set.seed(2025)

# ---------------------------------------------------------------------
# 1. FONCTIONS D'INDICATEURS
# ---------------------------------------------------------------------

# Gini pondéré (aire sous la courbe de Lorenz, règle des trapèzes,
# en partant bien de (0,0))
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

# Theil T pondéré (entropie généralisée, paramètre 1)
weighted_theil <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & x > 0 & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(NA_real_)
  w <- w / sum(w)
  rel <- x / sum(w * x)
  sum(w * rel * log(rel))
}

# Coefficient de variation pondéré
weighted_cv <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(NA_real_)
  m <- sum(w * x) / sum(w)
  if (m <= 0) return(NA_real_)
  sqrt(sum(w * (x - m)^2) / sum(w)) / m
}

# Quantiles pondérés (interpolation sur la fonction de répartition pondérée)
weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(rep(NA_real_, length(probs)))
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- (cumsum(w) - 0.5 * w) / sum(w)
  approx(cw, x, xout = probs, rule = 2, ties = "ordered")$y
}

# Indice de dissimilarité de Duncan entre deux groupes (a et b = effectifs par carreau)
dissimilarity <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 2 || sum(a) <= 0 || sum(b) <= 0) return(NA_real_)
  0.5 * sum(abs(a / sum(a) - b / sum(b)))
}

# Codes des arrondissements -> commune globale (Paris, Lyon, Marseille)
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
#    nb et lw sont déjà construits sur les carreaux de cette commune
# ---------------------------------------------------------------------
lisa_variable <- function(x, nb, lw, suffixe) {
  
  noms <- paste0(c("moran_I_", "moran_p_", "part_HH_", "part_LL_", "part_HL_",
                   "part_LH_", "part_NS_", "contraste_abs_", "contraste_rel_"), suffixe)
  vide <- setNames(as.list(rep(NA_real_, length(noms))), noms)
  
  if (!all(is.finite(x)) || sd(x) == 0) return(vide)
  
  # Standardisation DANS la commune
  z  <- as.numeric(scale(x))
  wz <- lag.listw(lw, z, zero.policy = TRUE)
  
  # Moran global (permutation)
  mi <- moran.mc(x, lw, nsim = NSIM, zero.policy = TRUE)
  
  # Moran local par permutation, p-values corrigées (FDR)
  loc <- localmoran_perm(x, lw, nsim = NSIM, zero.policy = TRUE, iseed = 1)
  col_p <- grep("Sim", colnames(loc))[1]
  if (is.na(col_p)) return(vide)
  p_adj <- p.adjust(loc[, col_p], method = "fdr")
  signif <- !is.na(p_adj) & p_adj < ALPHA
  
  quad <- rep("NS", length(x))
  quad[signif & z > 0 & wz > 0] <- "HH"
  quad[signif & z < 0 & wz < 0] <- "LL"
  quad[signif & z > 0 & wz < 0] <- "HL"
  quad[signif & z < 0 & wz > 0] <- "LH"
  
  # Contraste local : écart absolu à la moyenne des voisins
  voisins_moy <- lag.listw(lw, x, zero.policy = TRUE)
  contraste   <- abs(x - voisins_moy)
  
  out <- list(
    unname(mi$statistic), mi$p.value,
    mean(quad == "HH"), mean(quad == "LL"), mean(quad == "HL"),
    mean(quad == "LH"), mean(quad == "NS"),
    mean(contraste),
    ifelse(mean(x) > 0, mean(contraste) / mean(x), NA_real_)
  )
  setNames(out, noms)
}

# Indicateurs spatiaux pour UNE commune (data.frame de ses carreaux)
spatial_commune <- function(d) {
  
  vide <- tibble(
    part_carreaux_isoles = NA_real_, n_carreaux_spatial = NA_integer_,
    !!!setNames(as.list(rep(NA_real_, 18)),
                c(paste0(c("moran_I_", "moran_p_", "part_HH_", "part_LL_", "part_HL_",
                           "part_LH_", "part_NS_", "contraste_abs_", "contraste_rel_"), "rev"),
                  paste0(c("moran_I_", "moran_p_", "part_HH_", "part_LL_", "part_HL_",
                           "part_LH_", "part_NS_", "contraste_abs_", "contraste_rel_"), "pauv")))
  )
  
  # Carreaux exploitables : revenu et taux de pauvreté observés (pas d'imputation)
  d <- d[is.finite(d$revenu_moyen) & is.finite(d$taux_pauvrete) &
           !is.na(d$x_l93) & !is.na(d$y_l93) & d$ind > 0, ]
  n0 <- nrow(d)
  if (n0 < MIN_CARREAUX) return(vide)
  
  # 1er voisinage -> on retire les carreaux isolés (aucun voisin dans la commune)
  coords <- cbind(d$x_l93, d$y_l93)
  nb0 <- dnearneigh(coords, d1 = 0, d2 = DIST_VOISIN)
  connecte <- card(nb0) > 0
  if (sum(connecte) < MIN_CARREAUX) return(vide)
  d <- d[connecte, ]
  
  # Le voisinage est recalculé sur les carreaux restants
  # (retirer des carreaux isolés ne crée pas de nouveaux isolés)
  coords <- cbind(d$x_l93, d$y_l93)
  nb <- dnearneigh(coords, d1 = 0, d2 = DIST_VOISIN)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  res_rev  <- lisa_variable(d$revenu_moyen,   nb, lw, "rev")
  res_pauv <- lisa_variable(d$taux_pauvrete,  nb, lw, "pauv")
  
  tibble(
    part_carreaux_isoles = 1 - nrow(d) / n0,
    n_carreaux_spatial   = nrow(d),
    !!!res_rev, !!!res_pauv
  )
}

# ---------------------------------------------------------------------
# 3. RÉFÉRENTIEL DES COMMUNES (noms)
# ---------------------------------------------------------------------
message("Récupération de la table des communes...")
communes_ref <- fromJSON("https://geo.api.gouv.fr/communes?fields=nom,code&format=json") %>%
  rename(code_insee = code, nom_commune = nom)

# ---------------------------------------------------------------------
# 4. BOUCLE SUR LES DÉPARTEMENTS
# ---------------------------------------------------------------------
files <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*\\.parquet", full.names = TRUE)
n_files <- length(files)
results_list <- vector("list", n_files)
start_total <- Sys.time()

colonnes_requises <- c("lcog_geo", "x_l93", "y_l93", "ind", "men", "men_pauv", "revenu_moyen")

for (i in seq_along(files)) {
  t0 <- Sys.time()
  df <- read_parquet(files[i])
  manquantes <- setdiff(colonnes_requises, names(df))
  if (length(manquantes) > 0)
    stop("Colonnes manquantes dans ", basename(files[i]), " : ", paste(manquantes, collapse = ", "))
  
  df <- df %>%
    mutate(
      code_brut      = substr(as.character(lcog_geo), 1, 5),
      code_insee     = if (AGREGER_PLM) map_plm(code_brut) else code_brut,
      # taux de pauvreté : NA (et non 0) si aucun ménage
      taux_pauvrete  = ifelse(!is.na(men) & men > 0, men_pauv / men * 100, NA_real_),
      men_non_pauv   = pmax(men - men_pauv, 0)
    )
  
  # ---- 4a. Indicateurs non spatiaux (toutes communes) ----
  res_nonspatial <- df %>%
    group_by(code_insee) %>%
    summarise(
      nb_carreaux        = n(),
      pop_totale         = sum(ind, na.rm = TRUE),
      men_totale         = sum(men, na.rm = TRUE),
      part_carreaux_rev  = mean(is.finite(revenu_moyen) & ind > 0, na.rm = TRUE),  # qualité : couverture du revenu
      taux_pauv_commune  = ifelse(sum(men, na.rm = TRUE) > 0,
                                  sum(men_pauv, na.rm = TRUE) / sum(men, na.rm = TRUE) * 100, NA_real_),
      revenu_moyen_com   = weighted.mean(revenu_moyen, ind, na.rm = TRUE),
      
      # Hétérogénéité du revenu entre carreaux, pondérée par la population
      gini_revenu  = weighted_gini(revenu_moyen, ind),
      theil_revenu = weighted_theil(revenu_moyen, ind),
      cv_revenu    = weighted_cv(revenu_moyen, ind),
      ratio_p90_p10_rev = {
        q <- weighted_quantile(revenu_moyen, ind, c(0.1, 0.9))
        ifelse(is.finite(q[1]) && q[1] > 0, q[2] / q[1], NA_real_)
      },
      
      # Ségrégation résidentielle : pauvres vs non pauvres
      dissim_pauv = dissimilarity(men_pauv, men_non_pauv),
      
      # Concentration de la population (densité) : à NE PAS lire comme de l'hétérogénéité sociale
      gini_densite = ifelse(n() > 1, weighted_gini(ind, rep(1, n())), NA_real_),
      .groups = "drop"
    )
  
  # ---- 4b. Indicateurs spatiaux (communes assez peuplées seulement) ----
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
cat("Traitement terminé !\n")

# ---------------------------------------------------------------------
# 5. ASSEMBLAGE ET EXPORT
# ---------------------------------------------------------------------
df_final <- bind_rows(results_list) %>%
  left_join(communes_ref, by = "code_insee") %>%
  mutate(nom_commune = coalesce(nom_commune, code_insee)) %>%  # arrondissements : pas de nom dans l'API
  select(code_insee, nom_commune, everything())

# Heterogeneite synthétique de contraste local : "mixité" = parts de valeurs aberrantes spatiales
df_final <- df_final %>%
  mutate(
    part_mixte_rev   = part_HL_rev  + part_LH_rev,    # carreaux riches/pauvres entourés de l'inverse
    part_segreg_rev  = part_HH_rev  + part_LL_rev,    # carreaux appartenant à un cluster homogène
    part_mixte_pauv  = part_HL_pauv + part_LH_pauv,
    part_segreg_pauv = part_HH_pauv + part_LL_pauv
  )

write_parquet(df_final, file.path(DOSSIER_OUT, "df_final.parquet"))

df_poplarge <- df_final %>%
  filter(pop_totale > SEUIL_POP_LARGE) %>%
  arrange(desc(pop_totale))

write_parquet(df_poplarge, file.path(DOSSIER_OUT, "df_poplarge.parquet"))

# ---------------------------------------------------------------------
# 6. CONTRÔLES RAPIDES (à regarder avant toute régression)
# ---------------------------------------------------------------------
cat("\nNb de communes (pop > ", SEUIL_POP_LARGE, ") : ", nrow(df_poplarge), "\n", sep = "")
cat("Couverture du revenu (part de carreaux renseignés), résumé :\n")
print(summary(df_poplarge$part_carreaux_rev))
cat("Part de carreaux isolés (exclus du LISA), résumé :\n")
print(summary(df_poplarge$part_carreaux_isoles))
cat("Corrélations entre indicateurs principaux :\n")
print(round(cor(df_poplarge[, c("gini_revenu", "theil_revenu", "cv_revenu", "dissim_pauv",
                                "moran_I_rev", "contraste_rel_rev", "part_mixte_rev")],
                use = "pairwise.complete.obs"), 2))
cat("Corrélation avec la taille (log pop) :\n")
print(round(cor(log(df_poplarge$pop_totale),
                df_poplarge[, c("gini_revenu", "theil_revenu", "dissim_pauv", "moran_I_rev")],
                use = "pairwise.complete.obs"), 2))