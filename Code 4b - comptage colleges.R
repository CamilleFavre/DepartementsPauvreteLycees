# =====================================================================
# CODE 4b (v3) - COLLÈGES : Poisson / NB2 / Hurdle + FE RÉGION ACADÉMIQUE
# + indicateurs de ségrégation complémentaires (Code 2c)
# Entrée  : df_analyse_colleges.parquet + df_segregation_complementaire.parquet
# Sortie  : modeles_comptage.rds
# =====================================================================
library(arrow); library(dplyr); library(sandwich); library(lmtest)
library(pscl); library(AER); library(sf); library(stringi); library(stringr)

os_name <- Sys.info()[["sysname"]]
if (os_name == "Windows") {
  DOSSIER     <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
  CHEMIN_GEO  <- file.path(DOSSIER_OUT, "fr-en-contour-academies-2020.geojson")
  DOSSIER_ECO <- file.path(DOSSIER_OUT, "data_eco")
} else if (os_name == "Darwin") {
  DOSSIER     <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
  CHEMIN_GEO  <- file.path(DOSSIER_OUT, "fr-en-contour-academies-2020.geojson")
  DOSSIER_ECO <- file.path(DOSSIER_OUT, "data_eco")
} else stop("Système d'exploitation non pris en charge.")

VAR_Y <- "nb_colleges_prives"
set.seed(2025)

# ---------------------------------------------------------------------
# 1. LECTURE + JOINTURE INDICATEURS SÉGRÉGATION (Code 2c)
# ---------------------------------------------------------------------
d <- read_parquet(file.path(DOSSIER_OUT, "df_analyse_colleges.parquet"))
df_segreg <- read_parquet(file.path(DOSSIER_OUT, "df_segregation_complementaire.parquet"))
d <- d %>% left_join(df_segreg, by = "code_insee")
cat("Communes avec indicateurs 2c :", sum(!is.na(d$z_exposition_riches)), "/", nrow(d), "\n")

# ---------------------------------------------------------------------
# 2. ACADÉMIE (fine) + RÉGION ACADÉMIQUE (13, pour FE/cluster)
# ---------------------------------------------------------------------
norm_acad <- function(x) {
  x <- tolower(stri_trans_general(as.character(x), "Latin-ASCII"))
  x <- trimws(gsub("^academie (de |d')?", "", x))
  gsub("[^a-z0-9]", "", x)
}

DEMI_CARREAU <- 100
files_eco <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*[.]parquet$", full.names = TRUE)
centroides <- bind_rows(lapply(files_eco, function(f) {
  read_parquet(f, col_select = c("lcog_geo","x_l93","y_l93","ind")) %>%
    mutate(code_insee = substr(as.character(lcog_geo), 1, 5)) %>%
    filter(code_insee %in% d$code_insee, is.finite(ind), ind > 0) %>%
    group_by(code_insee) %>%
    summarise(x = weighted.mean(x_l93 + DEMI_CARREAU, ind),
              y = weighted.mean(y_l93 + DEMI_CARREAU, ind), .groups = "drop")
})) %>% st_as_sf(coords = c("x","y"), crs = 2154)

carte_acad <- st_read(CHEMIN_GEO, quiet = TRUE) %>% st_transform(2154)
cand  <- names(carte_acad)[grepl("acad|lib|nom", tolower(names(carte_acad)))]
col_acad <- cand[1]
if (is.na(col_acad)) stop("Colonne académie introuvable dans le geojson")
carte_acad <- carte_acad %>%
  mutate(academie = norm_acad(.data[[col_acad]])) %>%
  filter(!academie %in% c("guyane","guadeloupe","mayotte","lareunion","martinique")) %>%
  select(academie, geometry)

acad_com <- st_join(st_transform(centroides, st_crs(carte_acad)), carte_acad) %>%
  st_drop_geometry() %>% select(code_insee, academie) %>% distinct(code_insee, .keep_all = TRUE)

dept2academie <- c(
  "02"="amiens","60"="amiens","80"="amiens","59"="lille","62"="lille",
  "75"="paris","77"="creteil","93"="creteil","94"="creteil",
  "78"="versailles","91"="versailles","92"="versailles","95"="versailles",
  "14"="normandie","27"="normandie","50"="normandie","61"="normandie","76"="normandie",
  "21"="dijon","58"="dijon","71"="dijon","89"="dijon",
  "25"="besancon","39"="besancon","70"="besancon","90"="besancon",
  "54"="nancymetz","55"="nancymetz","57"="nancymetz","88"="nancymetz",
  "08"="reims","10"="reims","51"="reims","52"="reims",
  "67"="strasbourg","68"="strasbourg",
  "44"="nantes","49"="nantes","53"="nantes","72"="nantes","85"="nantes",
  "22"="rennes","29"="rennes","35"="rennes","56"="rennes",
  "16"="poitiers","17"="poitiers","79"="poitiers","86"="poitiers",
  "18"="orleanstours","28"="orleanstours","36"="orleanstours","37"="orleanstours","41"="orleanstours",
  "24"="bordeaux","33"="bordeaux","40"="bordeaux","47"="bordeaux","64"="bordeaux",
  "19"="limoges","23"="limoges","87"="limoges",
  "03"="clermontferrand","15"="clermontferrand","43"="clermontferrand","63"="clermontferrand",
  "01"="lyon","42"="lyon","69"="lyon",
  "07"="grenoble","26"="grenoble","38"="grenoble","73"="grenoble","74"="grenoble",
  "04"="aixmarseille","05"="aixmarseille","13"="aixmarseille","84"="aixmarseille",
  "06"="nice","83"="nice",
  "11"="montpellier","30"="montpellier","34"="montpellier","48"="montpellier","66"="montpellier",
  "09"="toulouse","12"="toulouse","31"="toulouse","32"="toulouse","46"="toulouse",
  "65"="toulouse","81"="toulouse","82"="toulouse",
  "2A"="corse","2B"="corse")

d <- d %>%
  left_join(acad_com, by = "code_insee") %>%
  mutate(academie = coalesce(academie, dept2academie[dept]),
         academie = factor(academie))

region_academique <- c(
  "clermontferrand"="auvergne_rhone_alpes", "grenoble"="auvergne_rhone_alpes", "lyon"="auvergne_rhone_alpes",
  "besancon"="bourgogne_franche_comte", "dijon"="bourgogne_franche_comte",
  "rennes"="bretagne",
  "orleanstours"="centre_val_de_loire",
  "corse"="corse",
  "nancymetz"="grand_est", "reims"="grand_est", "strasbourg"="grand_est",
  "amiens"="hauts_de_france", "lille"="hauts_de_france",
  "creteil"="ile_de_france", "paris"="ile_de_france", "versailles"="ile_de_france",
  "normandie"="normandie",
  "bordeaux"="nouvelle_aquitaine", "limoges"="nouvelle_aquitaine", "poitiers"="nouvelle_aquitaine",
  "montpellier"="occitanie", "toulouse"="occitanie",
  "nantes"="pays_de_la_loire",
  "aixmarseille"="provence_alpes_cote_azur", "nice"="provence_alpes_cote_azur")

d <- d %>% mutate(region_acad = factor(region_academique[as.character(academie)]))
cat("Communes sans académie :", sum(is.na(d$academie)), "\n")
print(table(d$region_acad))

# ---------------------------------------------------------------------
# 3. VARIABLES
# ---------------------------------------------------------------------
d <- d %>%
  mutate(y       = .data[[VAR_Y]],
         log_pop = log(pop_totale),
         log_rev = log(revenu_moyen_com),
         z_theil = as.numeric(scale(theil_revenu)),
         z_moran = as.numeric(scale(moran_I_rev)))

rhs <- function(x) paste(x, collapse = " + ")

X_TAILLE  <- "log_pop"
X_CTRL    <- c("log_rev", "taux_pauv_commune", "nb_colleges_publics",
               "dist_college_public_km", "nb_prives_5km_autres")
X_SEG_BASE <- c("z_theil", "z_moran")
X_SEG_NEW  <- c("z_exposition_riches", "z_polarisation", "z_part_pop_top10", "z_segregation_multidim")

formules <- list(
  M0  = reformulate("z_theil", response = "y"),
  M1  = reformulate(c("z_theil", X_TAILLE), response = "y"),
  M2  = reformulate(c(X_SEG_BASE, X_TAILLE, X_CTRL), response = "y"),
  M3  = reformulate(c("z_theil", "z_exposition_riches", "z_polarisation",
                      "z_part_pop_top10", X_TAILLE, X_CTRL), response = "y"),
  M3b = reformulate(c("z_theil", "z_getis_ord_riches", "z_polarisation",
                      "z_part_pop_top10", X_TAILLE, X_CTRL), response = "y"),
  M3c = reformulate(c("z_segregation_multidim", "z_polarisation",
                      "z_part_pop_top10", X_TAILLE, X_CTRL), response = "y"),
  M3d = reformulate(c("z_theil", "z_entropie_diversite", "z_diversite_binaire",
                      X_TAILLE, X_CTRL), response = "y"),
  M4  = reformulate(c("z_theil", "z_exposition_riches", "z_polarisation",
                      "z_part_pop_top10", "z_segregation_multidim",
                      X_TAILLE, X_CTRL, "region_acad"), response = "y")
)

vars_utiles <- unique(c("y", "dept", "academie", "region_acad",
                        unlist(lapply(formules, all.vars))))
vars_utiles <- intersect(vars_utiles, names(d))
d_mod <- d[complete.cases(d[, vars_utiles]), ]
cat(sprintf("N = %d (%d perdus) | régions académiques : %d\n",
            nrow(d_mod), nrow(d) - nrow(d_mod), nlevels(droplevels(d_mod$region_acad))))

cat("\n--- Distribution de", VAR_Y, "---\n")
print(table(d_mod$y))
cat(sprintf("Moyenne = %.3f | Var/Moy = %.2f | Zéros = %.1f%%\n",
            mean(d_mod$y), var(d_mod$y)/mean(d_mod$y), 100*mean(d_mod$y == 0)))

# ---------------------------------------------------------------------
# 4. POISSON ROBUSTE
# ---------------------------------------------------------------------
fit_pois <- lapply(formules, function(f) glm(f, family = poisson, data = d_mod))

cluster_de <- function(nm) {
  if ("region_acad" %in% all.vars(formules[[nm]])) d_mod$region_acad else d_mod$dept
}
coef_robuste <- function(m, cluster) {
  V <- tryCatch(vcovCL(m, cluster = cluster, type = "HC1"), error = function(e) vcov(m))
  coeftest(m, vcov. = V)
}

cat("\n=== POISSON M3 (indicateurs ségrégation, cluster dept) ===\n")
print(coef_robuste(fit_pois$M3, d_mod$dept))
cat("\n=== POISSON M4 (FE région académique, cluster région) ===\n")
print(coef_robuste(fit_pois$M4, d_mod$region_acad))

# ---------------------------------------------------------------------
# 5. STABILITÉ DES COEFFICIENTS
# ---------------------------------------------------------------------
vars_stab <- c("z_theil", "z_moran", "z_exposition_riches", "z_polarisation",
               "z_part_pop_top10", "z_segregation_multidim")
get_row <- function(ct, v) if (v %in% rownames(ct)) ct[v, ] else c(NA, NA, NA, NA)
tab_stab <- do.call(rbind, lapply(names(fit_pois), function(nm) {
  ct <- coef_robuste(fit_pois[[nm]], cluster_de(nm))
  do.call(rbind, lapply(vars_stab, function(v) {
    r <- get_row(ct, v)
    data.frame(modele = nm, variable = v, coef = r[1], se = r[2], p = r[4])
  }))
}))
print(tab_stab, digits = 3)

# ---------------------------------------------------------------------
# 6. WALD JOINT FE RÉGION ACADÉMIQUE
# ---------------------------------------------------------------------
wt <- tryCatch({
  V <- vcovCL(fit_pois$M4, cluster = d_mod$region_acad, type = "HC1")
  betas <- coef(fit_pois$M4)[grepl("^region_acad", names(coef(fit_pois$M4)))]
  W <- as.numeric(t(betas) %*% solve(V[names(betas), names(betas)]) %*% betas)
  cat(sprintf("\nWald joint FE région académique : W = %.1f, df = %d, p = %.4g\n",
              W, length(betas), pchisq(W, length(betas), lower.tail = FALSE)))
}, error = function(e) message("Wald FE : ", e$message))

# ---------------------------------------------------------------------
# 7. NB2 (M4)
# ---------------------------------------------------------------------
m_p <- fit_pois$M4
cat("\nCameron-Trivedi NB2 :"); print(dispersiontest(m_p, trafo = 2))
m_nb <- tryCatch(MASS::glm.nb(formules$M4, data = d_mod, control = glm.control(maxit = 200)),
                 error = function(e) { message("NB2 non convergent : ", e$message); NULL })
if (!is.null(m_nb)) {
  cat("\n=== NB2 M4 ===\n")
  print(coef_robuste(m_nb, d_mod$region_acad))
}

# ---------------------------------------------------------------------
# 8. HURDLE / ZINB (sur M4)
# ---------------------------------------------------------------------
X_COUNT <- setdiff(all.vars(formules$M4), "y")
X_ZERO  <- X_COUNT
f_hurdle <- as.formula(paste("y ~", rhs(X_COUNT), "|", rhs(X_ZERO)))

m_hnb  <- tryCatch(hurdle(f_hurdle, data = d_mod, dist = "negbin", zero.dist = "binomial"),
                   error = function(e) { message("Hurdle NB : ", e$message); NULL })
m_zinb <- tryCatch(zeroinfl(f_hurdle, data = d_mod, dist = "negbin", link = "logit"),
                   error = function(e) { message("ZINB : ", e$message); NULL })
m_hp   <- tryCatch(hurdle(f_hurdle, data = d_mod, dist = "poisson", zero.dist = "binomial"),
                   error = function(e) { message("Hurdle Poisson : ", e$message); NULL })
if (!is.null(m_hnb)) { cat("\n=== HURDLE NB ===\n"); print(summary(m_hnb)) }
if (!is.null(m_zinb)) { cat("\n=== ZINB ===\n"); print(summary(m_zinb)) }

# ---------------------------------------------------------------------
# 9. COMPARAISON
# ---------------------------------------------------------------------
modeles <- Filter(Negate(is.null), list(Poisson = m_p, NB2 = m_nb,
                                        Hurdle_NB = m_hnb, ZINB = m_zinb, Hurdle_Pois = m_hp))
cat("\n--- AIC / BIC ---\n")
print(data.frame(modele = names(modeles), AIC = sapply(modeles, AIC),
                 BIC = sapply(modeles, BIC)), row.names = FALSE)

K <- min(max(d_mod$y), 8)
proba_pred <- function(m, K) {
  if (inherits(m, "negbin")) { mu <- fitted(m); sapply(0:K, function(j) sum(dnbinom(j, size = m$theta, mu = mu)))
  } else if (inherits(m, "glm")) { mu <- fitted(m); sapply(0:K, function(j) sum(dpois(j, mu)))
  } else colSums(predict(m, type = "prob", at = 0:K))
}
tab_freq <- data.frame(y = 0:K, Observe = sapply(0:K, function(j) sum(d_mod$y == j)),
                       sapply(modeles, proba_pred, K = K))
cat("\n--- Fréquences observées vs prédites ---\n")
print(round(tab_freq, 1), row.names = FALSE)

# ---------------------------------------------------------------------
# 10. IRR + AME (modèle préféré : Poisson M4)
# ---------------------------------------------------------------------
irr_table <- function(m, cluster = d_mod$region_acad) {
  ct <- coef_robuste(m, cluster)
  data.frame(variable = rownames(ct), IRR = exp(ct[,1]),
             IC_bas = exp(ct[,1] - 1.96*ct[,2]), IC_haut = exp(ct[,1] + 1.96*ct[,2]),
             p = ct[,4], row.names = NULL)
}
cat("\n--- IRR (Poisson M4) ---\n")
print(irr_table(m_p), digits = 3)

ame_num <- function(m, data, var, h = 0.01) {
  d1 <- data; d0 <- data
  d1[[var]] <- d1[[var]] + h/2; d0[[var]] <- d0[[var]] - h/2
  mean((predict(m, newdata = d1, type = "response") -
          predict(m, newdata = d0, type = "response")) / h)
}
cat("\n--- AME des indicateurs de ségrégation (Poisson M4) ---\n")
print(sapply(vars_stab, function(v) ame_num(m_p, d_mod, v)))

# ---------------------------------------------------------------------
# 11. SAUVEGARDE
# ---------------------------------------------------------------------
saveRDS(list(
  d_mod = d_mod, modeles = modeles, formules = formules,
  tab_stab = tab_stab, tab_freq = tab_freq, wald_fe = wt,
  ct_nb2 = dispersiontest(m_p, trafo = 2)
), file.path(DOSSIER_OUT, "modeles_comptage.rds"))
cat("modeles_comptage.rds écrit.\n")

