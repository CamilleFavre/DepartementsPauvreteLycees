# =====================================================================
# CODE 5 - Diagnostics et robustesse (suite du Code 4)
#
# Ce script répond aux points ouverts de l'interprétation du Code 4 :
#   1. Quelles sont les 42 communes perdues ? (sélection ?)
#   2. Convergence du ZINB / du hurdle ; hurdle Poisson vs NB ; séparation
#   3. Écarts-types clusterisés pour le hurdle
#   4. Validation croisée : explication du RMSE anormal du NB2 (plusieurs graines)
#   5. Colinéarité : VIF, tests de Wald conjoints, ajout/retrait variable par variable
#   6. Robustesse : indicateurs alternatifs, non-linéarité, autres variables
#      dépendantes, offset, effet aléatoire département, observations influentes,
#      bootstrap par départements des effets marginaux
#
# TOUTES les sorties sont enregistrées dans diagnostics_code5.txt :
# c'est ce fichier qu'il faut me renvoyer.
# (si le script plante, tape sink() dans la console pour refermer le fichier)
# =====================================================================

library(arrow)
library(dplyr)
library(sandwich)
library(lmtest)
library(pscl)

DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
VAR_Y       <- "nb_lycees_prives"
N_SEEDS_CV  <- 5       # nb de graines pour la validation croisée
RUN_BOOT    <- TRUE    # bootstrap par départements (quelques minutes)
B_BOOT      <- 199

while (sink.number() > 0) sink()
sink(file.path(DOSSIER_OUT, "diagnostics_code5.txt"), split = TRUE)
set.seed(2025)
options(width = 140)

titre <- function(x) cat("\n\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

# ---------------------------------------------------------------------
# 0. DONNÉES (identique au Code 4)
# ---------------------------------------------------------------------
d <- read_parquet(file.path(DOSSIER_OUT, "df_analyse.parquet")) %>%
  mutate(
    y       = .data[[VAR_Y]],
    log_pop = log(pop_totale),
    log_rev = log(revenu_moyen_com),
    z_theil = as.numeric(scale(theil_revenu)),
    z_moran = as.numeric(scale(moran_I_rev))
  )

X_HET    <- "z_theil"
X_TAILLE <- "log_pop"
X_CTRL   <- c("log_rev", "taux_pauv_commune", "nb_lycees_publics",
              "dist_lycee_public_gt_km", "nb_prives_10km_autres")
X_HET2   <- c("z_theil", "z_moran")
X_ZERO   <- c(X_HET2, X_TAILLE, "log_rev", "taux_pauv_commune",
              "nb_lycees_publics", "nb_prives_10km_autres")
rhs <- function(x) paste(x, collapse = " + ")

formules <- list(
  M0 = reformulate(X_HET, response = "y"),
  M1 = reformulate(c(X_HET, X_TAILLE), response = "y"),
  M2 = reformulate(c(X_HET, X_TAILLE, X_CTRL), response = "y"),
  M3 = reformulate(c(X_HET2, X_TAILLE, X_CTRL), response = "y")
)
vars_utiles <- unique(c("y", "dept", all.vars(formules$M3)))

# ---------------------------------------------------------------------
# 1. LES COMMUNES PERDUES (sélection ?)
# ---------------------------------------------------------------------
titre("1. COMMUNES PERDUES POUR VALEURS MANQUANTES")

ok <- complete.cases(d[, vars_utiles])
cat("Nb de NA par variable (échantillon complet, N =", nrow(d), ") :\n")
print(colSums(is.na(d[, vars_utiles])))

dropped <- d[!ok, ]
dropped$manquants <- apply(is.na(dropped[, vars_utiles]), 1,
                           function(r) paste(vars_utiles[r], collapse = ","))
cat("\nVariable(s) manquante(s) pour les communes perdues :\n")
print(table(dropped$manquants))
cat("\nListe des communes perdues (triées par population) :\n")
print(as.data.frame(dropped %>% arrange(desc(pop_totale)) %>%
                      select(code_insee, nom_commune, pop_totale, y, manquants)), row.names = FALSE)

cat("\nComparaison communes retenues / perdues :\n")
comp <- d %>%
  mutate(groupe = ifelse(ok, "retenues", "perdues")) %>%
  group_by(groupe) %>%
  summarise(n = n(),
            pop_mediane = median(pop_totale),
            part_pop_100k = mean(pop_totale > 100000),
            y_moyen = mean(y, na.rm = TRUE),
            theil_moyen = mean(theil_revenu, na.rm = TRUE),
            revenu_moyen = mean(revenu_moyen_com, na.rm = TRUE))
print(as.data.frame(comp))

# Échantillon d'estimation ; variables continues CENTRÉES (ne change que les constantes,
# améliore la convergence des modèles à deux parties et rend les constantes interprétables)
d_mod <- d[ok, ]
d_mod <- d_mod %>% mutate(across(c(log_pop, log_rev, taux_pauv_commune), ~ . - mean(.)))
cat("\nN estimation :", nrow(d_mod), "| clusters (départements) :", length(unique(d_mod$dept)), "\n")

# ---------------------------------------------------------------------
# Fonctions utilitaires
# ---------------------------------------------------------------------
fit_cl <- function(formula, data, family = poisson(link = "log")) {
  vars <- unique(c(all.vars(formula), "dept"))
  dd <- data[complete.cases(data[, vars]), ]
  m <- glm(formula, family = family, data = dd)
  V <- tryCatch(vcovCL(m, cluster = dd$dept, type = "HC1"), error = function(e) vcov(m))
  list(m = m, V = V, n = nrow(dd), ct = coeftest(m, vcov. = V))
}
coef_row <- function(fit, term) {
  if (!term %in% rownames(fit$ct)) return(c(coef = NA, se = NA, p = NA))
  r <- fit$ct[term, c(1, 2, 4)]; names(r) <- c("coef", "se", "p"); r
}
wald_joint <- function(m, V, terms) {
  b <- coef(m)[terms]; Vs <- V[terms, terms, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(Vs) %*% b)
  c(W = W, df = length(terms), p = pchisq(W, length(terms), lower.tail = FALSE))
}
ame_num <- function(m, data, var, h = 0.01) {
  d1 <- data; d0 <- data
  d1[[var]] <- d1[[var]] + h / 2; d0[[var]] <- d0[[var]] - h / 2
  mean((predict(m, newdata = d1, type = "response") -
          predict(m, newdata = d0, type = "response")) / h)
}

f_count  <- formules$M3
f_hurdle <- as.formula(paste("y ~", rhs(c(X_HET2, X_TAILLE, X_CTRL)), "|", rhs(X_ZERO)))

# ---------------------------------------------------------------------
# 2. CONVERGENCE, HURDLE POISSON vs NB, SÉPARATION
# ---------------------------------------------------------------------
titre("2. CONVERGENCE ET SPÉCIFICATION DES MODÈLES À DEUX PARTIES")

m_hnb  <- tryCatch(hurdle(f_hurdle, data = d_mod, dist = "negbin", zero.dist = "binomial"),
                   error = function(e) { message(e$message); NULL })
m_hp   <- tryCatch(hurdle(f_hurdle, data = d_mod, dist = "poisson", zero.dist = "binomial"),
                   error = function(e) { message(e$message); NULL })
m_zinb <- tryCatch(zeroinfl(f_hurdle, data = d_mod, dist = "negbin",
                            control = zeroinfl.control(maxit = 5000, reltol = 1e-12)),
                   error = function(e) { message(e$message); NULL })
m_zip  <- tryCatch(zeroinfl(f_hurdle, data = d_mod, dist = "poisson"),
                   error = function(e) { message(e$message); NULL })

cat("Convergence (TRUE = ok) :\n")
cat(" Hurdle NB   :", if (!is.null(m_hnb))  m_hnb$converged  else NA, "\n")
cat(" Hurdle Pois :", if (!is.null(m_hp))   m_hp$converged   else NA, "\n")
cat(" ZINB        :", if (!is.null(m_zinb)) m_zinb$converged else NA,
    "| code optim :", if (!is.null(m_zinb)) m_zinb$optim$convergence else NA, "(0 = ok)\n")
cat(" ZIP         :", if (!is.null(m_zip))  m_zip$converged   else NA, "\n")

if (!is.null(m_hnb) && !is.null(m_hp)) {
  LR <- 2 * (as.numeric(logLik(m_hnb)) - as.numeric(logLik(m_hp)))
  cat(sprintf("\nHurdle NB vs Hurdle Poisson : LR = %.2f, p (chi-bar-carré) = %.4f\n",
              LR, 0.5 * pchisq(max(LR, 0), 1, lower.tail = FALSE)))
  cat("theta du hurdle NB :", m_hnb$theta, "(theta grand = quasi-Poisson)\n")
}
cat("\nAIC / BIC :\n")
mods2 <- Filter(Negate(is.null), list(Hurdle_NB = m_hnb, Hurdle_Pois = m_hp, ZINB = m_zinb, ZIP = m_zip))
print(data.frame(modele = names(mods2), logLik = sapply(mods2, function(m) as.numeric(logLik(m))),
                 AIC = sapply(mods2, AIC), BIC = sapply(mods2, BIC)), row.names = FALSE)

# Séparation (quasi-parfaite) dans la partie "présence" du hurdle
cat("\n--- Partie 'au moins un lycée privé' : logit simple, diagnostic de séparation ---\n")
d_mod$pos <- as.integer(d_mod$y > 0)
lg <- suppressWarnings(glm(reformulate(X_ZERO, response = "pos"), family = binomial, data = d_mod))
cat("Convergence du logit :", lg$converged, "\n")
cat("Probabilités ajustées : min =", round(min(fitted(lg)), 4), "| max =", round(max(fitted(lg)), 4), "\n")
cat("Nb de probabilités < 0.01 :", sum(fitted(lg) < 0.01), "| > 0.99 :", sum(fitted(lg) > 0.99), "\n")
cat("Part de communes avec au moins un lycée privé :", round(mean(d_mod$pos), 3), "\n")

# ---------------------------------------------------------------------
# 3. ÉCARTS-TYPES ROBUSTES (CLUSTERISÉS) POUR LE HURDLE
# ---------------------------------------------------------------------
titre("3. HURDLE : ÉCARTS-TYPES CLUSTERISÉS PAR DÉPARTEMENT")
for (nm in c("Hurdle_NB", "Hurdle_Pois")) {
  m <- mods2[[nm]]
  if (is.null(m)) next
  cat("\n---", nm, "---\n")
  V <- tryCatch(vcovCL(m, cluster = d_mod$dept), error = function(e) { cat("vcovCL a échoué :", e$message, "\n"); NULL })
  if (!is.null(V)) print(coeftest(m, vcov. = V))
}

# ---------------------------------------------------------------------
# 4. VALIDATION CROISÉE REPRODUITE SUR PLUSIEURS GRAINES
# ---------------------------------------------------------------------
titre("4. VALIDATION CROISÉE (plusieurs graines) ET DIAGNOSTIC DU RMSE DU NB2")

# Version tolérante aux avertissements (un bloc n'est écarté qu'en cas d'ERREUR)
cv_pred_safe <- function(modele, train, test) suppressWarnings(
  tryCatch({
    if (modele == "Poisson") {
      m <- glm(f_count, family = poisson, data = train)
      mu <- predict(m, test, type = "response"); lp <- dpois(test$y, mu, log = TRUE)
    } else if (modele == "NB2") {
      m <- MASS::glm.nb(f_count, data = train, control = glm.control(maxit = 100))
      mu <- predict(m, test, type = "response"); lp <- dnbinom(test$y, size = m$theta, mu = mu, log = TRUE)
    } else {
      m <- switch(modele,
                  Hurdle_NB   = hurdle(f_hurdle, data = train, dist = "negbin",  zero.dist = "binomial"),
                  Hurdle_Pois = hurdle(f_hurdle, data = train, dist = "poisson", zero.dist = "binomial"),
                  ZINB        = zeroinfl(f_hurdle, data = train, dist = "negbin"))
      mu <- predict(m, test, type = "response")
      P  <- predict(m, test, type = "prob", at = 0:max(test$y))
      lp <- log(P[cbind(seq_len(nrow(test)), test$y + 1)])
    }
    data.frame(id = test$.id, y = test$y, mu = as.numeric(mu), lp = as.numeric(lp))
  }, error = function(e) NULL))

d_cv <- d_mod; d_cv$.id <- seq_len(nrow(d_cv))
noms_cv <- c("Poisson", "NB2", "Hurdle_NB", "Hurdle_Pois", "ZINB")
res_cv <- list(); toutes_preds <- list()
for (s in seq_len(N_SEEDS_CV)) {
  set.seed(1000 + s)
  folds <- sample(rep(1:10, length.out = nrow(d_cv)))
  for (nm in noms_cv) {
    preds <- lapply(1:10, function(k) cv_pred_safe(nm, d_cv[folds != k, ], d_cv[folds == k, ]))
    P <- bind_rows(preds)
    toutes_preds[[paste(nm, s)]] <- cbind(modele = nm, graine = s, P)
    res_cv[[paste(nm, s)]] <- data.frame(
      modele = nm, graine = s,
      RMSE = sqrt(mean((P$y - P$mu)^2)),
      RMSE_median_err = median(abs(P$y - P$mu)),
      MAE = mean(abs(P$y - P$mu)),
      logscore = mean(P$lp[is.finite(P$lp)]),
      mu_max = max(P$mu),
      blocs_echoues = sum(sapply(preds, is.null)))
  }
}
res_cv <- bind_rows(res_cv)
cat("Résultats par graine :\n"); print(res_cv, row.names = FALSE, digits = 4)
cat("\nMoyennes sur les graines :\n")
print(res_cv %>% group_by(modele) %>%
        summarise(RMSE = mean(RMSE), MAE = mean(MAE), logscore = mean(logscore),
                  logscore_sd = sd(logscore), mu_max = max(mu_max)), digits = 4)

cat("\nLes 8 pires prédictions du NB2 (toutes graines) :\n")
pire <- bind_rows(toutes_preds) %>% filter(modele == "NB2") %>%
  mutate(err = abs(y - mu)) %>% arrange(desc(err)) %>% head(8)
pire$commune <- d_mod$nom_commune[pire$id]; pire$pop <- d_mod$pop_totale[pire$id]
print(pire[, c("graine", "commune", "pop", "y", "mu", "err")], row.names = FALSE, digits = 4)

# ---------------------------------------------------------------------
# 5. COLINÉARITÉ ET STABILITÉ DU COEFFICIENT D'HÉTÉROGÉNÉITÉ
# ---------------------------------------------------------------------
titre("5. COLINÉARITÉ ET STABILITÉ DE z_theil")

covs <- c(X_HET2, X_TAILLE, X_CTRL)
vif <- sapply(covs, function(v)
  1 / (1 - summary(lm(reformulate(setdiff(covs, v), response = v), data = d_mod))$r.squared))
cat("VIF (règle empirique : > 5 préoccupant, > 10 grave) :\n"); print(round(vif, 2))

fp3 <- fit_cl(formules$M3, d_mod)
fn3 <- {
  m <- MASS::glm.nb(formules$M3, data = d_mod, control = glm.control(maxit = 100))
  list(m = m, V = vcovCL(m, cluster = d_mod$dept, type = "HC1"))
}
cat("\nTests de Wald conjoints (variance clusterisée) :\n")
tests <- list(
  "z_theil = z_moran = 0 (hétérogénéité)"     = c("z_theil", "z_moran"),
  "log_rev = taux_pauv = 0"                   = c("log_rev", "taux_pauv_commune"),
  "lycées publics et distance publique = 0"   = c("nb_lycees_publics", "dist_lycee_public_gt_km"))
for (nm in names(tests)) {
  wp <- wald_joint(fp3$m, fp3$V, tests[[nm]]); wn <- wald_joint(fn3$m, fn3$V, tests[[nm]])
  cat(sprintf("  %-45s Poisson : W=%.2f p=%.4f | NB2 : W=%.2f p=%.4f\n",
              nm, wp["W"], wp["p"], wn["W"], wn["p"]))
}

cat("\n--- Coefficient de z_theil : AJOUT d'une variable à M1 (z_theil + log_pop) ---\n")
base <- c(X_HET, X_TAILLE)
ajout <- do.call(rbind, lapply(c(X_CTRL, "z_moran"), function(v) {
  f <- fit_cl(reformulate(c(base, v), response = "y"), d_mod)
  data.frame(variable_ajoutee = v, t(coef_row(f, "z_theil")))
}))
print(ajout, row.names = FALSE, digits = 3)

cat("\n--- Coefficient de z_theil : RETRAIT d'une variable de M3 ---\n")
retrait <- do.call(rbind, lapply(setdiff(covs, "z_theil"), function(v) {
  f <- fit_cl(reformulate(setdiff(covs, v), response = "y"), d_mod)
  data.frame(variable_retiree = v, t(coef_row(f, "z_theil")))
}))
print(retrait, row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------
# 6. ROBUSTESSE
# ---------------------------------------------------------------------
titre("6a. INDICATEURS D'HÉTÉROGÉNÉITÉ ALTERNATIFS (un à la fois, Poisson robuste)")
indicateurs <- intersect(c("theil_revenu", "gini_revenu", "cv_revenu", "ratio_p90_p10_rev",
                           "dissim_pauv", "moran_I_rev", "contraste_rel_rev", "contraste_rel_pauv",
                           "part_mixte_rev", "part_segreg_rev", "part_mixte_pauv", "part_segreg_pauv"),
                         names(d_mod))
alt <- do.call(rbind, lapply(indicateurs, function(v) {
  dd <- d_mod; dd$z_ind <- as.numeric(scale(dd[[v]]))
  f <- fit_cl(reformulate(c("z_ind", X_TAILLE, X_CTRL), response = "y"), dd)
  data.frame(indicateur = v, N = f$n, t(coef_row(f, "z_ind")))
}))
alt$p_holm <- p.adjust(alt$p, method = "holm")
print(alt, row.names = FALSE, digits = 3)
cat("(p_holm : correction de Holm pour tests multiples sur la famille d'indicateurs)\n")

titre("6b. NON-LINÉARITÉ ET INTERACTION")
dd <- d_mod
f_quad <- fit_cl(update(formules$M3, . ~ . + I(z_theil^2)), dd)
f_int  <- fit_cl(update(formules$M3, . ~ . + z_theil:z_moran), dd)
cat("Terme quadratique :\n"); print(f_quad$ct[c("z_theil", "I(z_theil^2)"), ], digits = 3)
cat("\nInteraction z_theil x z_moran :\n"); print(f_int$ct[c("z_theil", "z_moran", "z_theil:z_moran"), ], digits = 3)

titre("6c. AUTRES VARIABLES DÉPENDANTES (coefficients d'hétérogénéité, Poisson robuste)")
dvs <- intersect(c("nb_lycees_prives", "nb_prives_gt", "nb_prives_pro_seul",
                   "nb_prives_sous_contrat", "nb_prives_hors_contrat"), names(d_mod))
for (v in dvs) {
  cat(sprintf("\n%s : moyenne = %.2f, zéros = %.1f%%\n", v, mean(d_mod[[v]]), 100 * mean(d_mod[[v]] == 0)))
  f <- fit_cl(reformulate(c(X_HET2, X_TAILLE, X_CTRL), response = v), d_mod)
  print(f$ct[X_HET2, ], digits = 3)
}
if ("part_eff_prive_gt" %in% names(d_mod)) {
  cat("\npart_eff_prive_gt (part du privé dans les effectifs GT) : quasi-binomial, robuste\n")
  f <- tryCatch(fit_cl(reformulate(c(X_HET2, X_TAILLE, X_CTRL), response = "part_eff_prive_gt"),
                       d_mod, family = quasibinomial(link = "logit")), error = function(e) NULL)
  if (!is.null(f)) { cat("N =", f$n, "\n"); print(f$ct[X_HET2, ], digits = 3) }
}

titre("6d. OFFSET log(pop) À LA PLACE DU RÉGRESSEUR log_pop")
f_off <- fit_cl(reformulate(c(X_HET2, X_CTRL, "offset(log_pop)"), response = "y"), d_mod)
print(f_off$ct[X_HET2, ], digits = 3)
cat("(à comparer avec le Poisson M3 sans offset :)\n"); print(fp3$ct[X_HET2, ], digits = 3)

titre("6e. DÉPENDANCE INTRA-DÉPARTEMENT")
res_p <- residuals(fp3$m, type = "pearson")
an <- anova(lm(res_p ~ factor(d_mod$dept)))
k_moy <- nrow(d_mod) / length(unique(d_mod$dept))
msb <- an[1, "Mean Sq"]; msw <- an[2, "Mean Sq"]
cat(sprintf("Corrélation intra-classe des résidus de Pearson (par département) ~ %.3f (taille moyenne de grappe %.1f)\n",
            (msb - msw) / (msb + (k_moy - 1) * msw), k_moy))
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  cat("\nNB2 avec effet aléatoire département (glmmTMB) :\n")
  m_re <- tryCatch(glmmTMB::glmmTMB(update(formules$M3, . ~ . + (1 | dept)),
                                    family = glmmTMB::nbinom2, data = d_mod), error = function(e) NULL)
  if (!is.null(m_re)) print(summary(m_re))
} else cat("\n(glmmTMB non installé : install.packages('glmmTMB') pour l'effet aléatoire)\n")

titre("6f. OBSERVATIONS INFLUENTES")
cd <- cooks.distance(fp3$m)
top <- order(cd, decreasing = TRUE)[1:10]
print(data.frame(commune = d_mod$nom_commune[top], pop = d_mod$pop_totale[top],
                 y = d_mod$y[top], cooks = round(cd[top], 3)), row.names = FALSE)
f_sans <- fit_cl(formules$M3, d_mod[-top, ])
cat("\nCoefficients d'hétérogénéité après retrait des 10 communes les plus influentes (N =", f_sans$n, ") :\n")
print(f_sans$ct[X_HET2, ], digits = 3)
cat("Pour mémoire, échantillon complet :\n"); print(fp3$ct[X_HET2, ], digits = 3)

titre("6g. BOOTSTRAP PAR DÉPARTEMENTS DE L'EFFET MARGINAL MOYEN DE z_theil")
if (RUN_BOOT) {
  depts <- unique(d_mod$dept)
  res <- matrix(NA_real_, B_BOOT, 2, dimnames = list(NULL, c("Poisson", "NB2")))
  for (b in seq_len(B_BOOT)) {
    s <- sample(depts, length(depts), replace = TRUE)
    db <- bind_rows(lapply(s, function(g) d_mod[d_mod$dept == g, ]))
    mp <- suppressWarnings(tryCatch(glm(f_count, family = poisson, data = db), error = function(e) NULL))
    if (!is.null(mp)) res[b, 1] <- ame_num(mp, db, "z_theil")
    mn <- suppressWarnings(tryCatch(MASS::glm.nb(f_count, data = db, control = glm.control(maxit = 100)),
                                    error = function(e) NULL))
    if (!is.null(mn)) res[b, 2] <- ame_num(mn, db, "z_theil")
  }
  point <- c(Poisson = ame_num(fp3$m, d_mod, "z_theil"), NB2 = ame_num(fn3$m, d_mod, "z_theil"))
  cat("Effet marginal moyen (lycées privés en plus pour +1 écart-type de Theil) :\n")
  print(data.frame(modele = names(point), AME = round(point, 4),
                   IC95_bas = round(apply(res, 2, quantile, 0.025, na.rm = TRUE), 4),
                   IC95_haut = round(apply(res, 2, quantile, 0.975, na.rm = TRUE), 4),
                   replications_valides = colSums(!is.na(res)), row.names = NULL))
} else cat("(désactivé : RUN_BOOT = FALSE)\n")

titre("INFORMATIONS DE SESSION")
print(sessionInfo())
sink()
cat("\nTerminé. Fichier de sortie :", file.path(DOSSIER_OUT, "diagnostics_code5.txt"), "\n")