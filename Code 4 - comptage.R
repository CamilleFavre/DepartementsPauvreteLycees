# =====================================================================
# CODE 4 - Début de l'analyse économétrique : modèles de comptage
#
# Question : l'hétérogénéité sociale intra-communale est-elle associée
#            au nombre de lycées privés présents dans la commune ?
#
# Plan :
#  1. Données et statistiques descriptives du comptage
#  2. Poisson (PPML) avec écarts-types robustes clusterisés par département
#  3. Tests de surdispersion -> binomiale négative (NB2)
#  4. Excès de zéros : hurdle et zero-inflated (NB)
#  5. Comparaison : AIC/BIC, fréquences prédites, Vuong, validation croisée
#  6. Interprétation : IRR et effets marginaux moyens
#
# Packages à installer au besoin :
#   install.packages(c("arrow","dplyr","sandwich","lmtest","pscl","AER","MASS","ggplot2"))
# =====================================================================

library(arrow)
library(dplyr)
library(sandwich)
library(lmtest)
library(pscl)      # hurdle(), zeroinfl(), vuong()
library(AER)       # dispersiontest()
library(ggplot2)
# MASS est appelé via MASS::glm.nb (évite le conflit MASS::select / dplyr::select)

DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
VAR_Y       <- "nb_lycees_prives"     # alternatives : "nb_prives_gt", "nb_prives_sous_contrat"
set.seed(2025)

# ---------------------------------------------------------------------
# 1. DONNÉES
# ---------------------------------------------------------------------
d <- read_parquet(file.path(DOSSIER_OUT, "df_analyse.parquet")) %>%
  mutate(
    y       = .data[[VAR_Y]],
    log_pop = log(pop_totale),
    log_rev = log(revenu_moyen_com),
    # Hétérogénéité standardisée : le coefficient se lit "pour 1 écart-type de plus"
    z_theil = as.numeric(scale(theil_revenu)),
    z_moran = as.numeric(scale(moran_I_rev))
  )

# Spécifications (blocs de covariables)
X_HET   <- "z_theil"
X_TAILLE <- "log_pop"
X_CTRL  <- c("log_rev", "taux_pauv_commune", "nb_lycees_publics",
             "dist_lycee_public_gt_km", "nb_prives_10km_autres")
X_HET2  <- c("z_theil", "z_moran")      # deux dimensions : dispersion + structure spatiale

rhs <- function(x) paste(x, collapse = " + ")
formules <- list(
  M0 = reformulate(X_HET,                         response = "y"),
  M1 = reformulate(c(X_HET, X_TAILLE),            response = "y"),
  M2 = reformulate(c(X_HET, X_TAILLE, X_CTRL),    response = "y"),
  M3 = reformulate(c(X_HET2, X_TAILLE, X_CTRL),   response = "y")
)

# Échantillon commun à toutes les spécifications (comparabilité des vraisemblances)
vars_utiles <- unique(c("y", "dept", all.vars(formules$M3)))
d_mod <- d[complete.cases(d[, vars_utiles]), ]
cat(sprintf("N initial = %d | N estimation = %d (%d communes perdues pour NA)\n",
            nrow(d), nrow(d_mod), nrow(d) - nrow(d_mod)))
cat("Nombre de départements (clusters) :", length(unique(d_mod$dept)), "\n")

# ---------------------------------------------------------------------
# 1b. STATISTIQUES DESCRIPTIVES DU COMPTAGE
# ---------------------------------------------------------------------
cat("\n--- Distribution de", VAR_Y, "---\n")
print(table(d_mod$y))
cat(sprintf("Moyenne = %.3f | Variance = %.3f | Variance/Moyenne = %.2f | Zéros = %.1f%%\n",
            mean(d_mod$y), var(d_mod$y), var(d_mod$y) / mean(d_mod$y), 100 * mean(d_mod$y == 0)))
# Sous un Poisson, Var/Moyenne = 1 : un ratio >> 1 suggère surdispersion
# (attention : la surdispersion INCONDITIONNELLE peut venir seulement de l'hétérogénéité
#  des covariables ; les tests de la section 3 portent sur la surdispersion CONDITIONNELLE)

cat("\n--- Corrélations entre covariables (colinéarité) ---\n")
print(round(cor(d_mod[, c("z_theil", "z_moran", "log_pop", "log_rev",
                          "taux_pauv_commune", "nb_lycees_publics",
                          "dist_lycee_public_gt_km", "nb_prives_10km_autres")]), 2))

# ---------------------------------------------------------------------
# 2. POISSON / PPML
#    E[y|x] = exp(x'theta). Le Poisson est un PML : convergent dès que
#    l'espérance conditionnelle est bien spécifiée, même si y n'est pas
#    Poisson -> on utilise toujours un estimateur robuste de la variance
#    (sandwich), ici clusterisé par département (corrélation spatiale).
# ---------------------------------------------------------------------
fit_pois <- lapply(formules, function(f) glm(f, family = poisson(link = "log"), data = d_mod))

coef_robuste <- function(m, cluster = d_mod$dept) {
  V <- tryCatch(vcovCL(m, cluster = cluster, type = "HC1"),
                error = function(e) { warning("vcovCL a échoué : variance non robuste"); vcov(m) })
  coeftest(m, vcov. = V)
}

cat("\n=== POISSON M3 (écarts-types robustes, clusters = départements) ===\n")
print(coef_robuste(fit_pois$M3))

# Stabilité du coefficient d'hétérogénéité quand on change de spécification
cat("\n--- Coefficient de z_theil selon la spécification (Poisson, robuste) ---\n")
tab_stab <- do.call(rbind, lapply(names(fit_pois), function(nm) {
  ct <- coef_robuste(fit_pois[[nm]])
  data.frame(modele = nm, coef = ct["z_theil", 1], se = ct["z_theil", 2],
             p = ct["z_theil", 4])
}))
print(tab_stab, digits = 3)
# Biais attendu à discuter (variable omise) : un Theil élevé est corrélé à la taille et au revenu moyen
# -> M0 -> M1 -> M2 montre la direction et l'ampleur du biais de variable omise.

# Test : l'élasticité à la population est-elle égale à 1 (=> offset log(pop) justifié) ?
ct <- coef_robuste(fit_pois$M2)
z_pop <- (ct["log_pop", 1] - 1) / ct["log_pop", 2]
cat(sprintf("\nH0 : coef(log_pop) = 1 -> z = %.2f, p = %.3f\n", z_pop, 2 * pnorm(-abs(z_pop))))

# ---------------------------------------------------------------------
# 3. SURDISPERSION ET BINOMIALE NÉGATIVE (NB2)
#    Poisson : Var(y|x) = mu.   NB2 : Var(y|x) = mu (1 + mu/alpha).
#    (theta de MASS::glm.nb = alpha du cours ; alpha -> infini = Poisson)
# ---------------------------------------------------------------------
m_p <- fit_pois$M3
mu  <- fitted(m_p)

# (a) Test de Cameron-Trivedi fait "à la main" :
#     sous H0 (equidispersion), E[((y-mu)^2 - y)/mu | mu] = 0.
#     On régresse z_i = ((y-mu)^2 - y)/mu sur mu_i (sans constante) : t-test sur la pente.
#     Pente > 0 <=> Var = mu + c*mu^2 (NB2) ; on peut aussi tester c*mu (NB1).
aux <- lm(((d_mod$y - mu)^2 - d_mod$y) / mu ~ 0 + mu)
cat("\n--- Test de surdispersion de Cameron-Trivedi (forme NB2) ---\n")
print(summary(aux)$coefficients)

# (b) Même test via AER (trafo = 2 : alternative NB2 ; trafo = 1 : alternative NB1/quasi-Poisson)
cat("\nAER::dispersiontest, alternative NB2 :\n"); print(dispersiontest(m_p, trafo = 2))
cat("\nAER::dispersiontest, alternative NB1 :\n"); print(dispersiontest(m_p, trafo = 1))

# (c) Estimation NB2 et LR test Poisson vs NB
fit_nb <- lapply(formules, function(f) MASS::glm.nb(f, data = d_mod, control = glm.control(maxit = 100)))
m_nb <- fit_nb$M3
cat(sprintf("\nNB2 : alpha (theta) = %.3f (se = %.3f)\n", m_nb$theta, m_nb$SE.theta))

# H0 : 1/alpha = 0 se trouve sur le BORD de l'espace des paramètres ->
# la statistique LR suit un mélange 50:50 de chi2(0) et chi2(1) : on divise la p-value par 2
LR <- as.numeric(2 * (logLik(m_nb) - logLik(m_p)))
cat(sprintf("LR (NB vs Poisson) = %.2f | p-value (chi-bar-carré) = %.4g\n",
            LR, 0.5 * pchisq(LR, df = 1, lower.tail = FALSE)))

cat("\n=== NB2 M3 (robuste, clusters = départements) ===\n")
print(coef_robuste(m_nb))

# ---------------------------------------------------------------------
# 4. EXCÈS DE ZÉROS : HURDLE vs ZERO-INFLATED (binomiale négative)
#    Hurdle : P(y=0)=F(-f0(x)) ; y>0 suit une NB tronquée en 0.
#             (décision d'entrée du privé, puis combien)
#    ZI     : zéros "structurels" (contrainte) + zéros "d'échantillonnage" NB.
#    Dans pscl::hurdle la partie "zero" modélise P(y > 0) (logit).
#    Dans pscl::zeroinfl la partie "zero" modélise P(zéro structurel) (logit) :
#    les signes sont donc OPPOSÉS pour une même interprétation économique.
# ---------------------------------------------------------------------
f_count <- formules$M3
# Partie "zéro" : on met la taille, le revenu et l'hétérogénéité ; on peut en changer
X_ZERO  <- c(X_HET2, X_TAILLE, "log_rev", "taux_pauv_commune",
             "nb_lycees_publics", "nb_prives_10km_autres")
f_hurdle <- as.formula(paste("y ~", rhs(c(X_HET2, X_TAILLE, X_CTRL)), "|", rhs(X_ZERO)))

m_hnb <- tryCatch(hurdle(f_hurdle, data = d_mod, dist = "negbin", zero.dist = "binomial", link = "logit"),
                  error = function(e) { message("Hurdle NB non estimé : ", e$message); NULL })
m_zinb <- tryCatch(zeroinfl(f_hurdle, data = d_mod, dist = "negbin", link = "logit"),
                   error = function(e) { message("ZINB non estimé : ", e$message); NULL })

if (!is.null(m_hnb)) {
  cat("\n=== HURDLE NB ===\n")
  print(summary(m_hnb))
}
if (!is.null(m_zinb)) {
  cat("\n=== ZERO-INFLATED NB ===\n")
  print(summary(m_zinb))
}

# ---------------------------------------------------------------------
# 5. COMPARAISON DES MODÈLES
# ---------------------------------------------------------------------
modeles <- Filter(Negate(is.null), list(Poisson = m_p, NB2 = m_nb, Hurdle_NB = m_hnb, ZINB = m_zinb))

cat("\n--- AIC / BIC (même échantillon) ---\n")
print(data.frame(
  modele = names(modeles),
  logLik = sapply(modeles, function(m) as.numeric(logLik(m))),
  k      = sapply(modeles, function(m) attr(logLik(m), "df")),
  AIC    = sapply(modeles, AIC),
  BIC    = sapply(modeles, BIC)
), row.names = FALSE, digits = 6)

# Probabilités prédites P(y = j), j = 0..K, pour chaque modèle (moyennées sur l'échantillon)
K <- 6
proba_pred <- function(m, K) {
  if (inherits(m, "negbin")) {
    mu <- fitted(m); sapply(0:K, function(j) sum(dnbinom(j, size = m$theta, mu = mu)))
  } else if (inherits(m, "glm")) {
    mu <- fitted(m); sapply(0:K, function(j) sum(dpois(j, mu)))
  } else colSums(predict(m, type = "prob", at = 0:K))
}
obs <- sapply(0:K, function(j) sum(d_mod$y == j))
tab_freq <- data.frame(y = 0:K, Observe = obs,
                       sapply(modeles, proba_pred, K = K))
cat("\n--- Fréquences observées vs prédites ---\n")
print(round(tab_freq, 1), row.names = FALSE)

# Graphique : observé vs prédit
graph_df <- tidyr::pivot_longer(tab_freq, -y, names_to = "modele", values_to = "n")
p <- ggplot(graph_df, aes(x = factor(y), y = n, fill = modele)) +
  geom_col(position = "dodge") +
  labs(x = paste("Valeur de", VAR_Y), y = "Nombre de communes",
       title = "Fréquences observées et prédites par modèle") +
  theme_minimal()
ggsave(file.path(DOSSIER_OUT, "freq_observees_predites.png"), p, width = 8, height = 4.5, dpi = 150)

# Test de Vuong (modèles non emboîtés) - à interpréter avec prudence :
# NB vs ZINB n'est pas un cas régulier (paramètre de zéro-inflation sur le bord, cf. Wilson 2015)
if (!is.null(m_zinb)) {
  cat("\n--- Vuong : NB2 vs ZINB ---\n"); print(vuong(m_nb, m_zinb))
}
if (!is.null(m_hnb)) {
  cat("\n--- Vuong : NB2 vs Hurdle NB ---\n"); print(vuong(m_nb, m_hnb))
}

# ---------------------------------------------------------------------
# 5b. VALIDATION CROISÉE (K blocs) : performance HORS ÉCHANTILLON
#     Critères : RMSE, MAE et log-score prédictif moyen (plus haut = mieux).
#     Le log-score évalue toute la distribution prédite, pas seulement la moyenne.
# ---------------------------------------------------------------------
K_CV <- 10
folds <- sample(rep(seq_len(K_CV), length.out = nrow(d_mod)))

cv_fit_predict <- function(modele, train, test) {
  out <- tryCatch({
    if (modele == "Poisson") {
      m <- glm(f_count, family = poisson, data = train)
      mu <- predict(m, test, type = "response")
      lp <- dpois(test$y, mu, log = TRUE)
    } else if (modele == "NB2") {
      m <- MASS::glm.nb(f_count, data = train, control = glm.control(maxit = 100))
      mu <- predict(m, test, type = "response")
      lp <- dnbinom(test$y, size = m$theta, mu = mu, log = TRUE)
    } else {
      m <- if (modele == "Hurdle_NB")
        hurdle(f_hurdle, data = train, dist = "negbin", zero.dist = "binomial")
      else zeroinfl(f_hurdle, data = train, dist = "negbin")
      mu <- predict(m, test, type = "response")
      P  <- predict(m, test, type = "prob", at = 0:max(test$y))
      lp <- log(P[cbind(seq_len(nrow(test)), test$y + 1)])
    }
    data.frame(y = test$y, mu = as.numeric(mu), lp = as.numeric(lp))
  }, error = function(e) NULL)
  out
}

noms_cv <- intersect(c("Poisson", "NB2", "Hurdle_NB", "ZINB"), names(modeles))
cv_res <- lapply(noms_cv, function(nm) {
  preds <- lapply(seq_len(K_CV), function(k)
    cv_fit_predict(nm, d_mod[folds != k, ], d_mod[folds == k, ]))
  echecs <- sum(sapply(preds, is.null))
  P <- bind_rows(preds)
  data.frame(modele = nm,
             RMSE = sqrt(mean((P$y - P$mu)^2)),
             MAE  = mean(abs(P$y - P$mu)),
             logscore = mean(P$lp[is.finite(P$lp)]),
             blocs_echoues = echecs)
})
cat("\n--- Validation croisée (", K_CV, " blocs) ---\n", sep = "")
print(bind_rows(cv_res), row.names = FALSE, digits = 4)

# ---------------------------------------------------------------------
# 6. INTERPRÉTATION : IRR ET EFFETS MARGINAUX MOYENS
# ---------------------------------------------------------------------
# 6a. Incidence Rate Ratios (exp(coef)) pour Poisson/NB, IC robustes
irr_table <- function(m) {
  ct <- coef_robuste(m)
  data.frame(variable = rownames(ct),
             IRR = exp(ct[, 1]),
             IC_bas  = exp(ct[, 1] - 1.96 * ct[, 2]),
             IC_haut = exp(ct[, 1] + 1.96 * ct[, 2]),
             p = ct[, 4], row.names = NULL)
}
cat("\n--- IRR (NB2 M3) : effet multiplicatif sur le nombre attendu de lycées privés ---\n")
print(irr_table(m_nb), digits = 3)

# 6b. Effet marginal moyen par dérivée numérique de E[y|x] : valable pour TOUS les modèles
#     (dans le hurdle/ZI, E[y|x] combine les deux parties : l'effet marginal total
#      mêle marge extensive (entrée) et marge intensive (nombre)).
ame_num <- function(m, data, var, h = 0.01) {
  d1 <- data; d0 <- data
  d1[[var]] <- d1[[var]] + h / 2; d0[[var]] <- d0[[var]] - h / 2
  mean((predict(m, newdata = d1, type = "response") -
          predict(m, newdata = d0, type = "response")) / h)
}
cat("\n--- Effet marginal moyen d'1 écart-type de plus de z_theil sur E[y] ---\n")
ame_tab <- sapply(modeles, function(m) ame_num(m, d_mod, "z_theil"))
print(round(ame_tab, 4))
# Si ces valeurs sont proches entre modèles : le résultat ne dépend pas du choix de la distribution.
# Sinon : interpréter la différence en termes de biais (mauvaise spécification de la variance / des zéros).

# ---------------------------------------------------------------------
# 7. SUITE À PRÉVOIR (idées, non codées ici)
#  - Bootstrap par grappes (départements) des effets marginaux
#  - Effets aléatoires par département : MASS::glmmPQL / glmmTMB (family = nbinom2)
#    ou effets fixes de région (Mundlak) -> chapitre "panels multi-niveaux"
#  - Corrélation spatiale résiduelle : Moran des résidus de Pearson
#  - Alternatives pour la variable dépendante : nb_prives_gt, nb_prives_sous_contrat, part_eff_prive_gt
#  - Simulation de Monte Carlo : biais des estimateurs sous erreur de mesure du Theil
#  - Non-linéarité : terme quadratique en z_theil, interaction z_theil x z_moran
# ---------------------------------------------------------------------