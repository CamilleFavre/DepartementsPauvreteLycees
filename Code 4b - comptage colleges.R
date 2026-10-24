# =====================================================================
# CODE 4b - COLLÈGES : Poisson -> NB2 -> hurdle/ZI
# Identique au Code 4 lycées, seules les variables changent.
# =====================================================================
library(arrow); library(dplyr); library(sandwich); library(lmtest)
library(pscl); library(AER); library(ggplot2)

DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
VAR_Y       <- "nb_colleges_prives"
set.seed(2025)

d <- read_parquet(file.path(DOSSIER_OUT, "df_analyse_colleges.parquet")) %>%
  mutate(y       = .data[[VAR_Y]],
         log_pop = log(pop_totale),
         log_rev = log(revenu_moyen_com),
         z_theil = as.numeric(scale(theil_revenu)),
         z_moran = as.numeric(scale(moran_I_rev)))

X_HET    <- "z_theil"
X_TAILLE <- "log_pop"
X_CTRL   <- c("log_rev", "taux_pauv_commune", "nb_colleges_publics",
              "dist_college_public_km", "nb_prives_5km_autres")
X_HET2   <- c("z_theil", "z_moran")

rhs <- function(x) paste(x, collapse = " + ")
formules <- list(
  M0 = reformulate(X_HET, response = "y"),
  M1 = reformulate(c(X_HET, X_TAILLE), response = "y"),
  M2 = reformulate(c(X_HET, X_TAILLE, X_CTRL), response = "y"),
  M3 = reformulate(c(X_HET2, X_TAILLE, X_CTRL), response = "y"))

vars_utiles <- unique(c("y", "dept", all.vars(formules$M3)))
d_mod <- d[complete.cases(d[, vars_utiles]), ]
cat(sprintf("N = %d (%d perdus)\n", nrow(d_mod), nrow(d) - nrow(d_mod)))

cat("\n--- Distribution de", VAR_Y, "---\n")
print(table(d_mod$y))
cat(sprintf("Moyenne = %.3f | Var/Moy = %.2f | Zéros = %.1f%%\n",
            mean(d_mod$y), var(d_mod$y)/mean(d_mod$y), 100*mean(d_mod$y == 0)))

# ---- Poisson robuste (clusters départements) ----
fit_pois <- lapply(formules, function(f) glm(f, family = poisson, data = d_mod))
coef_robuste <- function(m, cluster = d_mod$dept) {
  V <- tryCatch(vcovCL(m, cluster = cluster, type = "HC1"),
                error = function(e) vcov(m))
  coeftest(m, vcov. = V)
}
cat("\n=== POISSON M3 ===\n"); print(coef_robuste(fit_pois$M3))

# Stabilité de z_theil et z_moran selon la spécification
get_row <- function(ct, v) if (v %in% rownames(ct)) ct[v, ] else c(NA, NA, NA, NA)

tab_stab <- do.call(rbind, lapply(names(fit_pois), function(nm) {
  ct <- coef_robuste(fit_pois[[nm]])
  rbind(data.frame(modele = nm, var = "z_theil",
                   coef = get_row(ct,"z_theil")[1], se = get_row(ct,"z_theil")[2],
                   p = get_row(ct,"z_theil")[4]),
        data.frame(modele = nm, var = "z_moran",
                   coef = get_row(ct,"z_moran")[1], se = get_row(ct,"z_moran")[2],
                   p = get_row(ct,"z_moran")[4]))
}))
print(tab_stab, digits = 3)

# ---- Surdispersion -> NB2 ----
m_p <- fit_pois$M3
cat("\nCameron-Trivedi NB2 :"); print(dispersiontest(m_p, trafo = 2))
m_nb <- MASS::glm.nb(formules$M3, data = d_mod, control = glm.control(maxit = 100))
cat("\n=== NB2 M3 ===\n"); print(coef_robuste(m_nb))

# ---- Hurdle NB (marge d'implantation + intensité) ----
X_ZERO <- c(X_HET2, X_TAILLE, "log_rev", "taux_pauv_commune",
            "nb_colleges_publics", "nb_prives_5km_autres")
f_hurdle <- as.formula(paste("y ~", rhs(c(X_HET2, X_TAILLE, X_CTRL)), "|", rhs(X_ZERO)))

m_hnb <- tryCatch(hurdle(f_hurdle, data = d_mod, dist = "negbin",
                         zero.dist = "binomial", link = "logit"),
                  error = function(e) { message(e$message); NULL })
m_zinb <- tryCatch(zeroinfl(f_hurdle, data = d_mod, dist = "negbin", link = "logit"),
                   error = function(e) { message(e$message); NULL })
if (!is.null(m_hnb)) { cat("\n=== HURDLE NB ===\n"); print(summary(m_hnb)) }
if (!is.null(m_zinb)) { cat("\n=== ZINB ===\n");     print(summary(m_zinb)) }

m_hp <- hurdle(f_hurdle, data = d_mod, dist = "poisson",
               zero.dist = "binomial", link = "logit")
summary(m_hp)


# ---- Comparaison ----
modeles <- Filter(Negate(is.null), list(Poisson = m_p, NB2 = m_nb,
                                        Hurdle_NB = m_hnb, ZINB = m_zinb))
cat("\n--- AIC / BIC ---\n")
print(data.frame(modele = names(modeles),
                 AIC = sapply(modeles, AIC), BIC = sapply(modeles, BIC)),
      row.names = FALSE)

# Fréquences observées vs prédites
K <- max(d_mod$y); K <- min(K, 8)
proba_pred <- function(m, K) {
  if (inherits(m, "negbin")) { mu <- fitted(m)
  sapply(0:K, function(j) sum(dnbinom(j, size = m$theta, mu = mu)))
  } else if (inherits(m, "glm")) { mu <- fitted(m)
  sapply(0:K, function(j) sum(dpois(j, mu)))
  } else colSums(predict(m, type = "prob", at = 0:K))
}
tab_freq <- data.frame(y = 0:K, Observe = sapply(0:K, function(j) sum(d_mod$y == j)),
                       sapply(modeles, proba_pred, K = K))
cat("\n--- Fréquences observées vs prédites ---\n")
print(round(tab_freq, 1), row.names = FALSE)

# ---- IRR et AME ----
irr_table <- function(m) {
  ct <- coef_robuste(m)
  data.frame(variable = rownames(ct), IRR = exp(ct[,1]),
             IC_bas = exp(ct[,1] - 1.96*ct[,2]), IC_haut = exp(ct[,1] + 1.96*ct[,2]),
             p = ct[,4], row.names = NULL)
}
cat("\n--- IRR (NB2 M3) ---\n"); print(irr_table(m_nb), digits = 3)

ame_num <- function(m, data, var, h = 0.01) {
  d1 <- data; d0 <- data
  d1[[var]] <- d1[[var]] + h/2; d0[[var]] <- d0[[var]] - h/2
  mean((predict(m, newdata = d1, type = "response") -
          predict(m, newdata = d0, type = "response")) / h)
}
cat("\n--- AME z_theil / z_moran (tous modèles) ---\n")
print(rbind(z_theil = sapply(modeles, ame_num, data = d_mod, var = "z_theil"),
            z_moran = sapply(modeles, ame_num, data = d_mod, var = "z_moran")))