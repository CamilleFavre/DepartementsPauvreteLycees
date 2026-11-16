# =====================================================================
# CODE 7 FINAL v2 - EXTENSIONS MARGOLIS (count data) - avec checkpoints
# Changements v2 : 4 indices de ségrégation ; MC sur z_moran (DGP = coefs M4) ;
#   f_cv depuis formules$M4 ; bootstrap aussi sur z_moran ; MASS::glm.nb.
# PATCHS VILLE (regroupement PLM 75056/69123/13055) :
#   - map_plm() applique aux centroïdes en §0  -> "Centroïdes manquants : 0"
#   - Vuong : terme de troncature -log(1-f(0)) ajouté pour y>0
#   - MC : Xv = design complet de M4 hors FE ; clusterEvalQ(stats)
#   - SAR : affichage de rho + SE + p
# En cas de crash, recharger le dernier checkpoint :
#   cp <- readRDS(file.path(DOSSIER_OUT, "code7_2_boot.rds"))
# =====================================================================
library(arrow); library(dplyr); library(sandwich); library(lmtest)

library(pscl); library(AER); library(MASS); library(sf); library(ggplot2)
library(spdep); library(parallel)
if (requireNamespace("spatialreg", quietly = TRUE)) library(spatialreg)
select <- dplyr::select          # MASS masque dplyr::select
set.seed(2026)

os_name <- Sys.info()[["sysname"]]
if (os_name == "Windows") {
  DOSSIER     <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
} else if (os_name == "Darwin") {
  DOSSIER     <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
} else stop("OS non pris en charge.")
DOSSIER_ECO <- file.path(DOSSIER_OUT, "data_eco")
CK <- function(nom, obj) saveRDS(obj, file.path(DOSSIER_OUT, paste0("code7_", nom, ".rds")))

# Fusion arrondissements PLM -> ville (identique Code 2 v2 / Code 2c / Code 3b)
map_plm <- function(code) {
  dplyr::case_when(
    grepl("^751(0[1-9]|1[0-9]|20)$", code) ~ "75056",
    grepl("^6938[1-9]$", code)             ~ "69123",
    grepl("^132(0[1-9]|1[0-6])$", code)    ~ "13055",
    TRUE ~ code
  )
}

charg    <- readRDS(file.path(DOSSIER_OUT, "modeles_comptage.rds"))
d_mod    <- charg$d_mod
modeles  <- charg$modeles
formules <- charg$formules
m_p      <- modeles$Poisson
cat("N =", nrow(d_mod), "| modeles :", paste(names(modeles), collapse = ", "), "\n")

# ---------------------------------------------------------------------
# 0. CENTROÏDES COMMUNAUX (niveau VILLE : map_plm sur lcog_geo)
# ---------------------------------------------------------------------
DEMI_CARREAU <- 100
files_eco <- list.files(DOSSIER_ECO, pattern = "eco_dept_.*[.]parquet$", full.names = TRUE)
centroides <- bind_rows(lapply(files_eco, function(f) {
  read_parquet(f, col_select = c("lcog_geo", "x_l93", "y_l93", "ind")) %>%
    mutate(code_insee = map_plm(substr(as.character(lcog_geo), 1, 5))) %>%
    filter(code_insee %in% d_mod$code_insee, is.finite(ind), ind > 0) %>%
    group_by(code_insee) %>%
    summarise(x_pop = weighted.mean(x_l93 + DEMI_CARREAU, ind),
              y_pop = weighted.mean(y_l93 + DEMI_CARREAU, ind), .groups = "drop")
}))
d_mod <- d_mod %>% left_join(centroides, by = "code_insee")
cat("Centroïdes manquants :", sum(is.na(d_mod$x_pop)), "\n")   # attendu : 0
CK("0_data", list(d_mod = d_mod, modeles = modeles, formules = formules))
cat("Checkpoint 0 OK\n")

# ---------------------------------------------------------------------
# 1. TESTS DE SPÉCIFICATION (CT, Van den Broek, Vuong manuel)
# ---------------------------------------------------------------------
mu  <- fitted(m_p)
ct2 <- dispersiontest(m_p, trafo = 2)
p0  <- mean(dpois(0, mu))
vdB <- (sum(d_mod$y == 0) - length(mu) * p0) / sqrt(length(mu) * p0 * (1 - p0))
p_vdB <- 2 * pnorm(-abs(vdB))

m_hp    <- modeles$Hurdle_Pois
ll_pois <- dpois(d_mod$y, mu, log = TRUE)
if (!is.null(m_hp)) {
  pi_i <- pmin(pmax(predict(m_hp, type = "zero"), 1e-10), 1 - 1e-10)   # anti-NaN
  mu_i <- predict(m_hp, type = "count")
  # log-vraisemblance hurdle-Poisson : y=0 -> log(pi) ; y>0 -> log(1-pi) + log f(y) - log(1 - f(0))
  ll_h <- ifelse(d_mod$y == 0, log(pi_i),
                 log(1 - pi_i) + dpois(d_mod$y, mu_i, log = TRUE) - log(1 - dpois(0, mu_i)))
  m_i  <- ll_h - ll_pois
  vuong_stat <- sum(m_i) / (sqrt(length(m_i)) * sd(m_i))
  p_vuong <- 2 * pnorm(-abs(vuong_stat))
} else { vuong_stat <- NA_real_; p_vuong <- NA_real_ }

tab_spec <- data.frame(
  test = c("Cameron-Trivedi alpha (trafo=2)", "Van den Broek (excès de zéros)", "Vuong Hurdle-Pois vs Poisson"),
  stat = c(ct2$estimate, vdB, vuong_stat),
  p    = c(ct2$p.value, p_vdB, p_vuong))
print(tab_spec)
write.csv(tab_spec, file.path(DOSSIER_OUT, "tab_spec.csv"), row.names = FALSE)
CK("1_spec", tab_spec)
cat("Checkpoint 1 OK\n")


# ---------------------------------------------------------------------
# 2. ÉCART beta(M4)-beta(M3) : delta method + bootstrap clusteré B=500
#    (sur z_exposition_riches pour tab_boot ; z_moran aussi sauvegardé)
# ---------------------------------------------------------------------
m3  <- glm(formules$M3, family = poisson, data = d_mod)
b3  <- coef(m3)["z_exposition_riches"]; b4 <- coef(m_p)["z_exposition_riches"]
V3  <- vcovCL(m3, cluster = d_mod$dept, type = "HC1")
V4  <- vcovCL(m_p, cluster = d_mod$region_acad, type = "HC1")
se3 <- sqrt(V3["z_exposition_riches", "z_exposition_riches"])
se4 <- sqrt(V4["z_exposition_riches", "z_exposition_riches"])
se_delta <- sqrt(se3^2 + se4^2)
diff_obs <- b4 - b3

B <- 500
clus4 <- split(seq_len(nrow(d_mod)), d_mod$region_acad)
boot_diff  <- rep(NA_real_, B)
boot_moran <- rep(NA_real_, B)
set.seed(7)
pb <- txtProgressBar(max = B, style = 3)
for (b in seq_len(B)) {
  db  <- d_mod[unlist(clus4[sample(names(clus4), replace = TRUE)]), ]
  m4b <- tryCatch(glm(formules$M4, family = poisson, data = db), error = function(e) NULL)
  m3b <- tryCatch(glm(formules$M3, family = poisson, data = db), error = function(e) NULL)
  if (!is.null(m4b) && !is.null(m3b) &&
      "z_exposition_riches" %in% names(coef(m4b)) && "z_exposition_riches" %in% names(coef(m3b)))
    boot_diff[b] <- coef(m4b)["z_exposition_riches"] - coef(m3b)["z_exposition_riches"]
  if (!is.null(m4b) && !is.null(m3b) &&
      "z_moran" %in% names(coef(m4b)) && "z_moran" %in% names(coef(m3b)))
    boot_moran[b] <- coef(m4b)["z_moran"] - coef(m3b)["z_moran"]
  setTxtProgressBar(pb, b)
}
close(pb)
boot_diff  <- boot_diff[!is.na(boot_diff)]
boot_moran <- boot_moran[!is.na(boot_moran)]
ic_boot <- quantile(boot_diff, c(.025, .975))
tab_boot <- data.frame(
  estimateur = c("delta method", "bootstrap cluster B=500"),
  diff = diff_obs, se = c(se_delta, sd(boot_diff)),
  ic_bas = c(diff_obs - 1.96 * se_delta, unname(ic_boot[1])),
  ic_haut = c(diff_obs + 1.96 * se_delta, unname(ic_boot[2])),
  p = c(2 * pnorm(-abs(diff_obs / se_delta)), 2 * min(mean(boot_diff > 0), mean(boot_diff < 0))))
print(tab_boot)
write.csv(tab_boot, file.path(DOSSIER_OUT, "tab_boot.csv"), row.names = FALSE)
p <- ggplot(data.frame(d = boot_diff), aes(d)) +
  geom_histogram(bins = 40, fill = "#6a51a3", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  geom_vline(xintercept = unname(ic_boot), color = "darkred", linetype = "dotted") +
  theme_minimal(base_size = 11) +
  labs(x = expression(beta[M4] - beta[M3]), title = "Bootstrap clusteré : écart de coef (z_exposition_riches)")
ggsave(file.path(DOSSIER_OUT, "fig_boot_diff.png"), p, width = 8, height = 5, dpi = 150)
CK("2_boot", list(tab_boot = tab_boot, boot_diff = boot_diff, boot_moran = boot_moran))
cat("Checkpoint 2 OK\n")

# ---------------------------------------------------------------------
# 3. MONTE CARLO M=200, parallèle - DGP = coefficients arrondis de M4 (hors FE)
# ---------------------------------------------------------------------
PARAM <- "z_moran"    # paramètre suivi (significatif en M4)
# design du MC = toutes les variables de M4 hors FE région académique
Xv <- model.matrix(reformulate(setdiff(all.vars(formules$M4), c("y", "region_acad")),
                               response = "y"), data = d_mod)
beta_vrai <- setNames(round(coef(m_p)[colnames(Xv)], 2), colnames(Xv))
beta_star <- unname(beta_vrai[PARAM])
cat("Paramètre suivi :", PARAM, "| beta* =", beta_star, "\n")   # attendu : z_moran | 0.13
n <- nrow(d_mod)
sim_once <- function(dgp_nb2, size_nb = 1.5) {
  mu <- exp(as.vector(Xv %*% beta_vrai))
  y  <- if (dgp_nb2) rnbinom(n, mu = mu, size = size_nb) else rpois(n, mu)
  df <- data.frame(y = y, Xv[, -1, drop = FALSE])
  mp <- tryCatch(glm(y ~ ., family = poisson, data = df), error = function(e) NULL)
  if (is.null(mp)) return(c(b_pois = NA_real_, se_pois = NA_real_, b_nb = NA_real_, se_nb = NA_real_))
  out <- c(b_pois  = unname(coef(mp)[PARAM]),
           se_pois = sqrt(diag(vcov(mp)))[PARAM],
           b_nb = NA_real_, se_nb = NA_real_)
  if (dgp_nb2) {
    mn <- tryCatch(MASS::glm.nb(y ~ ., data = df, control = glm.control(maxit = 100, epsilon = 1e-6)),
                   error = function(e) NULL)   # MASS:: car les workers n'ont pas MASS chargé
    if (!is.null(mn) && PARAM %in% names(coef(mn))) {
      out["b_nb"]  <- unname(coef(mn)[PARAM])
      out["se_nb"] <- sqrt(diag(vcov(mn)))[PARAM]
    }
  }
  out
}
M <- 200
ncores <- max(1, detectCores() - 1)
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(stats))                 # glm/vcov dispo dans les workers
parallel::clusterSetRNGStream(cl, 2026)          # reproductibilité du MC
clusterExport(cl, c("Xv", "beta_vrai", "n", "sim_once", "PARAM"))
mc_pois_dgp <- parSapply(cl, seq_len(M), function(i) sim_once(FALSE))
mc_nb_dgp   <- parSapply(cl, seq_len(M), function(i) sim_once(TRUE))
stopCluster(cl)
rownames(mc_pois_dgp) <- c("b_pois", "se_pois", "b_nb", "se_nb")
rownames(mc_nb_dgp)   <- c("b_pois", "se_pois", "b_nb", "se_nb")
g <- function(x) mean(x, na.rm = TRUE)
s <- function(x) sd(x, na.rm = TRUE)
tab_mc <- data.frame(
  dgp = rep(c("Poisson", "NB2 (size=1.5)"), each = 2),
  estimateur = rep(c("Poisson", "NB2"), 2),
  biais   = c(g(mc_pois_dgp["b_pois", ]) - beta_star, g(mc_pois_dgp["b_nb", ]) - beta_star,
              g(mc_nb_dgp["b_pois", ])   - beta_star, g(mc_nb_dgp["b_nb", ])   - beta_star),
  sd_emp  = c(s(mc_pois_dgp["b_pois", ]), s(mc_pois_dgp["b_nb", ]),
              s(mc_nb_dgp["b_pois", ]),   s(mc_nb_dgp["b_nb", ])),
  mean_se = c(g(mc_pois_dgp["se_pois", ]), g(mc_pois_dgp["se_nb", ]),
              g(mc_nb_dgp["se_pois", ]),   g(mc_nb_dgp["se_nb", ])),
  couv95  = c(mean(abs(mc_pois_dgp["b_pois", ] - beta_star) < 1.96 * mc_pois_dgp["se_pois", ], na.rm = TRUE),
              mean(abs(mc_pois_dgp["b_nb", ]   - beta_star) < 1.96 * mc_pois_dgp["se_nb", ],   na.rm = TRUE),
              mean(abs(mc_nb_dgp["b_pois", ]   - beta_star) < 1.96 * mc_nb_dgp["se_pois", ],   na.rm = TRUE),
              mean(abs(mc_nb_dgp["b_nb", ]     - beta_star) < 1.96 * mc_nb_dgp["se_nb", ],     na.rm = TRUE))
) %>% mutate(across(where(is.numeric), ~ round(., 4)))
print(tab_mc)
write.csv(tab_mc, file.path(DOSSIER_OUT, "tab_mc.csv"), row.names = FALSE)
mc_df <- data.frame(b = c(mc_pois_dgp["b_pois", ], mc_pois_dgp["b_nb", ],
                          mc_nb_dgp["b_pois", ], mc_nb_dgp["b_nb", ]),
                    dgp = rep(c("DGP Poisson", "DGP Poisson", "DGP NB2", "DGP NB2"), each = M),
                    est = rep(c("Poisson", "NB2"), 2, each = M))
p_mc <- ggplot(mc_df, aes(b, fill = est)) +
  geom_histogram(alpha = 0.5, bins = 40, position = "identity") +
  geom_vline(xintercept = beta_star, color = "red", linetype = "dashed") +
  facet_wrap(~ dgp) + theme_minimal(base_size = 10) +
  labs(x = paste0("hat(beta) de ", PARAM), title = paste("Monte Carlo M=200 : vraie valeur =", beta_star))
ggsave(file.path(DOSSIER_OUT, "fig_mc_count.png"), p_mc, width = 10, height = 5, dpi = 150)
CK("3_mc", list(tab_mc = tab_mc, pois = mc_pois_dgp, nb = mc_nb_dgp, beta_star = beta_star, PARAM = PARAM))
cat("Checkpoint 3 OK\n")

# ---------------------------------------------------------------------
# 4. SPATIAL COMMUNAL : Moran + SAR sous 2 poids (8-NN, 1/d 200km) ;
#    Moran des résidus Poisson ET Hurdle ; CV blocs vs aléatoire sous Poisson ET Hurdle ;
#    baselines naïves + MAE dans le même protocole
# ---------------------------------------------------------------------
# library(pscl) en haut du script (hurdle)

d_sp   <- d_mod[complete.cases(d_mod[, c("x_pop", "y_pop")]), ]
coords <- as.matrix(d_sp[, c("x_pop", "y_pop")])          # L93, mètres

# --- poids principal : 8-NN, row-standardisé ---
nb8 <- knn2nb(knearneigh(coords, k = 8))
lw  <- nb2listw(nb8, style = "W")

# --- poids robustesse : 1/d_ij tronqué à 200 km, row-standardisé ---
CUTOFF_IDW <- 200
if (all(c("lon", "lat") %in% names(d_sp))) {
  D_km <- tryCatch(geosphere::distm(cbind(d_sp$lon, d_sp$lat), fun = geosphere::distHaversine) / 1000,
                   error = function(e) as.matrix(dist(coords)) / 1000)
} else {
  D_km <- as.matrix(dist(coords)) / 1000
}
W_idw <- 1 / (D_km + 1e-6); W_idw[D_km > CUTOFF_IDW] <- 0; diag(W_idw) <- 0
dimnames(W_idw) <- list(d_sp$code_insee, d_sp$code_insee)
lw_idw <- mat2listw(W_idw, style = "W")

# --- 2 modèles "full" (avec FE) : Poisson M4 et Hurdle M4 ---
m_p_sp <- glm(formules$M4, family = poisson, data = d_sp)
m_h_sp <- tryCatch(hurdle(formules$M4, data = d_sp, dist = "poisson", zero.dist = "poisson"),
                   error = function(e) NULL)

# --- pipeline Moran + SAR sous une matrice de poids (Poisson M4) ---
run_moran_sar <- function(lw_x, tag) {
  mt <- lm.morantest(m_p_sp, lw_x)
  cat(sprintf("[%s] Moran I residus M4 (Poisson) : I = %.4f, p = %.4g\n", tag, mt$estimate[1], mt$p.value))
  m_sar <- tryCatch(
    if (requireNamespace("spatialreg", quietly = TRUE))
      spatialreg::lagsarlm(formules$M4, data = d_sp, listw = lw_x, method = "eigen", quiet = TRUE)
    else spdep::lagsarlm(formules$M4, data = d_sp, listw = lw_x, method = "eigen", quiet = TRUE),
    error = function(e) { message("SAR non estimé (", tag, ") : ", e$message); NULL })
  rho_v <- se_rho <- p_rho <- NA_real_
  if (!is.null(m_sar)) {
    cf <- coef(m_sar)
    i  <- grep("^rho$", names(cf)); if (!length(i)) i <- grep("rho", names(cf))
    if (length(i)) rho_v <- unname(cf[i[1]])
    V_sar <- tryCatch(vcov(m_sar), error = function(e) NULL)
    if (!is.null(V_sar)) {
      j <- grep("^rho$", rownames(V_sar)); if (!length(j)) j <- grep("rho", rownames(V_sar))
      if (length(j)) {
        se_rho <- sqrt(V_sar[j[1], j[1]])
        p_rho  <- if (is.finite(se_rho) && se_rho > 0) 2 * pnorm(-abs(rho_v / se_rho)) else NA_real_
      }
    }
  }
  cat(sprintf("[%s] SAR rho = %.4f | SE = %.4f | p = %.4g\n", tag, rho_v, se_rho, p_rho))
  list(moran = mt, m_sar = m_sar, sar_rho = c(rho = rho_v, se = se_rho, p = p_rho))
}

res_knn <- run_moran_sar(lw,     "8-NN")
res_idw <- run_moran_sar(lw_idw, "1/d trq 200km")

mt <- res_knn$moran; m_sar <- res_knn$m_sar; sar_rho <- res_knn$sar_rho   # alias principal = 8-NN
moran_idw <- res_idw$moran; m_sar_idw <- res_idw$m_sar; sar_rho_idw <- res_idw$sar_rho

# --- Moran des résidus du HURDLE sous les 2 poids (moran.test car hurdle n'est pas un lm) ---
if (!is.null(m_h_sp)) {
  e_h <- residuals(m_h_sp, type = "response")             # y - E[y]
  moran_h_knn <- moran.test(e_h, lw)
  moran_h_idw <- moran.test(e_h, lw_idw)
  cat(sprintf("[hurdle 8-NN] I = %.4f, p = %.4g | [hurdle 1/d 200km] I = %.4f, p = %.4g\n",
              moran_h_knn$estimate[1], moran_h_knn$p.value,
              moran_h_idw$estimate[1], moran_h_idw$p.value))
} else { moran_h_knn <- moran_h_idw <- NULL }

# ================= CV : 2 modèles x 2 découpages =================
# /!\ ordre des arguments : (te, fitfun) -> sapply passe l'élément en 1er argument
f_cv <- reformulate(setdiff(all.vars(formules$M4), c("y", "region_acad")), response = "y")  # M4 sans FE
set.seed(1)
bloc <- kmeans(scale(coords), centers = 10)$cluster

fit_pois <- function(train) glm(f_cv, family = poisson, data = train)
fit_hurd <- function(train) tryCatch(hurdle(f_cv, data = train, dist = "poisson", zero.dist = "poisson"),
                                     error = function(e) NULL)

mse_fold <- function(te, fitfun) {
  m <- fitfun(d_sp[!te, ])
  if (is.null(m)) return(NA_real_)                          # fold non convergé (hurdle)
  mean((d_sp$y[te] - predict(m, d_sp[te, ], type = "response"))^2)   # hurdle: response = E[y]
}
mae_fold <- function(te, fitfun) {
  m <- fitfun(d_sp[!te, ])
  if (is.null(m)) return(NA_real_)
  mean(abs(d_sp$y[te] - predict(m, d_sp[te, ], type = "response")))
}
naive_fold <- function(te, mode = c("mean", "zero")) {
  if (mode == "mean") mean((d_sp$y[te] - mean(d_sp$y[!te]))^2)   # moyenne du train
  else               mean(d_sp$y[te]^2)                            # toujours 0
}

blocks_te <- lapply(sort(unique(bloc)), function(b) bloc == b)      # 10 folds spatiaux
rand_te   <- lapply(1:20, function(r) {                              # 20 folds aléatoires 10 %
  set.seed(r)
  te <- logical(nrow(d_sp)); te[sample(nrow(d_sp), round(0.10 * nrow(d_sp)))] <- TRUE
  te
})
paris_rd <- sapply(rand_te, function(te) any(d_sp$y[te] == max(d_sp$y)))

err_sp_pois <- sapply(blocks_te, mse_fold, fit_pois)
err_rd_pois <- sapply(rand_te,   mse_fold, fit_pois)
err_sp_hurd <- sapply(blocks_te, mse_fold, fit_hurd)
err_rd_hurd <- sapply(rand_te,   mse_fold, fit_hurd)

mae_sp_pois <- sapply(blocks_te, mae_fold, fit_pois)
mae_rd_pois <- sapply(rand_te,   mae_fold, fit_pois)
mae_sp_hurd <- sapply(blocks_te, mae_fold, fit_hurd)
mae_rd_hurd <- sapply(rand_te,   mae_fold, fit_hurd)

# baselines naïves dans le MÊME protocole (mêmes folds)
naive_tab <- data.frame(
  decoupage     = c("blocs k=10", "aléatoire 10% x20"),
  naive_moyenne = c(mean(sapply(blocks_te, naive_fold, "mean")), mean(sapply(rand_te, naive_fold, "mean"))),
  naive_zero    = c(mean(sapply(blocks_te, naive_fold, "zero")), mean(sapply(rand_te, naive_fold, "zero"))))
print(naive_tab)

mse_bloc_pois <- mean(err_sp_pois, na.rm = TRUE)
mse_bloc_hurd <- mean(err_sp_hurd, na.rm = TRUE)
tab_cv <- data.frame(
  modele    = rep(c("Poisson (f_cv)", "Hurdle-Pois (f_cv)"), each = 2),
  decoupage = rep(c("CV blocs spatiaux (k-means k=10)", "CV aléatoire (10%, 20 rép.)"), 2),
  mse       = c(mse_bloc_pois, mean(err_rd_pois, na.rm = TRUE),
                mse_bloc_hurd, mean(err_rd_hurd, na.rm = TRUE)),
  mae       = c(mean(mae_sp_pois, na.rm = TRUE), mean(mae_rd_pois, na.rm = TRUE),
                mean(mae_sp_hurd, na.rm = TRUE), mean(mae_rd_hurd, na.rm = TRUE)),
  n_folds_NA = c(sum(is.na(err_sp_pois)), sum(is.na(err_rd_pois)),
                 sum(is.na(err_sp_hurd)), sum(is.na(err_rd_hurd))))
tab_cv$ratio_vs_blocs   <- tab_cv$mse / rep(c(mse_bloc_pois, mse_bloc_hurd), each = 2)
tab_cv$mse_naive_moy    <- rep(naive_tab$naive_moyenne, 2)     # lignes 1,2 = blocs ; 3,4 = aléa
tab_cv$skill_vs_naive   <- 1 - tab_cv$mse / tab_cv$mse_naive_moy
# --- CV pondérée par le nombre de villes du fold (MSE poolée par ville) ---
# taille_bloc alignée avec blocks_te (levels triés 1..10)
taille_bloc <- as.integer(table(factor(bloc, levels = sort(unique(bloc)))))
wtd <- function(err, w) weighted.mean(err, w, na.rm = TRUE)
tab_cv$mse_ponderee <- c(wtd(err_sp_pois, taille_bloc), mean(err_rd_pois, na.rm = TRUE),
                         wtd(err_sp_hurd, taille_bloc), mean(err_rd_hurd, na.rm = TRUE))
tab_cv$mae_ponderee <- c(wtd(mae_sp_pois, taille_bloc), mean(mae_rd_pois, na.rm = TRUE),
                         wtd(mae_sp_hurd, taille_bloc), mean(mae_rd_hurd, na.rm = TRUE))
# NB : en CV aléatoire tous les folds ont la même taille (round(0.10*n) = 92 villes),
#     donc la pondération y est sans effet (mse_ponderee = mse) ; c'est en blocs qu'elle compte.

# baseline naïve pondérée (blocs) pour le skill pondéré
naive_sp_moy  <- sapply(blocks_te, naive_fold, "mean")
naive_sp_zero <- sapply(blocks_te, naive_fold, "zero")
naive_tab$naive_moyenne_pond <- c(wtd(naive_sp_moy, taille_bloc), NA_real_)
naive_tab$naive_zero_pond    <- c(wtd(naive_sp_zero, taille_bloc), NA_real_)
tab_cv$skill_vs_naive_pond <- 1 - tab_cv$mse_ponderee /
  rep(c(naive_tab$naive_moyenne_pond[1], naive_tab$naive_moyenne[2]), 2)

print(tab_cv)
write.csv(tab_cv, file.path(DOSSIER_OUT, "tab_cv.csv"), row.names = FALSE)
write.csv(naive_tab, file.path(DOSSIER_OUT, "tab_cv_naive.csv"), row.names = FALSE)

# --- diagnostic : écart hurdle - Poisson par bloc (Paris ?) ---
cmp <- data.frame(bloc = sort(unique(bloc)), n = as.integer(table(bloc)),
                  mse_pois = err_sp_pois, mse_hurd = err_sp_hurd)
cmp$diff <- cmp$mse_hurd - cmp$mse_pois
print(cmp[order(-cmp$diff), ])
i_paris <- which.max(d_sp$y); cat("bloc de Paris :", bloc[i_paris], "\n")

# décomposition dans le fold du bloc de Paris : zéros vs positifs, ligne Paris
b <- bloc[i_paris]; te <- bloc == b
m_pois <- glm(f_cv, family = poisson, data = d_sp[!te, ])
m_hurd <- hurdle(f_cv, data = d_sp[!te, ], dist = "poisson", zero.dist = "poisson")
tst <- data.frame(y = d_sp$y[te],
                  yhat_pois = predict(m_pois, d_sp[te, ], type = "response"),
                  yhat_hurd = predict(m_hurd, d_sp[te, ], type = "response"))
tst$sq_pois <- (tst$y - tst$yhat_pois)^2; tst$sq_hurd <- (tst$y - tst$yhat_hurd)^2
cat("MSE zéros (y=0)   pois/hurd :", mean(tst$sq_pois[tst$y == 0]), "/", mean(tst$sq_hurd[tst$y == 0]), "\n")
cat("MSE positifs (y>0) pois/hurd :", mean(tst$sq_pois[tst$y > 0]), "/", mean(tst$sq_hurd[tst$y > 0]), "\n")
print(tst[order(-tst$sq_hurd), ][1:8, ])
cat("MAE blocs  Poisson :", mean(mae_sp_pois, na.rm = TRUE), " | Hurdle :", mean(mae_sp_hurd, na.rm = TRUE), "\n")
cat("MAE aléa   Poisson :", mean(mae_rd_pois, na.rm = TRUE), " | Hurdle :", mean(mae_rd_hurd, na.rm = TRUE), "\n")

# --- figure : MSE par modèle x découpage + baseline naïve ---
plot_df <- rbind(
  tab_cv[, c("modele", "decoupage", "mse")],
  data.frame(modele = "Naïf (moy. train)",
             decoupage = c("CV blocs spatiaux (k-means k=10)", "CV aléatoire (10%, 20 rép.)"),
             mse = naive_tab$naive_moyenne))
p_cv <- ggplot(plot_df, aes(decoupage, mse, fill = decoupage)) +
  geom_col(width = 0.6) + geom_text(aes(label = sprintf("%.3f", mse)), vjust = -0.4) +
  facet_wrap(~ modele) + theme_minimal() + theme(legend.position = "none") +
  labs(x = NULL, y = "MSE",
       title = "CV par blocs spatiaux vs aléatoire, sous Poisson, Hurdle et modèle naïf",
       subtitle = "La CV aléatoire est dominée par Paris (y = 103) dans le test + fuite spatiale ; la CV par blocs est la référence")
ggsave(file.path(DOSSIER_OUT, "fig_cv_communes.png"), p_cv, width = 11, height = 5, dpi = 150)

CK("4_spatial", list(tab_cv = tab_cv, naive_tab = naive_tab, cmp = cmp,
                     bloc = bloc, paris_rd = paris_rd, i_paris = i_paris, tst = tst,
                     err_sp_pois = err_sp_pois, err_rd_pois = err_rd_pois,
                     err_sp_hurd = err_sp_hurd, err_rd_hurd = err_rd_hurd,
                     mae_sp_pois = mae_sp_pois, mae_rd_pois = mae_rd_pois,
                     mae_sp_hurd = mae_sp_hurd, mae_rd_hurd = mae_rd_hurd,
                     moran = mt, m_sar = m_sar, sar_rho = sar_rho,
                     moran_idw = moran_idw, m_sar_idw = m_sar_idw, sar_rho_idw = sar_rho_idw, lw_idw = lw_idw,
                     moran_h_knn = moran_h_knn, moran_h_idw = moran_h_idw, m_h_sp = m_h_sp))
cat("Checkpoint 4 OK\n")

# =====================================================================
# 4b. SENSIBILITÉ DE LA CV SPATIALE : portée, nb de blocs, forme, buffer
# =====================================================================
f3 <- function(x) sprintf("%.3f", x)
D_km    <- as.matrix(dist(coords)) / 1000          # km (L93 en mètres)
i_paris <- which.max(d_sp$y)

# ---------- (A) Portée de l'autocorrélation résiduelle ----------
m_cv   <- glm(f_cv, family = poisson, data = d_sp)  # modèle validé (sans FE)
e_pear <- residuals(m_cv, type = "pearson")

moran_band <- function(e, Dk, d1, d2, nperm = 199) {
  W  <- (Dk > d1 & Dk <= d2) * 1                    # poids binaires de la bande
  S0 <- sum(W)
  if (S0 == 0) return(c(I = NA_real_, p = NA_real_, npairs = 0))
  z  <- e - mean(e); n <- length(z)
  stat <- function(zz) (n / S0) * sum(zz * (W %*% zz)) / sum(zz^2)
  I0 <- stat(z)
  Ip <- replicate(nperm, stat(sample(z)))
  c(I = I0, p = (1 + sum(Ip >= I0)) / (nperm + 1), npairs = S0 / 2)
}

set.seed(11)
breaks <- c(0, 10, 25, 50, 75, 100, 150, 200, 300, 400)
corr <- do.call(rbind, lapply(seq_len(length(breaks) - 1), function(j) {
  a <- moran_band(e_pear, D_km, breaks[j], breaks[j + 1])
  b <- moran_band(e_pear[-i_paris], D_km[-i_paris, -i_paris], breaks[j], breaks[j + 1])
  data.frame(d1 = breaks[j], d2 = breaks[j + 1], npairs = unname(a["npairs"]),
             I = unname(a["I"]), p = unname(a["p"]),
             I_sP = unname(b["I"]), p_sP = unname(b["p"]))
}))
print(corr)

ns <- which(corr$p > 0.05)                          # 1re bande non significative
RANGE_KM <- if (length(ns)) corr$d1[min(ns)] else max(breaks)
RANGE_KM <- max(RANGE_KM, 25)                       # <- À VÉRIFIER / ajuster à la main sur la figure
cat("Portée estimée de l'autocorrélation résiduelle :", RANGE_KM, "km\n")

cl <- rbind(data.frame(mid = (corr$d1 + corr$d2) / 2, I = corr$I,    p = corr$p,    serie = "Toutes villes"),
            data.frame(mid = (corr$d1 + corr$d2) / 2, I = corr$I_sP, p = corr$p_sP, serie = "Sans Paris"))
p_corr <- ggplot(cl, aes(mid, I, colour = serie)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_vline(xintercept = RANGE_KM, linetype = "dashed", colour = "red") +
  geom_line() + geom_point(aes(shape = p < 0.05), size = 2.5) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "p < 0,05") +
  scale_x_continuous(breaks = breaks) +
  labs(x = "Distance entre villes (km, milieu de bande)", y = "I de Moran (résidus de Pearson)",
       colour = NULL, title = "Corrélogramme des résidus Poisson (M4 sans FE)",
       subtitle = "Ligne rouge = portée retenue pour la taille des blocs / le buffer") +
  theme_minimal(base_size = 11)
ggsave(file.path(DOSSIER_OUT, "fig_correlogramme_residus.png"), p_corr, width = 8, height = 5, dpi = 150)
print(p_corr)

# ---------- (B) Moteur de CV : plis quelconques + buffer optionnel ----------
safe_pred <- function(m, nd) {
  if (is.null(m)) return(rep(NA_real_, nrow(nd)))
  tryCatch(as.numeric(predict(m, nd, type = "response")), error = function(e) rep(NA_real_, nrow(nd)))
}
cv_oof <- function(fold, buffer_km = 0) {
  n   <- nrow(d_sp)
  oof <- data.frame(pois = rep(NA_real_, n), hurd = rep(NA_real_, n), naive = rep(NA_real_, n))
  n_train <- integer(0)
  for (f in sort(unique(fold))) {
    te <- fold == f
    tr <- !te
    if (buffer_km > 0) {                            # retire du train les villes proches du test
      dmin <- apply(D_km[, te, drop = FALSE], 1, min)
      tr   <- tr & dmin >= buffer_km
    }
    n_train <- c(n_train, sum(tr))
    trd <- d_sp[tr, ]
    mp  <- tryCatch(suppressWarnings(glm(f_cv, family = poisson, data = trd)), error = function(e) NULL)
    mh  <- tryCatch(suppressWarnings(hurdle(f_cv, data = trd, dist = "poisson", zero.dist = "poisson")),
                    error = function(e) NULL)
    oof$pois[te]  <- safe_pred(mp, d_sp[te, ])
    oof$hurd[te]  <- safe_pred(mh, d_sp[te, ])
    oof$naive[te] <- mean(trd$y)
  }
  attr(oof, "n_train") <- n_train
  oof
}

y_obs <- d_sp$y; non_paris <- seq_along(y_obs) != i_paris
cv_summary <- function(p, fold) {
  se <- (y_obs - p)^2
  c(mse = mean(se, na.rm = TRUE),                                  # poolée (1 ville = 1 poids)
    mse_fold = mean(tapply(se, fold, mean, na.rm = TRUE), na.rm = TRUE),
    mse_sP = mean(se[non_paris], na.rm = TRUE),                    # hors Paris
    mae = mean(abs(y_obs - p), na.rm = TRUE),
    n_na = sum(is.na(p)))
}

# ---------- (C) Définition des plans de CV ----------
mk <- function(nom, type, fold, buffer = 0)
  list(nom = nom, type = type, fold = as.integer(factor(fold)), buffer = buffer)
km_fold   <- function(k) { set.seed(1); kmeans(scale(coords), centers = k)$cluster }  # même appel que §4
grid_fold <- function(side_km) paste(floor(coords[, 1] / (side_km * 1000)),
                                     floor(coords[, 2] / (side_km * 1000)))
BUFFERS <- sort(unique(c(50, 100, round(RANGE_KM))))

designs <- c(
  lapply(c(5, 10, 20, 30, 50), function(k) mk(sprintf("kmeans k=%d", k), "kmeans", km_fold(k))),
  lapply(c(100, 200, 300),     function(s) mk(sprintf("grille %d km", s), "grille", grid_fold(s))),
  list(mk("région académique", "admin", d_sp$region_acad)),
  unlist(lapply(c(10, 30), function(k) lapply(BUFFERS, function(r)
    mk(sprintf("kmeans k=%d + buffer %d km", k, r), "kmeans+buffer", km_fold(k), r))),
    recursive = FALSE)
)
cat("k=10 identique à `bloc` de §4 :", all(designs[[2]]$fold == as.integer(factor(bloc))), "\n")

# ---------- (D) Exécution (quelques minutes : ~300 plis x 2 modèles) ----------
res <- do.call(rbind, lapply(designs, function(dz) {
  t0  <- Sys.time()
  oof <- cv_oof(dz$fold, dz$buffer)
  cat(sprintf("%-34s %5.0f s\n", dz$nom, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  sp <- cv_summary(oof$pois, dz$fold); sh <- cv_summary(oof$hurd, dz$fold); sn <- cv_summary(oof$naive, dz$fold)
  data.frame(design = dz$nom, type = dz$type,
             n_plis = length(unique(dz$fold)), buffer_km = dz$buffer,
             taille_test = round(mean(table(dz$fold)), 1),
             n_train_moy = round(mean(attr(oof, "n_train"))),
             mse_pois = sp[["mse"]], mse_hurd = sh[["mse"]], mse_naif = sn[["mse"]],
             mse_pois_sansParis = sp[["mse_sP"]], mse_hurd_sansParis = sh[["mse_sP"]],
             mse_naif_sansParis = sn[["mse_sP"]],
             mse_pli_pois = sp[["mse_fold"]], mse_pli_hurd = sh[["mse_fold"]],
             mae_pois = sp[["mae"]], mae_hurd = sh[["mae"]],
             n_NA_pois = sp[["n_na"]], n_NA_hurd = sh[["n_na"]])
}))
res$skill_pois <- 1 - res$mse_pois / res$mse_naif
res$skill_hurd <- 1 - res$mse_hurd / res$mse_naif
res$gagnant            <- ifelse(res$mse_hurd < res$mse_pois, "Hurdle", "Poisson")
res$gagnant_sansParis  <- ifelse(res$mse_hurd_sansParis < res$mse_pois_sansParis, "Hurdle", "Poisson")

cat("Check k=10 : MSE poolée Poisson =", f3(res$mse_pois[res$design == "kmeans k=10"]),
    "| tab_cv$mse_ponderee[1] =", f3(tab_cv$mse_ponderee[1]), "\n")
print(res %>% mutate(across(where(is.numeric), ~ round(., 3))))
write.csv(res,  file.path(DOSSIER_OUT, "tab_cv_sensibilite.csv"), row.names = FALSE)
write.csv(corr, file.path(DOSSIER_OUT, "tab_correlogramme.csv"),  row.names = FALSE)

# ---------- (E) Figure ----------
lg <- do.call(rbind, lapply(c("pois", "hurd", "naif"), function(m) rbind(
  data.frame(design = res$design, modele = m, panel = "Avec Paris", mse = res[[paste0("mse_", m)]]),
  data.frame(design = res$design, modele = m, panel = "Sans Paris", mse = res[[paste0("mse_", m, "_sansParis")]]))))
lg$design <- factor(lg$design, levels = res$design)
lg$modele <- factor(lg$modele, levels = c("pois", "hurd", "naif"),
                    labels = c("Poisson", "Hurdle", "Naïf (moy. train)"))
p_sens <- ggplot(lg, aes(design, mse, colour = modele, group = modele)) +
  geom_line(alpha = 0.4) + geom_point(size = 2.2) +
  facet_wrap(~ panel, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "MSE poolée (1 ville = 1 poids)", colour = NULL,
       title = "Sensibilité de la CV spatiale au découpage",
       subtitle = "Buffer = villes du train à moins de r km du pli test retirées")
ggsave(file.path(DOSSIER_OUT, "fig_cv_sensibilite.png"), p_sens, width = 11, height = 8, dpi = 150)
print(p_sens)

CK("4b_cv_sensibilite", list(res = res, corr = corr, RANGE_KM = RANGE_KM))
cat("Checkpoint 4b OK\n")



# ============ Van den Broek : versions cluster-robustes (spatial) ============
# (m_p = Poisson M4 ; adapte le nom si ton Code 7 l'appelle autrement)
m_p <- modeles$Poisson                              # M4 avec FE (cohérent avec les tests CT)
p0  <- dpois(0, fitted(m_p))
n0  <- sum(d_mod$y == 0)
z_vdb_iid <- (n0 - sum(p0)) / sqrt(sum(p0 * (1 - p0)))

u_score <- ifelse(d_mod$y == 0, 1/p0 - 1, -1)
vdB_z <- function(cl, n = length(u_score)) {
  if (is.null(cl)) return(list(z = NA_real_, p = NA_real_, G = 0))
  if (inherits(cl, "kmeans")) cl <- cl$cluster
  cl <- as.character(cl)
  if (length(cl) != n) return(list(z = NA_real_, p = NA_real_, G = length(cl)))
  Ug <- tapply(u_score, cl, sum); G <- length(Ug)
  if (G < 3) return(list(z = NA_real_, p = NA_real_, G = G))
  z <- sum(u_score) / sqrt(G/(G-1) * sum(Ug^2))
  list(z = z, p = 2 * pnorm(-abs(z)), G = G)
}

cl_acad  <- as.character(d_mod$region_acad)
cl_dept  <- as.character(d_mod$dept)
cl_block <- as.character(bloc)                      # ← le bon vecteur (length 925)

vdB_rob <- list(
  iid   = list(z = z_vdb_iid, p = 2*pnorm(-abs(z_vdb_iid)), G = nrow(d_mod)),
  acad  = vdB_z(cl_acad),
  dept  = vdB_z(cl_dept),
  block = vdB_z(cl_block)
)

# wild cluster bootstrap (Rademacher) sur les 3 partitions
wild_boot <- function(cl, B = 2000) {
  Ug <- as.numeric(tapply(d_mod$y == 0, cl, sum) - tapply(p0, cl, sum))
  G <- length(Ug); T_obs <- sum(Ug); V_cl <- G/(G-1) * sum(Ug^2)
  set.seed(2026)
  W  <- matrix(sample(c(-1, 1), G * B, replace = TRUE), nrow = G)
  Tb <- as.numeric(t(W) %*% Ug)
  list(z = T_obs/sqrt(V_cl), p_normal = 2*pnorm(-abs(T_obs/sqrt(V_cl))),
       p_wild = mean(abs(Tb/sqrt(V_cl)) >= abs(T_obs/sqrt(V_cl))), B = B, G = G)   # ← list(), pas c()
}

vdB_rob$wild_acad  <- wild_boot(cl_acad)
vdB_rob$wild_dept  <- wild_boot(cl_dept)
vdB_rob$wild_block <- wild_boot(cl_block)

print(vdB_rob)
sapply(vdB_rob[c("iid","acad","dept","block")], function(x) round(c(x$z, x$p, x$G), 4))

# ============ Vuong : versions cluster-robustes (spatial) ============
# même convention que §1 : m_i = ll_Hurdle - ll_Poisson
h   <- modeles$Hurdle_Pois
m_p <- modeles$Poisson          # déjà défini plus haut, sécurité :
if (!exists("m_p")) m_p <- modeles$Poisson
n   <- nrow(d_mod)

# log-vraisemblances par observation (termes clampés ; troncature -log(1-f_c(0)) incluse)
mu_p <- fitted(m_p)
mu_h <- pmax(predict(h, type = "count"), 1e-10)
pi_h <- pmin(pmax(predict(h, type = "zero"), 1e-10), 1 - 1e-10)
ll_p <- dpois(d_mod$y, mu_p, log = TRUE)
ll_h <- ifelse(d_mod$y == 0, log(pi_h),
               log(1 - pi_h) + dpois(d_mod$y, mu_h, log = TRUE) - log(1 - dpois(0, mu_h)))
m_i   <- ll_h - ll_p
V_iid <- sum(m_i) / (sqrt(n) * sd(m_i))
p_iid <- 2 * pnorm(-abs(V_iid))

vuong_z <- function(cl, n = length(m_i)) {
  if (is.null(cl)) return(list(V = NA_real_, p = NA_real_, G = 0))
  if (inherits(cl, "kmeans")) cl <- cl$cluster
  cl <- as.character(cl)
  if (length(cl) != n) return(list(V = NA_real_, p = NA_real_, G = length(cl)))
  mg <- tapply(m_i, cl, sum); G <- length(mg)
  if (G < 3) return(list(V = NA_real_, p = NA_real_, G = G))
  V <- sum(m_i) / sqrt(G/(G-1) * sum(mg^2))          # variance cluster de sum(m_i)
  list(V = V, p = 2 * pnorm(-abs(V)), G = G)
}

vuong_rob <- list(
  iid   = list(V = V_iid, p = p_iid, G = n),
  acad  = vuong_z(cl_acad),
  dept  = vuong_z(cl_dept),
  block = vuong_z(cl_block)
)

# wild cluster bootstrap (Rademacher, B = 2000) sur les 3 partitions
wild_boot_v <- function(cl, B = 2000) {
  mg <- as.numeric(tapply(m_i, cl, sum))
  G <- length(mg); T_obs <- sum(mg); V_cl <- G/(G-1) * sum(mg^2)
  set.seed(2026)
  W  <- matrix(sample(c(-1, 1), G * B, replace = TRUE), nrow = G)
  Vb <- as.numeric(t(W) %*% mg) / sqrt(V_cl)
  list(V = T_obs/sqrt(V_cl), p_normal = 2*pnorm(-abs(T_obs/sqrt(V_cl))),
       p_wild = mean(abs(Vb) >= abs(T_obs/sqrt(V_cl))), B = B, G = G)   # list(), pas c()
}
vuong_rob$wild_acad  <- wild_boot_v(cl_acad)
vuong_rob$wild_dept  <- wild_boot_v(cl_dept)
vuong_rob$wild_block <- wild_boot_v(cl_block)

print(vuong_rob)


# ---------------------------------------------------------------------
# 5. SAUVEGARDE FINALE (étendue : tous les champs spatiaux, y compris IDW et hurdle)
# ---------------------------------------------------------------------
saveRDS(list(tab_spec = tab_spec, tab_boot = tab_boot, tab_mc = tab_mc, tab_cv = tab_cv,
             naive_tab = naive_tab,                       # <-- AJOUT : baselines naïves (mêmes folds CV)
             vuong_rob = vuong_rob, vdB_rob = vdB_rob,
             moran = mt, m_sar = m_sar, sar_rho = sar_rho,
             moran_idw = moran_idw, m_sar_idw = m_sar_idw, sar_rho_idw = sar_rho_idw, lw_idw = lw_idw,
             moran_h_knn = moran_h_knn, moran_h_idw = moran_h_idw, m_h_sp = m_h_sp,   # <-- AJOUT
             bloc = bloc, paris_rd = paris_rd,
             cmp = cmp, tst = tst, i_paris = i_paris,      # <-- AJOUT : diagnostic bloc de Paris
             err_sp_pois = err_sp_pois, err_rd_pois = err_rd_pois,
             err_sp_hurd = err_sp_hurd, err_rd_hurd = err_rd_hurd,
             mae_sp_pois = mae_sp_pois, mae_rd_pois = mae_rd_pois,   # <-- AJOUT : MAE
             mae_sp_hurd = mae_sp_hurd, mae_rd_hurd = mae_rd_hurd,   # <-- AJOUT : MAE
             boot_diff = boot_diff, boot_moran = boot_moran,
             mc_raw = list(pois = mc_pois_dgp, nb = mc_nb_dgp),
             PARAM = PARAM, beta_star = beta_star),
        file.path(DOSSIER_OUT, "code7.rds"))
cat("Code 7 v2 terminé\n")