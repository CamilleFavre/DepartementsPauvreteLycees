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
# >>> saveRDS(code7.rds) en ajoutant vdB_rob = vdB_rob
# >>> AJOUTER dans le saveRDS de code7.rds :  vdB_rob = vdB_rob

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
# 4. SPATIAL COMMUNAL : Moran I des résidus, SAR (tryCatch), CV par blocs
# ---------------------------------------------------------------------
d_sp   <- d_mod[complete.cases(d_mod[, c("x_pop", "y_pop")]), ]
coords <- as.matrix(d_sp[, c("x_pop", "y_pop")])
nb8 <- knn2nb(knearneigh(coords, k = 8))
lw  <- nb2listw(nb8, style = "W")
m_p_sp <- glm(formules$M4, family = poisson, data = d_sp)
mt <- lm.morantest(m_p_sp, lw)
cat(sprintf("Moran I residus M4 : I = %.4f, p = %.4g\n", mt$estimate[1], mt$p.value))

m_sar <- tryCatch(
  if (requireNamespace("spatialreg", quietly = TRUE))
    spatialreg::lagsarlm(formules$M4, data = d_sp, listw = lw, method = "eigen", quiet = TRUE)
  else spdep::lagsarlm(formules$M4, data = d_sp, listw = lw, method = "eigen", quiet = TRUE),
  error = function(e) { message("SAR non estimé : ", e$message); NULL })
sar_rho <- NULL
if (!is.null(m_sar) && "rho" %in% names(coef(m_sar))) {
  rho_v  <- unname(coef(m_sar)["rho"])
  V_sar  <- tryCatch(vcov(m_sar), error = function(e) NULL)
  se_rho <- if (!is.null(V_sar) && "rho" %in% rownames(V_sar)) sqrt(V_sar["rho", "rho"]) else NA_real_
  p_rho  <- if (!is.na(se_rho) && se_rho > 0) 2 * pnorm(-abs(rho_v / se_rho)) else NA_real_
  sar_rho <- c(rho = rho_v, se = se_rho, p = p_rho)
  cat(sprintf("SAR rho = %.4f | SE = %.4f | p = %.4g\n", rho_v, se_rho, p_rho))
}

f_cv <- reformulate(setdiff(all.vars(formules$M4), c("y", "region_acad")), response = "y")  # inclut z_moran
set.seed(1)
bloc <- kmeans(scale(coords), centers = 10)$cluster
err_sp <- sapply(sort(unique(bloc)), function(b) {
  te <- bloc == b
  m  <- glm(f_cv, family = poisson, data = d_sp[!te, ])
  mean((d_sp$y[te] - predict(m, d_sp[te, ], type = "response"))^2)
})
err_rd <- sapply(1:20, function(r) {
  set.seed(r)
  te <- sample(nrow(d_sp), round(0.10 * nrow(d_sp)))
  m  <- glm(f_cv, family = poisson, data = d_sp[-te, ])
  mean((d_sp$y[te] - predict(m, d_sp[te, ], type = "response"))^2)
})
tab_cv <- data.frame(methode = c("CV blocs spatiaux (k-means k=10)", "CV aléatoire (10%, 20 rép.)"),
                     mse = c(mean(err_sp), mean(err_rd)))
print(tab_cv)
write.csv(tab_cv, file.path(DOSSIER_OUT, "tab_cv.csv"), row.names = FALSE)
p_cv <- ggplot(tab_cv, aes(methode, mse, fill = methode)) +
  geom_col(width = 0.5) + geom_text(aes(label = round(mse, 3)), vjust = -0.5) +
  theme_minimal(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = "MSE", title = "CV spatial vs aléatoire : le CV aléatoire sous-estime l'erreur")
ggsave(file.path(DOSSIER_OUT, "fig_cv_communes.png"), p_cv, width = 7, height = 5, dpi = 150)
CK("4_spatial", list(tab_cv = tab_cv, moran = mt, m_sar = m_sar, bloc = bloc))
cat("Checkpoint 4 OK\n")

# ---------------------------------------------------------------------
# 5. SAUVEGARDE FINALE
# ---------------------------------------------------------------------
saveRDS(list(tab_spec = tab_spec, tab_boot = tab_boot, tab_mc = tab_mc, tab_cv = tab_cv,
             moran = mt, m_sar = m_sar, sar_rho = sar_rho, boot_diff = boot_diff, boot_moran = boot_moran, vdB_rob = vdB_rob,
             mc_raw = list(pois = mc_pois_dgp, nb = mc_nb_dgp),
             PARAM = PARAM, beta_star = beta_star),
        file.path(DOSSIER_OUT, "code7.rds"))
cat("Code 7 v2 terminé : code7.rds + tab_spec/tab_boot/tab_mc/tab_cv.csv + 3 figs + checkpoints\n")

