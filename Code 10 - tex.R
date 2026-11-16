# =====================================================================
# CODE 9 v3 - TABLES LATEX POUR LE .TEX (sortie console uniquement)
# A lancer APRES Code 4b v3 (modeles_comptage.rds, N=925 villes) et Code 7 v2 (code7.rds, Vuong corrigé)
# PATCHS v3 : en-têtes explicites (tabs 2, 9, 12, 13) ; régions en clair (tab 4) ;
#   y en entiers (tab 9) ; note CV spatiale en ratio (tab 14) ; vars_stab depuis rds ;
#   AME avec z_moran ; MC beta_star/PARAM dynamiques ; SAR conditionnel ; note PLM.
# Copie chaque bloc entre les "% ====" dans ton .tex, compile 2x.
# =====================================================================
library(arrow); library(dplyr); library(tidyr); library(sandwich)
library(lmtest); library(AER)          # PAS de MASS (masque dplyr::select)
set.seed(2026)

os_name <- Sys.info()[["sysname"]]
if (os_name == "Windows") {
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
} else if (os_name == "Darwin") {
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
} else stop("OS non pris en charge.")

charg <- readRDS(file.path(DOSSIER_OUT, "modeles_comptage.rds"))
c7    <- tryCatch(readRDS(file.path(DOSSIER_OUT, "code7.rds")), error = function(e) NULL)

d_mod <- charg$d_mod; modeles <- charg$modeles; formules <- charg$formules
tab_stab <- charg$tab_stab; tab_freq <- charg$tab_freq
m_p <- modeles$Poisson
vars_stab <- if (!is.null(charg$vars_stab)) charg$vars_stab else
  c("z_theil", "z_moran", "z_exposition_riches", "z_polarisation")
cat("N =", nrow(d_mod), "| vars_stab :", paste(vars_stab, collapse = ", "), "\n")
if (nrow(d_mod) != 925) cat("!! ATTENTION : rds pas au niveau ville (attendu N=925) — relance le 4b patché d'abord.\n")

H <- function(t) cat("\n\n############################################################\n## ", t,
                     "\n############################################################\n", sep = "")
section <- function(t, code) { H(t); tryCatch(code, error = function(e)
  cat("!! ERREUR section :", conditionMessage(e), "\n")) }

# ---------- helpers LaTeX ----------
esc  <- function(x) { x <- as.character(x)
x <- gsub("%", "\\\\%", x); x <- gsub("&", "\\\\&", x); x <- gsub("_", "\\\\_", x); x }
f0 <- function(x) ifelse(is.na(x), "--", sprintf("%.0f", x))
f1 <- function(x) ifelse(is.na(x), "--", sprintf("%.1f", x))
f2 <- function(x) ifelse(is.na(x), "--", sprintf("%.2f", x))
f3 <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
f4 <- function(x) ifelse(is.na(x), "--", sprintf("%.4f", x))
fmtp  <- function(p) ifelse(is.na(p), "--", ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p)))
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                                                          ifelse(p < 0.05, "*", ""))))

tex_tab <- function(df, col.names, align, cap, lab, note = NULL) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  cat("\n% ============================================================\n% TABLEAU : ", cap,
      "\n% ============================================================\n", sep = "")
  cat("\\begin{table}[htbp]\n  \\centering\n  \\small\n")
  cat("  \\caption{", cap, "}\n  \\label{", lab, "}\n", sep = "")
  cat("  \\begin{tabular}{", align, "}\n    \\toprule\n", sep = "")
  cat("    ", paste(col.names, collapse = " & "), " \\\\\n", sep = "")
  cat("    \\midrule\n")
  for (i in seq_len(nrow(df)))
    cat("    ", paste(as.character(unlist(df[i, ], use.names = FALSE)), collapse = " & "), " \\\\\n", sep = "")
  cat("    \\bottomrule\n")
  if (!is.null(note))
    cat("    \\multicolumn{", ncol(df), "}{l}{\\footnotesize ", note, "} \\\\\n", sep = "")
  cat("  \\end{tabular}\n\\end{table}\n")
  cat("\n% -- vérif rapide (texte) --\n"); print(df, row.names = FALSE)
}

lab_z <- c(y = "$y$",
           z_theil = "$z_{\\text{Theil}}$", z_moran = "$z_{\\text{Moran}}$",
           z_exposition_riches = "$z_{\\text{expo}}$", z_getis_ord_riches = "$z_{\\text{Getis}}$",
           z_polarisation = "$z_{\\text{pol}}$", z_part_pop_top10 = "$z_{\\text{top10}}$",
           z_segregation_multidim = "$z_{\\text{multi}}$", z_entropie_diversite = "$z_{\\text{entropie}}$",
           z_diversite_binaire = "$z_{\\text{div}}$",
           log_pop = "$\\log pop$", log_rev = "$\\log rev$",
           taux_pauv_commune = "$\\text{tauxPauv}$", nb_colleges_publics = "$\\text{nbPub}$",
           dist_college_public_km = "$\\text{distPub}$", nb_prives_5km_autres = "$\\text{priv5km}$")
var_label <- function(v) {
  v <- as.character(v); out <- rep(NA_character_, length(v))
  hit <- v %in% names(lab_z); out[hit] <- unname(lab_z[v[hit]])
  r <- v[!hit]
  out[!hit] <- ifelse(r == "(Intercept)", "$\\beta_0$ (intercept)",
                      ifelse(startsWith(r, "region_acad"),
                             paste0("\\textit{FE ", esc(sub("^region_acad", "", r)), "}"), esc(r)))
  out
}

mods <- intersect(c("M0","M1","M2","M3","M3b","M4"), names(formules))
note_plm <- "Paris, Lyon et Marseille regroupées en une seule ville (75056 / 69123 / 13055)."

# ================= TABLEAU 1 : spécifications (généré depuis formules) ==========
section("TABLEAU 1 - Spécifications M0-M4", {
  allv <- setdiff(unique(unlist(lapply(formules[mods], all.vars))), c("y", "region_acad"))
  ord <- c("z_theil","z_moran","z_exposition_riches","z_polarisation",
           "log_pop","log_rev","taux_pauv_commune","nb_colleges_publics",
           "dist_college_public_km","nb_prives_5km_autres")
  allv <- c(intersect(ord, allv), setdiff(allv, ord))
  spec <- data.frame(Variable = var_label(allv), stringsAsFactors = FALSE)
  for (m in mods) spec[[m]] <- ifelse(allv %in% all.vars(formules[[m]]), "$\\times$", "--")
  aic_mods <- sapply(mods, function(m) tryCatch(AIC(glm(formules[[m]], family = poisson, data = d_mod)),
                                                error = function(e) NA_real_))
  aic_row <- data.frame(Variable = "AIC (Poisson)", stringsAsFactors = FALSE)
  fe_row  <- data.frame(Variable = "$\\alpha_{r(i)}$ (FE région académique)", stringsAsFactors = FALSE)
  for (m in mods) { aic_row[[m]] <- f1(aic_mods[[m]]); fe_row[[m]] <- if (m == "M4") "$\\times$" else "--" }
  spec <- rbind(spec, aic_row, fe_row)
  tex_tab(spec, c("Variable", mods), paste0("l", strrep("c", length(mods))),
          "Spécifications M0--M4 (Poisson, log-link)", "tab:spec",
          note = paste0("N = ", nrow(d_mod), " villes. FE = effets fixes région académique (13 modalités). ", note_plm))
})

# ================= TABLEAU 2 : descriptives ======================================
section("TABLEAU 2 - Descriptives", {
  vars_m4 <- setdiff(all.vars(formules$M4), c("y", "region_acad"))
  desc <- data.frame(Variable = var_label(vars_m4),
                     Moyenne = f2(sapply(d_mod[vars_m4], mean, na.rm = TRUE)),
                     `Écart-type` = f2(sapply(d_mod[vars_m4], sd, na.rm = TRUE)),
                     Min = f2(sapply(d_mod[vars_m4], min, na.rm = TRUE)),
                     Médiane = f2(sapply(d_mod[vars_m4], median, na.rm = TRUE)),
                     Max = f2(sapply(d_mod[vars_m4], max, na.rm = TRUE)), stringsAsFactors = FALSE)
  # PATCH v3 : en-têtes explicites (sinon check.names produit "Écart.type")
  tex_tab(desc, c("Variable", "Moyenne", "Écart-type", "Min", "Médiane", "Max"), "lrrrrr",
          "Statistiques descriptives (variables de M4)", "tab:desc",
          note = paste0("N = ", nrow(d_mod), " villes. ", note_plm))
})

# ================= TABLEAU 3 : distribution de y ================================
section("TABLEAU 3 - Distribution de y", {
  tb <- table(d_mod$y)
  dist <- data.frame(`$y$` = names(tb), Effectif = as.integer(tb), stringsAsFactors = FALSE)
  tex_tab(dist, c("$y$", "Effectif"), "lr",
          "Distribution du nombre de collèges privés par ville", "tab:dist",
          note = paste0("N = ", nrow(d_mod), " ; moyenne = ", f2(mean(d_mod$y)),
                        " ; variance = ", f2(var(d_mod$y)), " ; Var/Moy = ", f2(var(d_mod$y)/mean(d_mod$y)),
                        " ; zéros = ", f1(100*mean(d_mod$y == 0)), "\\% ; min = ", min(d_mod$y), " ; max = ", max(d_mod$y),
                        " (Paris). ", note_plm))
})

# ================= TABLEAU 4 : régions (noms en clair) ==========================
section("TABLEAU 4 - Régions académiques", {
  reg <- as.data.frame(table(d_mod$region_acad), stringsAsFactors = FALSE)
  names(reg) <- c("code", "n")
  pretty <- c(auvergne_rhone_alpes = "Auvergne-Rhône-Alpes",
              bourgogne_franche_comte = "Bourgogne-Franche-Comté",
              bretagne = "Bretagne", centre_val_de_loire = "Centre-Val de Loire",
              corse = "Corse", grand_est = "Grand Est",
              hauts_de_france = "Hauts-de-France", ile_de_france = "Île-de-France",
              normandie = "Normandie", nouvelle_aquitaine = "Nouvelle-Aquitaine",
              occitanie = "Occitanie", pays_de_la_loire = "Pays de la Loire",
              provence_alpes_cote_azur = "Provence-Alpes-Côte d'Azur")
  stopifnot(all(as.character(reg$code) %in% names(pretty)))
  reg$Région <- unname(pretty[as.character(reg$code)])
  tex_tab(reg[, c("Région", "n")], c("Région académique", "n"), "lr",
          "Répartition par région académique", "tab:reg",
          note = paste0("N = ", nrow(d_mod), " villes."))
})

# ================= TABLEAUX 5-6 : coefficients M4 + FE + LRT ====================
section("TABLEAUX 5-6 - Coefficients M4 (SE cluster, IRR, IC95) et FE", {
  V_cl <- vcovCL(m_p, cluster = d_mod$region_acad, type = "HC1")
  ct_p <- coeftest(m_p, vcov. = V_cl)
  ct_num <- data.frame(Variable = rownames(ct_p), beta = ct_p[,1], se = ct_p[,2],
                       z = ct_p[,3], p = ct_p[,4], irr = exp(ct_p[,1]),
                       icb = exp(ct_p[,1] - 1.96*ct_p[,2]), ich = exp(ct_p[,1] + 1.96*ct_p[,2]))
  fmt_coef <- function(d) data.frame(Variable = var_label(d$Variable), `$\\beta$` = f4(d$beta),
                                     SE = f4(d$se), z = f2(d$z), p = fmtp(d$p), IRR = f3(d$irr),
                                     `IC95 (IRR)` = paste0("$[", f3(d$icb), " ; ", f3(d$ich), "]$"), stringsAsFactors = FALSE)
  m3  <- glm(formules$M3, family = poisson, data = d_mod)
  lrt <- anova(m3, m_p, test = "Chisq")
  lr_stat <- lrt$Deviance[2]; lr_df <- lrt$Df[2]; lr_p <- lrt[["Pr(>Chi)"]][2]
  note5 <- paste0("SE clusterisées par région académique (HC1), N = ", nrow(d_mod),
                  ". LRT des FE (M3 vs M4) : $\\chi^2(", lr_df, ") = ", f1(lr_stat), "$, p ", fmtp(lr_p), ".")
  hf <- ct_num[!grepl("^region_acad", ct_num$Variable), ]
  fe <- ct_num[ grepl("^region_acad", ct_num$Variable), ]; fe <- fe[order(-fe$irr), ]
  tex_tab(fmt_coef(hf), c("Variable","$\\beta$","SE","z","p","IRR","IC95 (IRR)"), "lrrrrrr",
          "Poisson M4 : coefficients (hors effets fixes régionaux)", "tab:coefM4", note = note5)
  tex_tab(fmt_coef(fe), c("Variable","$\\beta$","SE","z","p","IRR","IC95 (IRR)"), "lrrrrrr",
          "Poisson M4 : effets fixes régionaux (triés par IRR décroissant)", "tab:feM4",
          note = "Région de référence absorbée par l'intercept.")
})

# ================= TABLEAU 7 : stabilité M0 -> M4 ===============================
section("TABLEAU 7 - Stabilité des coefficients de ségrégation", {
  mods_stab <- intersect(mods, unique(tab_stab$modele))
  stab_w <- data.frame(Variable = var_label(vars_stab), stringsAsFactors = FALSE)
  for (m in mods_stab) {
    sel <- tab_stab$modele == m
    cf <- tab_stab$coef[sel][match(vars_stab, tab_stab$variable[sel])]
    pp <- tab_stab$p[sel][match(vars_stab, tab_stab$variable[sel])]
    stab_w[[m]] <- ifelse(is.na(cf), "--", paste0(f4(cf), stars(pp)))
  }
  get_stab <- function(m, v) { sel <- tab_stab$modele == m & tab_stab$variable == v
  if (!any(sel)) return(NA_real_) else tab_stab$coef[sel][1] }
  note7 <- paste0("*** p<0.001, ** p<0.01, * p<0.05. z_Theil : ", f4(get_stab("M0","z_theil")), stars(tab_stab$p[tab_stab$modele=="M0" & tab_stab$variable=="z_theil"]),
                  " en M0, ", f4(get_stab("M4","z_theil")), stars(tab_stab$p[tab_stab$modele=="M4" & tab_stab$variable=="z_theil"]),
                  " en M4 (effet de composition absorbé par log pop). z_Moran : ", f4(get_stab("M4","z_moran")),
                  stars(tab_stab$p[tab_stab$modele=="M4" & tab_stab$variable=="z_moran"]), " en M4.")
  tex_tab(stab_w, c("Variable", mods_stab), paste0("l", strrep("c", length(mods_stab))),
          "Coefficients de ségrégation selon le modèle", "tab:stab", note = note7)
})

# ================= TABLEAU 8 : AIC / BIC ========================================
section("TABLEAU 8 - Comparaison des modèles", {
  aics <- sapply(modeles, AIC)
  nb2_deg <- !is.null(modeles$NB2) && (is.null(modeles$NB2$theta) || !is.finite(modeles$NB2$theta) || modeles$NB2$theta > 1e6)
  zinb_sing <- !is.null(modeles$ZINB) && any(is.na(coef(modeles$ZINB)))
  aics_ok <- aics; if (nb2_deg && "NB2" %in% names(aics_ok)) aics_ok["NB2"] <- NA_real_
  best <- names(which.min(aics_ok))
  ta <- data.frame(Modèle = esc(names(modeles)), AIC = f1(aics),
                   BIC = f1(sapply(modeles, BIC)), logLik = f1(sapply(modeles, logLik)),
                   stringsAsFactors = FALSE)
  note8 <- paste0(
    if (nb2_deg) "NB2 dégénéré (theta -> inf, coefficients = Poisson) : son AIC est un artefact numérique, modèle écarté. " else "",
    if (zinb_sing) "ZINB numériquement singulier (SE = NA), modèle écarté. " else "",
    "Meilleur modèle retenu : ", best, " (AIC = ", f1(aics_ok[best]), ").")
  tex_tab(ta, c("Modèle", "AIC", "BIC", "logLik"), "lrrr", "Critères d'information (AIC / BIC)", "tab:aic", note = note8)
})

# ================= TABLEAU 9 : fréquences observées vs prédites (y entiers) ====
section("TABLEAU 9 - Fréquences observées vs prédites", {
  tf <- as.data.frame(tab_freq, stringsAsFactors = FALSE)
  if (!"y" %in% names(tf)) tf <- data.frame(y = rownames(tab_freq), tf, stringsAsFactors = FALSE)
  z_obs <- tf$Observe[which(as.numeric(tf$y) == 0)]
  z_pois <- if ("Poisson" %in% names(tf)) tf$Poisson[which(as.numeric(tf$y) == 0)] else NA
  z_hurp <- if ("Hurdle_Pois" %in% names(tf)) tf$Hurdle_Pois[which(as.numeric(tf$y) == 0)] else NA
  y_int <- as.integer(as.numeric(tf$y))                       # PATCH v3 : y en entiers
  tf[] <- lapply(tf, function(x) if (is.numeric(x)) f1(x) else as.character(x))
  tf$y <- as.character(y_int)
  tex_tab(tf, c("$y$", esc(sub("^Observe$", "Observé", names(tf)[-1]))), paste0("l", strrep("r", ncol(tf)-1)),
          "Fréquences observées et prédites (y = 0 à 8+)", "tab:freq",
          note = paste0("Zéros observés : ", f0(z_obs), " ; prédits par Hurdle_Pois : ", f0(z_hurp),
                        " ; par Poisson : ", f0(z_pois), "."))
})

# ================= TABLEAU 10 : AME (4 indices, dont z_moran) ===================
section("TABLEAU 10 - AME des indicateurs de ségrégation (M4)", {
  ame_num <- function(m, data, var, h = 0.01) {
    d1 <- data; d0 <- data
    d1[[var]] <- d1[[var]] + h/2; d0[[var]] <- d0[[var]] - h/2
    mean((predict(m, d1, type = "response") - predict(m, d0, type = "response"))/h) }
  vars_ame <- intersect(vars_stab, all.vars(formules$M4))   # inclut z_moran
  tame <- data.frame(Variable = var_label(vars_ame),
                     AME = f4(sapply(vars_ame, ame_num, m = m_p, data = d_mod)), stringsAsFactors = FALSE)
  tex_tab(tame, c("Variable", "AME"), "lr", "Effets marginaux moyens, M4 (différence finie h = 0.01)", "tab:ame",
          note = "AME en nombre de collèges privés par ville pour +1 écart-type de l'indicateur (z-score).")
})

# ================= TABLEAUX 11-14 : Code 7 (tests, bootstrap, MC, spatial) ======
section("TABLEAU 14 - Spatial : Moran I (Poisson & Hurdle) et rho SAR (2 poids)", {
  # --- repli sur le checkpoint 4_spatial si code7.rds n'a pas encore les champs IDW/hurdle ---
  ck4 <- tryCatch(readRDS(file.path(DOSSIER_OUT, "code7_4_spatial.rds")), error = function(e) NULL)
  pick <- function(a, b) if (!is.null(a)) a else if (!is.null(b)) b else NULL
  moran_idw   <- pick(c7$moran_idw,   if (!is.null(ck4)) ck4$moran_idw   else NULL)
  m_sar_idw   <- pick(c7$m_sar_idw,   if (!is.null(ck4)) ck4$m_sar_idw   else NULL)
  sar_rho_idw <- pick(c7$sar_rho_idw, if (!is.null(ck4)) ck4$sar_rho_idw else NULL)
  mh_knn      <- pick(c7$moran_h_knn, if (!is.null(ck4)) ck4$moran_h_knn else NULL)
  mh_idw      <- pick(c7$moran_h_idw, if (!is.null(ck4)) ck4$moran_h_idw else NULL)
  has_idw <- !is.null(moran_idw)
  
  # extraction robuste de rho (grep) + p depuis sar_rho si dispo
  get_rho <- function(m_sar, sar_rho) {
    if (!is.null(sar_rho) && all(c("rho", "p") %in% names(sar_rho)))
      return(c(rho = unname(sar_rho["rho"]), p = unname(sar_rho["p"])))
    if (is.null(m_sar)) return(c(rho = NA_real_, p = NA_real_))
    cf <- coef(m_sar)
    i <- grep("^rho$", names(cf)); if (!length(i)) i <- grep("rho", names(cf))
    if (!length(i)) return(c(rho = NA_real_, p = NA_real_))
    rho_v <- unname(cf[i[1]])
    V <- tryCatch(vcov(m_sar), error = function(e) NULL)
    se <- NA_real_
    if (!is.null(V)) {
      j <- grep("^rho$", rownames(V)); if (!length(j)) j <- grep("rho", rownames(V))
      if (length(j)) se <- sqrt(V[j[1], j[1]])
    }
    p <- if (is.finite(se) && se > 0) 2 * pnorm(-abs(rho_v / se)) else NA_real_
    c(rho = rho_v, p = p)
  }
  r_knn <- get_rho(c7$m_sar, c7$sar_rho)
  r_idw <- if (has_idw) get_rho(m_sar_idw, sar_rho_idw) else c(rho = NA_real_, p = NA_real_)
  
  w_idw <- if (has_idw) "1/d, $\\leq$ 200 km, W" else "--"
  mor <- data.frame(
    Indicateur = c("Moran $I$ résidus Poisson M4", "Moran $I$ résidus Poisson M4",
                   "Moran $I$ résidus Hurdle M4", "Moran $I$ résidus Hurdle M4",
                   "$\\rho$ SAR (lagsarlm)", "$\\rho$ SAR (lagsarlm)"),
    Poids = c("8-NN, W", w_idw, "8-NN, W", w_idw, "8-NN, W", w_idw),
    Valeur = c(f4(c7$moran$estimate[1]),
               if (has_idw) f4(moran_idw$estimate[1]) else "--",
               if (!is.null(mh_knn)) f4(mh_knn$estimate[1]) else "--",
               if (has_idw && !is.null(mh_idw)) f4(mh_idw$estimate[1]) else "--",
               f4(r_knn["rho"]),
               if (has_idw) f4(r_idw["rho"]) else "--"),
    p = c(fmtp(c7$moran$p.value),
          if (has_idw) fmtp(moran_idw$p.value) else "--",
          if (!is.null(mh_knn)) fmtp(mh_knn$p.value) else "--",
          if (has_idw && !is.null(mh_idw)) fmtp(mh_idw$p.value) else "--",
          fmtp(r_knn["p"]),
          if (has_idw) fmtp(r_idw["p"]) else "--"),
    stringsAsFactors = FALSE)
  
  note14 <- paste0(
    "Moran $I$ sur les résidus (Poisson : $y-\\hat\\mu$ via lm.morantest ; Hurdle : résidus response $y-\\hat{E}[y]$ via moran.test). ",
    "Poisson : 8-NN $I$ = ", f4(c7$moran$estimate[1]), " (p = ", sprintf("%.1e", c7$moran$p.value), ")",
    if (has_idw) paste0(" ; 1/d $\\leq$ 200 km $I$ = ", f4(moran_idw$estimate[1]),
                        " (p = ", sprintf("%.1e", moran_idw$p.value), ")") else "",
    if (!is.null(mh_knn)) paste0(". Hurdle : 8-NN $I$ = ", f4(mh_knn$estimate[1]), " (p = ", fmtp(mh_knn$p.value), ")",
                                 if (has_idw && !is.null(mh_idw)) paste0(" ; 1/d $I$ = ", f4(mh_idw$estimate[1]), " (p = ", fmtp(mh_idw$p.value), ")") else "") else "",
    " : autocorrélation résiduelle positive sous les deux poids, d'où SE clusterisées par département et CV par blocs. ",
    "SAR $\\rho$ : 8-NN $\\rho$ = ", f4(r_knn["rho"]), " (p = ", fmtp(r_knn["p"]), ")",
    if (has_idw) paste0(" ; 1/d $\\leq$ 200 km $\\rho$ = ", f4(r_idw["rho"]), " (p = ", fmtp(r_idw["p"]), ")") else "",
    " : en 8-NN la dépendance de $y$ est absorbée par les FE régionales ; sous 1/d dense, $\\rho<0$ = ajustement intra-régional (non causal).")
  tex_tab(mor, c("Indicateur", "Poids $w_{ij}$", "Valeur", "p"), "llrr",
          "Autocorrélation spatiale des résidus (Moran) et SAR, 2 poids", "tab:spatial", note = note14)
})

section("TABLEAU 15 - CV blocs vs aléatoire : Poisson/Hurdle, moyenne des folds ET pondérée par ville", {
  cv <- c7$tab_cv
  gcol <- function(nm) if (nm %in% names(cv)) cv[[nm]] else rep(NA_real_, nrow(cv))
  has_w   <- "mse_ponderee" %in% names(cv)     # pondération par taille de bloc (nouveau code7.rds)
  has_mae <- "mae" %in% names(cv)
  nt  <- c7$naive_tab
  gnt <- function(nm, k) if (!is.null(nt) && nm %in% names(nt)) nt[[nm]][k] else NA_real_
  
  cv2 <- data.frame(
    Modele    = cv$modele,
    Decoupage = cv$decoupage,
    MSE       = f3(cv$mse),
    MSE_pond  = if (has_w) f3(gcol("mse_ponderee")) else "--",
    MAE       = if (has_mae) f3(gcol("mae")) else "--",
    ratio     = if ("ratio_vs_blocs" %in% names(cv)) f2(gcol("ratio_vs_blocs")) else "--",
    skill     = if ("skill_vs_naive" %in% names(cv)) f3(gcol("skill_vs_naive")) else "--",
    skill_p   = if (has_w && "skill_vs_naive_pond" %in% names(cv)) f3(gcol("skill_vs_naive_pond")) else "--",
    stringsAsFactors = FALSE)
  
  note15 <- paste0(
    "MSE = moyenne des MSE par fold (10 blocs k-means k=10, ou 20 folds aléatoires 10\\%) ; ",
    "MSE pond. = MSE poolée par ville, $\\sum_b n_b\\,\\mathrm{MSE}_b/\\sum_b n_b$ (poids = taille du fold) : ",
    "sans effet en CV aléatoire (folds de taille constante 92), décisif en blocs (tailles 34--299). ",
    "Skill = $1-\\mathrm{MSE}/\\mathrm{MSE}_{naive}$, baseline naïve (moyenne du train, ou prédiction 0) évaluée dans les MÊMES folds : ",
    "naive moyenne = ", f3(gnt("naive_moyenne", 1)), " (blocs) / ", f3(gnt("naive_moyenne", 2)), " (aléa)",
    if ("naive_moyenne_pond" %in% names(nt)) paste0(" [pond. blocs = ", f3(gnt("naive_moyenne_pond", 1)), "]") else "",
    " ; naive zéro = ", f3(gnt("naive_zero", 1)), " / ", f3(gnt("naive_zero", 2)), ". ",
    "Ratio aléa/blocs (Poisson) = ", f2(cv$mse[2]/cv$mse[1]),
    " : la CV aléatoire est dominée par l'outlier Paris (y = 103) présent dans le test + fuite spatiale ; ",
    "la CV par blocs spatiaux reste la référence.")
  if (!is.null(c7$cmp)) {
    cb <- c7$cmp[which.max(c7$cmp$diff), ]
    note15 <- paste0(note15, " Diagnostic par bloc (champ cmp) : écart Hurdle$-$Poisson maximal dans le bloc ",
                     cb$bloc, " (n = ", cb$n, ", diff = ", f3(cb$diff), ") = bloc de Paris",
                     if (!is.null(c7$tst) && nrow(c7$tst) > 0)
                       paste0(" ; dans ce fold la ligne Paris (y = 103) pèse sq_pois = ", f2(c7$tst$sq_pois[1]),
                              " vs sq_hurd = ", f2(c7$tst$sq_hurd[1]), " (tst trié par sq_hurd décroissant).") else ".")
  }
  tex_tab(cv2, c("Modèle", "Découpage", "MSE", "MSE pond.", "MAE", "ratio", "skill", "skill pond."), "llrrrrrr",
          "CV par blocs spatiaux vs aléatoire, Poisson et Hurdle (moyenne des folds et pondérée par ville)", "tab:cv", note = note15)
})