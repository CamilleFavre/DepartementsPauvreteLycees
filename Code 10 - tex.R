# =====================================================================
# CODE 9 - TABLES LATEX POUR LE .TEX (sortie console uniquement)
# A lancer APRES Code 4b v3 (modeles_comptage.rds) et Code 7 (code7.rds)
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

mods <- intersect(c("M0","M1","M2","M3","M3b","M3c","M3d","M4"), names(formules))

# ================= TABLEAU 1 : spécifications (généré depuis formules) ==========
section("TABLEAU 1 - Spécifications M0-M4", {
  allv <- setdiff(unique(unlist(lapply(formules[mods], all.vars))), c("y", "region_acad"))
  ord <- c("z_theil","z_moran","z_exposition_riches","z_getis_ord_riches","z_polarisation",
           "z_part_pop_top10","z_segregation_multidim","z_entropie_diversite","z_diversite_binaire",
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
          note = paste0("N = ", nrow(d_mod), " communes. FE = effets fixes région académique (13 modalités)."))
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
  tex_tab(desc, names(desc), "lrrrrr", "Statistiques descriptives (variables de M4)", "tab:desc")
})

# ================= TABLEAU 3 : distribution de y ================================
section("TABLEAU 3 - Distribution de y", {
  tb <- table(d_mod$y)
  dist <- data.frame(`$y$` = names(tb), Effectif = as.integer(tb), stringsAsFactors = FALSE)
  tex_tab(dist, c("$y$", "Effectif"), "lr",
          "Distribution du nombre de collèges privés par commune", "tab:dist",
          note = paste0("N = ", nrow(d_mod), " ; moyenne = ", f2(mean(d_mod$y)),
                        " ; variance = ", f2(var(d_mod$y)), " ; Var/Moy = ", f2(var(d_mod$y)/mean(d_mod$y)),
                        " ; zéros = ", f1(100*mean(d_mod$y == 0)), "\\% ; min = ", min(d_mod$y), " ; max = ", max(d_mod$y)))
})

# ================= TABLEAU 4 : régions ==========================================
section("TABLEAU 4 - Régions académiques", {
  reg <- as.data.frame(table(d_mod$region_acad), stringsAsFactors = FALSE)
  names(reg) <- c("Région", "n"); reg$Région <- esc(reg$Région)
  tex_tab(reg, c("Région académique", "n"), "lr", "Répartition par région académique", "tab:reg")
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
  vars_stab <- c("z_theil","z_moran","z_exposition_riches","z_polarisation",
                 "z_part_pop_top10","z_segregation_multidim")
  mods_stab <- intersect(mods, unique(tab_stab$modele))
  stab_w <- data.frame(Variable = var_label(vars_stab), stringsAsFactors = FALSE)
  for (m in mods_stab) {
    sel <- tab_stab$modele == m
    cf <- tab_stab$coef[sel][match(vars_stab, tab_stab$variable[sel])]
    pp <- tab_stab$p[sel][match(vars_stab, tab_stab$variable[sel])]
    stab_w[[m]] <- ifelse(is.na(cf), "--", paste0(f4(cf), stars(pp)))
  }
  tex_tab(stab_w, c("Variable", mods_stab), paste0("l", strrep("c", length(mods_stab))),
          "Coefficients de ségrégation selon le modèle", "tab:stab",
          note = "*** p<0.001, ** p<0.01, * p<0.05. z_Theil : +0.323*** en M0, non significatif en M4 (effet de composition).")
})

# ================= TABLEAU 8 : AIC / BIC ========================================
section("TABLEAU 8 - Comparaison des modèles", {
  ta <- data.frame(Modèle = esc(names(modeles)), AIC = f1(sapply(modeles, AIC)),
                   BIC = f1(sapply(modeles, BIC)), logLik = f1(sapply(modeles, logLik)),
                   stringsAsFactors = FALSE)
  tex_tab(ta, names(ta), "lrrr", "Critères d'information (AIC / BIC)", "tab:aic",
          note = "NB2 non convergent : AIC = 50 est un artefact. Meilleur modèle : Hurdle-Poisson (AIC = 2260.9).")
})

# ================= TABLEAU 9 : fréquences observées vs prédites =================
section("TABLEAU 9 - Fréquences observées vs prédites", {
  tf <- as.data.frame(tab_freq, stringsAsFactors = FALSE)
  if (!"y" %in% names(tf)) tf <- data.frame(y = rownames(tab_freq), tf, stringsAsFactors = FALSE)
  tf[] <- lapply(tf, function(x) if (is.numeric(x)) f1(x) else as.character(x))
  tex_tab(tf, c("$y$", esc(names(tf)[-1])), paste0("l", strrep("r", ncol(tf)-1)),
          "Fréquences observées et prédites (y = 0 à 8+)", "tab:freq",
          note = "Le hurdle prédit exactement les 353 zéros observés ; le Poisson en prévoit 392.")
})

# ================= TABLEAU 10 : AME ==============================================
section("TABLEAU 10 - AME des indicateurs de ségrégation (M4)", {
  ame_num <- function(m, data, var, h = 0.01) {
    d1 <- data; d0 <- data
    d1[[var]] <- d1[[var]] + h/2; d0[[var]] <- d0[[var]] - h/2
    mean((predict(m, d1, type = "response") - predict(m, d0, type = "response"))/h) }
  vars_ame <- intersect(c("z_theil","z_exposition_riches","z_polarisation",
                          "z_part_pop_top10","z_segregation_multidim"), all.vars(formules$M4))
  tame <- data.frame(Variable = var_label(vars_ame),
                     AME = f4(sapply(vars_ame, ame_num, m = m_p, data = d_mod)), stringsAsFactors = FALSE)
  tex_tab(tame, names(tame), "lr", "Effets marginaux moyens, M4 (différence finie h = 0.01)", "tab:ame",
          note = "AME en nombre de collèges privés par commune pour +1 écart-type de l'indicateur.")
})

# ================= TABLEAUX 11-14 : Code 7 (tests, bootstrap, MC, spatial) ======
if (!is.null(c7)) {
  
  section("TABLEAU 11 - Tests de spécification", {
    ct1 <- dispersiontest(m_p, trafo = 1)
    ts <- data.frame(
      Test = c("Cameron-Trivedi $\\alpha$ (trafo=2)", "Cameron-Trivedi (trafo=1, recalculé)",
               "Van den Broek (excès de zéros)", "Vuong Hurdle-Pois vs Poisson"),
      Statistique = c(f3(c7$tab_spec$stat[1]), f3(ct1$statistic), f3(c7$tab_spec$stat[2]), f2(c7$tab_spec$stat[3])),
      p = c(fmtp(c7$tab_spec$p[1]), fmtp(ct1$p.value), fmtp(c7$tab_spec$p[2]), fmtp(c7$tab_spec$p[3])),
      stringsAsFactors = FALSE)
    tex_tab(ts, names(ts), "lrr", "Tests de spécification", "tab:spectests",
            note = "CT : pas de surdispersion (Poisson OK). Van den Broek : excès de zéros (hurdle justifié). Vuong : hurdle > Poisson.")
  })
  
  section("TABLEAU 12 - Bootstrap cluster de l'écart beta(M4)-beta(M3)", {
    tb <- c7$tab_boot
    tb2 <- data.frame(Estimateur = esc(tb$estimateur), diff = f4(tb$diff), SE = f4(tb$se),
                      `IC95 bas` = f4(tb$ic_bas), `IC95 haut` = f4(tb$ic_haut), p = fmtp(tb$p),
                      stringsAsFactors = FALSE)
    tex_tab(tb2, names(tb2), "lrrrrr", "Écart de coefficient z_expo entre M3 et M4 (delta method vs bootstrap)", "tab:boot",
            note = "Bootstrap clusteré B = 500 (clusters = régions académiques). Écart non significatif : robuste aux FE.")
  })
  
  section("TABLEAU 13 - Monte Carlo (M = 200, beta vrai = 0.035)", {
    P <- c7$mc_raw$pois; N <- c7$mc_raw$nb
    s_mean <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else mean(x) }
    s_sd   <- function(x) { x <- x[!is.na(x)]; if (length(x) < 2) NA_real_ else sd(x) }
    s_cov  <- function(b, se) { ok <- !is.na(b) & !is.na(se)
    if (!any(ok)) NA_real_ else mean(abs(b[ok] - 0.035) < 1.96*se[ok]) }
    sfun <- function(M, dgp) data.frame(DGP = dgp, Estimateur = c("Poisson", "NB2"),
                                        Biais = c(s_mean(M[1,]) - 0.035, s_mean(M[3,]) - 0.035),
                                        SD_emp = c(s_sd(M[1,]), s_sd(M[3,])), SE_moy = c(s_mean(M[2,]), s_mean(M[4,])),
                                        Couv95 = c(s_cov(M[1,], M[2,]), s_cov(M[3,], M[4,])), stringsAsFactors = FALSE)
    tm <- rbind(sfun(P, "Poisson"), sfun(N, "NB2 (size = 1.5)"))
    n_na <- sum(is.na(P[3,])) + sum(is.na(N[3,]))
    tm2 <- data.frame(DGP = tm$DGP, Estimateur = tm$Estimateur, Biais = f4(tm$Biais),
                      `SD emp.` = f4(tm$SD_emp), `SE moyen` = f4(tm$SE_moy), `Couv. 95\\%` = f3(tm$Couv95),
                      stringsAsFactors = FALSE, check.names = FALSE)
    tex_tab(tm2, names(tm2), "llrrrr", "Monte Carlo : biais, SD, SE moyen, couverture (M = 200)", "tab:mc",
            note = paste0(n_na, " réplications NB2 non convergées (NA). Sous DGP Poisson : biais nul, couverture = 0.975."))
  })
  
  section("TABLEAU 14 - Spatial : Moran I + CV", {
    mor <- data.frame(
      Indicateur = c("Moran I des résidus de M4 (8-NN, poids W)",
                     "MSE CV par blocs spatiaux (k-means, k = 10)",
                     "MSE CV aléatoire (10\\%, 20 répétitions)"),
      Valeur = c(f4(c7$moran$estimate[1]), f3(c7$tab_cv$mse[1]), f3(c7$tab_cv$mse[2])),
      p = c(fmtp(c7$moran$p.value), "--", "--"), stringsAsFactors = FALSE)
    tex_tab(mor, names(mor), "lrr", "Autocorrélation spatiale des résidus et validation croisée", "tab:spatial",
            note = paste0("Le CV aléatoire sous-estime l'erreur de ",
                          f0(100*(c7$tab_cv$mse[1] - c7$tab_cv$mse[2])/c7$tab_cv$mse[1]),
                          "\\% (fuite spatiale). SAR non estimé (spatialreg absent) : Moran I + CV spatial suffisent."))
  })
} else cat("\n!! code7.rds illisible : relance Code 7 pour les tableaux 11-14.\n")

cat("\n\n>>> Copie chaque bloc entre '% ====' dans ton .tex, puis compile 2x (pdflatex) pour les \\ref / \\eqref.\n")