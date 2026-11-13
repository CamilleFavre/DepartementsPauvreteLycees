# =====================================================================
# FIG - Nuage : taux de collèges privés pour 10 000 hab. (Y)
#        vs indice de Moran brut du revenu (X), COULEUR = classe de
#        polarisation (W) en intervalles D1 / Q1 / Méd / Q3 / D9
#        Forme imitée du modèle : fond blanc, axes noirs + ticks, pas de
#        grille, pas de cadre, pas de titre, points ronds alpha<1,
#        légende verticale à droite, palette bleu->rouge du modèle
# =====================================================================
rm(list = ls())
library(dplyr); library(ggplot2); library(RColorBrewer)

os_name <- Sys.info()[["sysname"]]
if (os_name == "Windows") {
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
} else {
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
}

select <- dplyr::select

d_mod <- readRDS(file.path(DOSSIER_OUT, "modeles_comptage.rds"))$d_mod


d_fig <- d_mod[d_mod$y > 0 & d_mod$pop_age_11_17 > 100, ] %>%
  filter(is.finite(moran_I_rev), is.finite(y)) %>%
  mutate(taux_prive_1k = y * 1000 / pop_age_11_17)

# --- seuils D1, Q1, Méd, Q3, D9 de la polarisation ---
q  <- quantile(d_fig$polarisation, probs = c(.10, .25, .50, .75, .90), na.rm = TRUE)
D1 <- q[1]; Q1 <- q[2]; MED <- q[3]; Q3 <- q[4]; D9 <- q[5]
cat("D1 =", round(D1,3), "| Q1 =", round(Q1,3), "| Méd =", round(MED,3),
    "| Q3 =", round(Q3,3), "| D9 =", round(D9,3), "\n")

# --- 6 classes (ordre croissant), puis levels inversés pour que "> D9" soit EN HAUT de la légende ---
lab <- c("< D1", "D1 - Q1", "Q1 - Méd", "Méd - Q3", "Q3 - D9", "> D9")
d_fig$pol_classe <- cut(d_fig$polarisation,
                        breaks = c(-Inf, D1, Q1, MED, Q3, D9, Inf),
                        labels = lab, include.lowest = TRUE)
d_fig$pol_classe <- factor(d_fig$pol_classe, levels = rev(lab))   # haut de légende = classe haute
print(table(d_fig$pol_classe))

# --- palette du modèle (haut = bleu foncé, bas = rouge), 6 couleurs ---
cols <- c("#08519c", "#4393c3", "#92c5de", "#fee090", "#f46d43", "#a50026")  # >D9 -> ... -> <D1
names(cols) <- levels(d_fig$pol_classe)

p <- ggplot(d_fig, aes(x = moran_I_rev, y = taux_prive_1k, colour = pol_classe)) +
  geom_point(shape = 16, size = 2.2, alpha = 0.65, stroke = 0) +
  scale_colour_manual(values = cols,
                      name = "Polarisation (indice W),\nintervalles D1 / Q1 / Méd / Q3 / D9 :") +
  labs(x = "Indice de Moran du revenu (brut)",
       y = "Collèges privés sous contrat\npour 100 habitants de 11 à 17 ans") +
  theme_classic(base_size = 12) +          # fond blanc, axes noirs + ticks, pas de grille/cadre
  theme(legend.position = "right",          # légende verticale à droite, clés = points ronds
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

print(p)




# =====================================================================
# FIG 2 - INVERSÉ : x = polarisation (brut), couleur = Moran (brut, dégradé continu)
# =====================================================================
p2 <- ggplot(d_fig, aes(x = polarisation, y = taux_prive_1k, colour = moran_I_rev)) +
  geom_point(shape = 16, size = 2.2, alpha = 0.65, stroke = 0) +
  scale_colour_gradient2(low = "#a50026", mid = "#f7f7f7", high = "#08519c",
                         midpoint = 0,
                         name = "Indice de Moran\ndu revenu (brut)",
                         guide = "colourbar") +          # légende à droite en barre continue
  labs(x = "Polarisation (indice W, brut)",
       y = "Collèges privés sous contrat\npour 1 000 habitants de 11 à 17 ans") +
  theme_classic(base_size = 12) +          # même forme que le modèle / fig 1
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

print(p2)



# --- régression y ~ x (une seule droite, tous points confondus) ---
m_reg <- lm(taux_prive_1k ~ moran_I_rev + I(moran_I_rev)^2, data = d_fig)
print(summary(m_reg)$coefficients)     # pente, SE, t, p-value
cat("R² =", round(summary(m_reg)$r.squared, 4), " | N =", nrow(d_fig), "\n")

p <- ggplot(d_fig, aes(x = moran_I_rev, y = taux_prive_1k, colour = pol_classe)) +
  geom_point(shape = 16, size = 2.2, alpha = 0.65, stroke = 0) +
  geom_smooth(aes(colour = NULL),                          # <- UNE seule courbe, pas par classe
              method = "lm", formula = y ~ x + I(x^2),     # <- même spec quadratique que le lm
              se = TRUE, colour = "black", fill = "grey60",
              linewidth = 0.6, alpha = 0.15) +
  scale_colour_manual(values = cols,
                      name = "Polarisation (indice W),\nintervalles D1 / Q1 / Méd / Q3 / D9 :") +
  labs(x = "Indice de Moran du revenu (brut)",
       y = "Collèges privés sous contrat\npour 100 habitants de 11 à 17 ans") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

print(p)



# =====================================================================
# FIG 2 - INVERSÉ : x = polarisation (brut), COULEUR = classes du Moran brut
#                 (mêmes classes D1/Q1/Méd/Q3/D9 et mêmes couleurs que la fig 1)
#                 + régression y ~ x (une seule droite)
# =====================================================================
qm   <- quantile(d_fig$moran_I_rev, probs = c(.10, .25, .50, .75, .90), na.rm = TRUE)
mD1 <- qm[1]; mQ1 <- qm[2]; mMED <- qm[3]; mQ3 <- qm[4]; mD9 <- qm[5]
cat("Moran : D1 =", round(mD1,3), "| Q1 =", round(mQ1,3), "| Méd =", round(mMED,3),
    "| Q3 =", round(mQ3,3), "| D9 =", round(mD9,3), "\n")

d_fig$moran_classe <- cut(d_fig$moran_I_rev,
                          breaks = c(-Inf, mD1, mQ1, mMED, mQ3, mD9, Inf),
                          labels = lab, include.lowest = TRUE)   # mêmes labels que pol_classe
d_fig$moran_classe <- factor(d_fig$moran_classe, levels = rev(lab))  # classe haute EN HAUT
print(table(d_fig$moran_classe))




# --- régression y ~ x (ici x = polarisation) ---
m_reg2 <- lm(taux_prive_1k ~ polarisation, data = d_fig)
print(summary(m_reg2)$coefficients)
cat("R² =", round(summary(m_reg2)$r.squared, 4), " | N =", nrow(d_fig), "\n")

p2 <- ggplot(d_fig, aes(x = polarisation, y = taux_prive_1k, colour = moran_classe)) +
  geom_point(shape = 16, size = 2.2, alpha = 0.65, stroke = 0) +
  geom_smooth(aes(colour = NULL),                       # <- UNE seule droite, pas par classe
              method = "lm", se = TRUE,
              colour = "black", fill = "grey60",
              linewidth = 0.6, alpha = 0.15) +
  scale_colour_manual(values = cols,                    # <- MÊMES 6 couleurs que la fig 1
                      name = "Moran du revenu (brut),\nintervalles D1 / Q1 / Méd / Q3 / D9 :") +
  labs(x = "Polarisation (indice W, brut)",
       y = "Collèges privés sous contrat\npour 1 000 habitants de 11 à 17 ans") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 10),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA))

print(p2)
