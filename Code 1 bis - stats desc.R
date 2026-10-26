library(arrow)
library(dplyr)

# Détection automatique du système d'exploitation
os_name <- Sys.info()[["sysname"]]

# Attribution des chemins selon l'OS
if (os_name == "Windows") {
  DOSSIER           <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT       <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
  CHEMIN_CENTROIDES <- "C:/Users/camil/Downloads/centroides-communes-population.geojson"
  
} else if (os_name == "Darwin") { # macOS
  DOSSIER           <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT       <- "/Users/camillefavre/DepartementsPauvreteLycees"
  CHEMIN_CENTROIDES <- "/Users/camillefavre/Downloads/centroides-communes-population.geojson"
  
} else {
  stop("Système d'exploitation non pris en charge.")
}


table_lycees <- read_parquet(file.path(DOSSIER_OUT, "df_analyse.parquet")) %>%
  transmute(
    code_insee,
    nom_commune,
    nb_lycees_prives,
    nb_lycees_publics,
    habitants = pop_totale,
    # Effectifs de lycéens (voie générale et technologique)
    eleves_prives = eff_gt_prive,
    eleves_publics = eff_gt_public,
    # Ratports offre
    prive_par_hab     = ifelse(habitants > 0, nb_lycees_prives / habitants, NA_real_),
    prive_pour_10k_hab = ifelse(habitants > 0, 1e4 * nb_lycees_prives / habitants, NA_real_),
    prive_par_public  = ifelse(nb_lycees_publics > 0,
                               nb_lycees_prives / nb_lycees_publics, NA_real_),
    part_privee_offre = ifelse(nb_lycees_prives + nb_lycees_publics > 0,
                               nb_lycees_prives / (nb_lycees_prives + nb_lycees_publics), NA_real_),
    # Ratios élèves
    eleve_prive_par_public = ifelse(eleves_publics > 0,
                                    eleves_prives / eleves_publics, NA_real_),
    part_eleves_prives     = ifelse(eleves_prives + eleves_publics > 0,
                                    eleves_prives / (eleves_prives + eleves_publics), NA_real_),
    eleves_prives_pour_10k_hab = ifelse(habitants > 0,
                                        1e4 * eleves_prives / habitants, NA_real_)
  ) %>%
  arrange(desc(nb_lycees_prives))

head(table_lycees, 20)
