# =====================================================================
# app.R v2 — Dashboard pro : collèges privés & ségrégation urbaine
# Lancer dans RStudio : bouton "Run App"  (ou shiny::runApp("."))
# Packages : shiny bslib bsicons leaflet DT plotly sf arrow dplyr tidyr
#            tibble ggplot2 sandwich lmtest cachem
# =====================================================================
library(shiny); library(bslib); library(bsicons)
library(dplyr); library(tidyr); library(tibble)
library(sf); library(arrow); library(ggplot2)
library(leaflet); library(DT); library(plotly)
library(sandwich); library(lmtest); library(cachem)

# ---- Chemins ----
os_name <- Sys.info()[["sysname"]]
if (os_name == "Windows") {
  DOSSIER     <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "C:/Users/camil/Documents/DepartementsPauvreteLycees"
} else if (os_name == "Darwin") {
  DOSSIER     <- "/Users/camillefavre/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics"
  DOSSIER_OUT <- "/Users/camillefavre/DepartementsPauvreteLycees"
} else stop("OS non pris en charge")

CHEMIN_CARREAUX <- file.path(DOSSIER, "france_metro", "carreaux_200m_met.gpkg")

# ---- Données (chargées une seule fois) ----
master   <- read_parquet(file.path(DOSSIER_OUT, "df_dashboard.parquet"))
charg    <- readRDS(file.path(DOSSIER_OUT, "modeles_comptage.rds"))
annuaire <- read_parquet(file.path(DOSSIER, "fr-en-annuaire-education.parquet"))
d_mod <- charg$d_mod; modeles <- charg$modeles; m_p <- modeles$Poisson

codes_ech <- unique(master$code_insee)

colleges_sf <- annuaire %>%
  filter(grepl("^Coll", type_etablissement),
         toupper(coalesce(as.character(etat), "")) == "OUVERT",
         !is.na(latitude), !is.na(longitude)) %>%
  mutate(code_insee = str_pad(as.character(code_commune), 5, pad = "0"),
         secteur = case_when(statut_public_prive == "Public" ~ "public",
                             grepl("^Priv", statut_public_prive) ~ "prive",
                             TRUE ~ NA_character_)) %>%
  filter(!is.na(secteur), code_insee %in% codes_ech) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(4326)

# Centroïdes des communes en WGS84, en data.frame simple (lng/lat explicites)
master_geo <- master %>%
  filter(!is.na(x_pop), !is.na(y_pop)) %>%
  mutate(lng = x_pop, lat = y_pop) %>%
  st_as_sf(coords = c("x_pop", "y_pop"), crs = 2154) %>%
  st_transform(4326) %>%
  mutate(lng = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2]) %>%
  st_drop_geometry() %>%
  filter(is.finite(lng), is.finite(lat))

# ---- IRR / AIC précalculés ----
V_cl <- tryCatch(vcovCL(m_p, cluster = d_mod$region_acad, type = "HC1"),
                 error = function(e) vcov(m_p))
ct   <- coeftest(m_p, vcov. = V_cl)
irr_tbl <- data.frame(variable = rownames(ct), coef = ct[, 1], se = ct[, 2],
                      IRR = exp(ct[, 1]),
                      IC_bas = exp(ct[, 1] - 1.96 * ct[, 2]),
                      IC_haut = exp(ct[, 1] + 1.96 * ct[, 2]),
                      p = ct[, 4], row.names = NULL)
aic_tbl <- data.frame(modele = names(modeles),
                      AIC = sapply(modeles, AIC), BIC = sapply(modeles, BIC))

# ---- Cache des carreaux par ville (test key_missing, pas is.null) ----
cache_carreaux <- cachem::cache_mem(max_size = 500e6)
lire_carreaux <- function(code) {
  key <- paste0("car_", code)
  hit <- cache_carreaux$get(key)
  if (!cachem::is.key_missing(hit)) return(hit)
  car <- st_read(CHEMIN_CARREAUX,
                 query = sprintf("SELECT * FROM carreaux_200m_met WHERE lcog_geo LIKE '%s%%'", code),
                 quiet = TRUE) %>%
    filter(is.finite(ind), ind > 0, is.finite(ind_snv)) %>%
    mutate(revenu = ind_snv / ind) %>%
    select(revenu, ind, geometry) %>%
    st_make_valid() %>%
    st_transform(4326)
  cache_carreaux$set(key, car)
  car
}

pal_reg <- colorFactor("turbo", domain = levels(master$region_acad))

# ---- Thème sombre ----
css <- "
body { background:#0f172a; color:#e2e8f0; font-family:'Inter',system-ui,sans-serif; }
.nav-tabs .nav-link { color:#94a3b8; border:none; }
.nav-tabs .nav-link.active { background:#1e293b !important; color:#f8fafc !important; border-bottom:2px solid #38bdf8 !important; }
.card { background:#1e293b !important; border:1px solid #334155 !important; border-radius:12px !important; }
.card-header { background:#1e293b !important; color:#f8fafc !important; border-bottom:1px solid #334155 !important; font-weight:600; }
.bslib-value-box { background:#1e293b !important; border:1px solid #334155 !important; border-radius:12px !important; }
.bslib-value-box .value-box-value { color:#38bdf8 !important; font-size:1.8rem !important; }
.bslib-value-box .value-box-title { color:#94a3b8 !important; }
.sidebar { background:#1e293b !important; border-right:1px solid #334155 !important; }
.form-control, select, .selectize-input { background:#0f172a !important; color:#e2e8f0 !important; border:1px solid #334155 !important; }
.irs-bar, .irs-bar-edge, .irs-single, .irs-from, .irs-to { background:#38bdf8 !important; border-color:#38bdf8 !important; }
.leaflet-container { border-radius:12px; }
.dataTables_wrapper, table.dataTable, table.dataTable td { color:#e2e8f0; }
table.dataTable thead th { background:#334155 !important; color:#f8fafc !important; }
table.dataTable tbody tr { background:#1e293b !important; }
.dataTables_paginate .paginate_button { color:#94a3b8 !important; }
.dataTables_info, .dataTables_length, .dataTables_filter { color:#94a3b8 !important; }
h3, h4 { color:#f8fafc; }
.help-block, .form-text { color:#64748b; }
"

# ============================ UI ============================
ui <- page_sidebar(
  title = div(icon("graduation-cap"), " Collèges privés & ségrégation urbaine"),
  theme = bs_theme(version = 5, bg = "#0f172a", fg = "#e2e8f0",
                   primary = "#38bdf8", base_font = font_google("Inter")),
  tags$head(tags$style(HTML(css))),
  sidebar = sidebar(
    width = 290,
    selectInput("region", "Région académique",
                choices = c("Toutes les régions" = "Toutes", levels(master$region_acad))),
    sliderInput("popmin", "Population min.",
                min = 10000, max = max(master$pop_totale, na.rm = TRUE),
                value = 100000, step = 10000),
    selectizeInput("ville", "Ville", choices = NULL, width = "100%",
                   options = list(placeholder = "Choisir une ville…")),
    hr(),
    div(class = "help-block",
        "Sources : Filosofi 200 m, annuaire Éducation nationale, modèles Code 4b (M4). Localhost uniquement.")
  ),
  navset_card_tab(
    nav_panel("Vue d'ensemble", icon("earth-americas"),
              layout_columns(
                value_box("Communes", format(nrow(master), big.mark = " "), showcase = bs_icon("map"), theme = "primary"),
                value_box("Collèges privés", format(sum(master$nb_colleges_prives, na.rm = TRUE), big.mark = " "), showcase = bs_icon("mortarboard"), theme = "danger"),
                value_box("Sans collège privé", paste0(round(100 * mean(master$nb_colleges_prives == 0, na.rm = TRUE), 1), " %"), showcase = bs_icon("slash-circle"), theme = "warning"),
                value_box("Régions académiques", nlevels(master$region_acad), showcase = bs_icon("flag"), theme = "info")
              ),
              card(full_screen = TRUE, leafletOutput("carte_france", height = "68vh"))
    ),
    nav_panel("Résultats globaux", icon("chart-line"),
              layout_columns(
                card(full_screen = TRUE, plotlyOutput("scatter_ind", height = "46vh")),
                card(full_screen = TRUE, plotOutput("heatmap", height = "46vh"))
              ),
              layout_columns(
                card(full_screen = TRUE, plotlyOutput("pred_obs", height = "42vh")),
                card(full_screen = TRUE, plotOutput("resid_region", height = "42vh"))
              ),
              layout_columns(
                card(full_screen = TRUE, DTOutput("tbl_irr")),
                card(full_screen = TRUE, DTOutput("tbl_aic"))
              )
    ),
    nav_panel("Ville sélectionnée", icon("city"),
              uiOutput("ville_header"),
              layout_columns(
                card(full_screen = TRUE, leafletOutput("carte_ville", height = "60vh")),
                card(full_screen = TRUE, plotOutput("barres_indices", height = "60vh"))
              ),
              layout_columns(
                card(full_screen = TRUE, DTOutput("tbl_colleges_ville")),
                card(full_screen = TRUE, DTOutput("fiche_ville"))
              )
    ),
    nav_panel("Données", icon("table"),
              card(full_screen = TRUE, DTOutput("tbl_master"))
    )
  )
)

# ============================ SERVER ============================
server <- function(input, output, session) {
  
  villes_dispo <- reactive({
    df <- master
    if (input$region != "Toutes") df <- df %>% filter(as.character(region_acad) == input$region)
    df %>%
      filter(!is.na(pop_totale), pop_totale >= input$popmin) %>%
      arrange(desc(pop_totale)) %>%
      mutate(label = paste0(nom_commune, " — ", format(pop_totale, big.mark = " "),
                            " hab., ", nb_colleges_prives, " privé(s)"))
  })
  
  observe({
    v <- villes_dispo()
    updateSelectizeInput(session, "ville", choices = setNames(v$code_insee, v$label), server = TRUE)
  })
  
  ville_sel <- reactive({ req(input$ville); master %>% filter(code_insee == input$ville) })
  
  # ---- Carte de France : lng/lat explicites, NA filtrés ----
  output$carte_france <- renderLeaflet({
    mg  <- master_geo
    pal <- colorNumeric("YlOrRd", domain = mg$nb_colleges_prives)
    cc  <- st_coordinates(colleges_sf)
    keep <- is.finite(cc[, 1]) & is.finite(cc[, 2])
    cc  <- cc[keep, , drop = FALSE]
    sect <- colleges_sf$secteur[keep]
    leaflet() %>%
      addProviderTiles(providers$CartoDB.DarkMatter) %>%
      addCircleMarkers(lng = mg$lng, lat = mg$lat,
                       radius = ~sqrt(pmax(nb_colleges_prives, 0)) * 2 + 2,
                       fillColor = ~pal_reg(as.character(region_acad)), fillOpacity = 0.8,
                       color = "#94a3b8", weight = 0.5, group = "Communes",
                       label = ~paste0(nom_commune, " : ", nb_colleges_prives, " privé(s)")) %>%
      addCircleMarkers(lng = cc[, 1], lat = cc[, 2],
                       color = ifelse(sect == "prive", "#f87171", "#60a5fa"),
                       radius = 3, fillOpacity = 0.85, stroke = FALSE, group = "Collèges",
                       label = ~paste0(sect, " — ", colleges_sf$code_insee[keep])) %>%
      addLayersControl(overlayGroups = c("Communes", "Collèges"),
                       options = layersControlOptions(collapsed = FALSE)) %>%
      addLegend("bottomright", pal = pal_reg, values = levels(master$region_acad),
                title = "Région académique", opacity = 0.9) %>%
      addLegend("bottomleft", pal = pal, values = mg$nb_colleges_prives,
                title = "Nb privés (commune)", opacity = 0.9)
  })
  
  # ---- Scatter indices (caché) ----
  output$scatter_ind <- renderPlotly({
    viz <- master %>%
      select(nom_commune, nb_colleges_prives, region_acad,
             z_theil, z_moran, z_exposition_riches, z_polarisation,
             z_part_pop_top10, z_getis_ord_riches, z_segregation_multidim) %>%
      pivot_longer(cols = starts_with("z_"), names_to = "indice", values_to = "z") %>%
      filter(!is.na(z), !is.na(nb_colleges_prives))
    p <- ggplot(viz, aes(z, nb_colleges_prives, color = as.character(region_acad), text = nom_commune)) +
      geom_point(alpha = 0.45) +
      geom_smooth(method = "lm", se = FALSE, color = "#e2e8f0", linewidth = 0.5) +
      facet_wrap(~indice, scales = "free_x", ncol = 3) +
      scale_y_sqrt() + scale_color_viridis_d(option = "turbo") +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "#1e293b", color = NA),
            plot.background = element_rect(fill = "#1e293b", color = NA),
            text = element_text(color = "#e2e8f0"),
            axis.text = element_text(color = "#94a3b8"),
            strip.text = element_text(color = "#f8fafc")) +
      labs(x = "Indice (z-score)", y = "Nb collèges privés (racine)", color = "Région")
    ggplotly(p, tooltip = c("text", "x", "y")) %>%
      layout(paper_bgcolor = "#1e293b", plot_bgcolor = "#1e293b", font = list(color = "#e2e8f0"))
  }) %>% bindCache(master)
  
  output$heatmap <- renderPlot({
    inds <- c("theil_revenu", "moran_I_rev", "exposition_riches", "polarisation",
              "part_pop_top10", "getis_ord_riches", "segregation_multidim",
              "taux_pauv_commune", "revenu_moyen_com")
    cm <- cor(master[, inds], use = "pairwise.complete.obs")
    df <- as.data.frame(cm) %>% rownames_to_column("v1") %>%
      pivot_longer(-v1, names_to = "v2", values_to = "r")
    ggplot(df, aes(v2, v1, fill = r)) +
      geom_tile() + geom_text(aes(label = round(r, 2)), size = 3.2, color = "white") +
      scale_fill_gradient2(low = "#f87171", mid = "#1e293b", high = "#60a5fa",
                           midpoint = 0, limits = c(-1, 1)) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "#94a3b8"),
            axis.text.y = element_text(color = "#94a3b8"),
            panel.background = element_rect(fill = "#1e293b", color = NA),
            plot.background = element_rect(fill = "#1e293b", color = NA),
            text = element_text(color = "#e2e8f0"),
            plot.title = element_text(color = "#f8fafc")) +
      labs(x = NULL, y = NULL, title = "Corrélations entre indicateurs")
  }) %>% bindCache(master)
  
  output$pred_obs <- renderPlotly({
    pv <- master %>% filter(!is.na(y_pred), !is.na(nb_colleges_prives))
    p <- ggplot(pv, aes(y_pred, nb_colleges_prives, color = as.character(region_acad), text = nom_commune)) +
      geom_point(alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#64748b") +
      scale_x_sqrt() + scale_y_sqrt() + scale_color_viridis_d(option = "turbo") +
      theme_minimal(base_size = 10) +
      theme(panel.background = element_rect(fill = "#1e293b", color = NA),
            plot.background = element_rect(fill = "#1e293b", color = NA),
            text = element_text(color = "#e2e8f0"),
            axis.text = element_text(color = "#94a3b8"),
            plot.title = element_text(color = "#f8fafc")) +
      labs(x = "Prédit (M4, racine)", y = "Observé (racine)", color = "Région",
           title = "Modèle M4 : prédit vs observé")
    ggplotly(p, tooltip = c("text", "x", "y")) %>%
      layout(paper_bgcolor = "#1e293b", plot_bgcolor = "#1e293b", font = list(color = "#e2e8f0"))
  }) %>% bindCache(master)
  
  output$resid_region <- renderPlot({
    master %>% filter(!is.na(resid)) %>%
      ggplot(aes(region_acad, resid, fill = region_acad)) +
      geom_boxplot(outlier.alpha = 0.3, color = "#94a3b8") + coord_flip() +
      scale_fill_viridis_d(option = "turbo") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none",
            panel.background = element_rect(fill = "#1e293b", color = NA),
            plot.background = element_rect(fill = "#1e293b", color = NA),
            text = element_text(color = "#e2e8f0"),
            axis.text = element_text(color = "#94a3b8"),
            plot.title = element_text(color = "#f8fafc")) +
      labs(x = NULL, y = "Résidu de Pearson (M4)", title = "Résidus par région académique")
  }) %>% bindCache(master)
  
  dt_opts <- list(pageLength = 15, scrollY = "380px", dom = "Bfrtip",
                  buttons = c("copy", "csv", "excel"),
                  language = list(url = "//cdn.datatables.net/plug-ins/1.13.6/i18n/fr-FR.json"))
  
  output$tbl_irr <- renderDT({
    irr_tbl %>% mutate(across(where(is.numeric), ~round(., 4))) %>%
      datatable(rownames = FALSE, options = dt_opts, style = "bootstrap",
                caption = "Modèle M4 — IRR et IC 95 % (cluster région académique)")
  })
  output$tbl_aic <- renderDT({
    aic_tbl %>% mutate(across(c(AIC, BIC), ~round(., 1))) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10), style = "bootstrap",
                caption = "Comparaison des modèles (AIC / BIC)")
  })
  
  # ---- Ville sélectionnée ----
  output$ville_header <- renderUI({
    v <- ville_sel(); req(nrow(v) > 0)
    h3(icon("city"), paste0(v$nom_commune, " (", v$code_insee, ") — ", as.character(v$region_acad)))
  })
  
  output$carte_ville <- renderLeaflet({
    v <- ville_sel(); req(nrow(v) > 0)
    car  <- lire_carreaux(v$code_insee)   # en cache après le 1er appel
    cols <- colleges_sf %>% filter(code_insee == v$code_insee)
    pal  <- colorNumeric("magma", domain = car$revenu, na.color = "#334155")
    m <- leaflet() %>%
      addProviderTiles(providers$CartoDB.DarkMatter) %>%
      addPolygons(data = car, fillColor = ~pal(revenu), fillOpacity = 0.75,
                  color = "#475569", weight = 0.3, group = "Carreaux 200 m",
                  label = ~paste0(round(revenu), " €/uc — pop ", ind))
    if (nrow(cols) > 0) {
      cc <- st_coordinates(cols)
      keep <- is.finite(cc[, 1]) & is.finite(cc[, 2])
      m <- m %>% addCircleMarkers(lng = cc[keep, 1], lat = cc[keep, 2],
                                  fillColor = ifelse(cols$secteur[keep] == "prive", "#f87171", "#60a5fa"),
                                  color = "#f8fafc", weight = 1.5, radius = 8,
                                  fillOpacity = 0.9, group = "Collèges",
                                  label = paste0(cols$secteur[keep], " — UAI ",
                                                 cols$identifiant_de_l_etablissement[keep]))
    }
    m %>%
      addLayersControl(overlayGroups = c("Carreaux 200 m", "Collèges"),
                       options = layersControlOptions(collapsed = FALSE)) %>%
      addLegend("bottomright", pal = pal, values = car$revenu, title = "Revenu (€/uc)", opacity = 0.9)
  })
  
  output$barres_indices <- renderPlot({
    v <- ville_sel(); req(nrow(v) > 0)
    inds <- c("z_theil", "z_moran", "z_exposition_riches", "z_polarisation",
              "z_part_pop_top10", "z_getis_ord_riches", "z_segregation_multidim")
    df <- tibble(indice = inds, z = as.numeric(v[1, inds]))
    ggplot(df, aes(reorder(indice, z), z)) +
      geom_col(fill = "#38bdf8", alpha = 0.9) + coord_flip() +
      geom_hline(yintercept = 0, color = "#64748b") +
      geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "#475569") +
      theme_minimal(base_size = 11) +
      theme(panel.background = element_rect(fill = "#1e293b", color = NA),
            plot.background = element_rect(fill = "#1e293b", color = NA),
            text = element_text(color = "#e2e8f0"),
            axis.text = element_text(color = "#94a3b8"),
            plot.title = element_text(color = "#f8fafc")) +
      labs(x = NULL, y = "z-score", title = paste("Indices de ségrégation —", v$nom_commune))
  })
  
  output$tbl_colleges_ville <- renderDT({
    v <- ville_sel(); req(nrow(v) > 0)
    cols <- colleges_sf %>% st_drop_geometry() %>% filter(code_insee == v$code_insee)
    nom_col <- intersect(c("appellation_officielle", "patronyme", "nom_etablissement",
                           "denomination_principale"), names(cols))[1]
    cols %>% select(code_insee, secteur, type_contrat_prive, all_of(nom_col)) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10), style = "bootstrap",
                caption = paste("Collèges à", v$nom_commune, "(", nrow(cols), ")"))
  })
  
  output$fiche_ville <- renderDT({
    v <- ville_sel(); req(nrow(v) > 0)
    data.frame(
      indicateur = c("Population", "Revenu moyen (€/uc)", "Taux de pauvreté (%)",
                     "Nb collèges publics", "Nb collèges privés", "  dont sous contrat",
                     "  dont hors contrat", "Dist. collège public + proche (km)",
                     "Privés à 5 km (hors commune)", "Theil", "Moran I",
                     "Exposition riches (Pxx)", "Polarisation (W)", "Part pop top 10 %",
                     "Ségrégation multidim.", "Prédit M4 (attendu)", "Observé", "Résidu de Pearson"),
      valeur = c(format(v$pop_totale, big.mark = " "), round(v$revenu_moyen_com),
                 round(v$taux_pauv_commune, 1), v$nb_colleges_publics, v$nb_colleges_prives,
                 v$nb_prives_sous_contrat, v$nb_prives_hors_contrat,
                 round(v$dist_college_public_km, 2), v$nb_prives_5km_autres,
                 round(v$theil_revenu, 4), round(v$moran_I_rev, 3),
                 round(v$exposition_riches, 3), round(v$polarisation, 3),
                 round(v$part_pop_top10, 3), round(v$segregation_multidim, 3),
                 round(v$y_pred, 2), v$nb_colleges_prives, round(v$resid, 2))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20),
                    style = "bootstrap", caption = paste("Fiche —", v$nom_commune))
  })
  
  # ---- Données brutes ----
  output$tbl_master <- renderDT({
    master %>%
      select(code_insee, nom_commune, region_acad, pop_totale, revenu_moyen_com,
             taux_pauv_commune, nb_colleges_publics, nb_colleges_prives, y_pred, resid,
             theil_revenu, moran_I_rev, exposition_riches, polarisation,
             part_pop_top10, segregation_multidim,
             z_theil, z_moran, z_exposition_riches, z_polarisation) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      datatable(filter = "top", rownames = FALSE, style = "bootstrap",
                options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
                               buttons = c("copy", "csv", "excel"),
                               language = list(url = "//cdn.datatables.net/plug-ins/1.13.6/i18n/fr-FR.json")),
                caption = paste("Master —", nrow(master), "communes"))
  })
}

shinyApp(ui, server)