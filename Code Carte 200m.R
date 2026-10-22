library(leaflet)
library(sf)
library(dplyr)
library(htmltools)
library(data.table)
library(stringr)
library(jsonlite)

# 1. Dossier de ton projet local et dossier des données brutes
DOSSIER <- "C:/Users/camil/Documents/DepartementsPauvreteLycees/"
dossier_donnees <- "G:/Mon Drive/Education/Master APE 25-27/Mémoire/Données MicroMetrics/france_metro/"

shp_path  <- paste0(dossier_donnees, "carreaux_200m_met.shp")
gpkg_path <- paste0(dossier_donnees, "carreaux_200m_met.gpkg")

# 2. Conversion unique en GeoPackage si non présent (optimisation vitesse)
if (!file.exists(gpkg_path)) {
  message("⏳ Conversion du Shapefile 200m en GeoPackage optimisé (patientez 30-40 secondes)...")
  carreaux_200m_brut <- st_read(shp_path, quiet = TRUE)
  st_write(carreaux_200m_brut, gpkg_path, delete_dsn = TRUE, quiet = TRUE)
  rm(carreaux_200m_brut)
  gc()
  message("✅ Conversion terminée avec succès !")
}

# 3. Récupérer la liste des départements
deps_df <- fromJSON("https://geo.api.gouv.fr/departements") %>%
  select(code, nom) %>%
  arrange(code)

centroides_l93_trans <- st_transform(centroides_l93, 2154)
dep_files <- data.frame(code = character(), nom = character(), filename = character(), stringsAsFactors = FALSE)

# 4. Boucle de génération : un fichier HTML par département avec les carreaux 200m et salaires/revenus
# 3. Boucle de génération : un fichier HTML par département avec TOUTES les variables
for(i in 1:nrow(deps_df)) {
  code_dept <- deps_df$code[i]
  nom_dept <- deps_df$nom[i]
  
  message(sprintf("Génération du département %s - %s (%d/101)...", code_dept, nom_dept, i))
  
  query_sql <- sprintf("SELECT * FROM carreaux_200m_met WHERE lcog_geo LIKE '%s%%'", code_dept)
  carreaux_sf <- tryCatch({
    st_read(gpkg_path, query = query_sql, quiet = TRUE)
  }, error = function(e) NULL)
  
  if(is.null(carreaux_sf) || nrow(carreaux_sf) == 0) next 
  
  carreaux_sf <- carreaux_sf %>%
    st_transform(4326) %>%
    mutate(
      taux_pauvrete = ifelse(men > 0, (men_pauv / men) * 100, 0),
      part_proprio  = ifelse(men > 0, (men_prop / men) * 100, 0),
      revenu_moyen  = ifelse(ind > 0, ind_snv / ind, 0),
      part_maison   = ifelse(men > 0, (men_mais / men) * 100, 0),
      part_collectif= ifelse(men > 0, (men_coll / men) * 100, 0)
    )
  
  # Palettes de couleurs par thématique
  pal_pop   <- colorNumeric("Blues", domain = carreaux_sf$ind, na.color = "transparent")
  pal_pauv  <- colorNumeric(c("#22c55e", "#fef08a", "#ef4444"), domain = carreaux_sf$taux_pauvrete, na.color = "transparent")
  pal_prop  <- colorNumeric("YlGn", domain = carreaux_sf$part_proprio, na.color = "transparent")
  pal_rev   <- colorNumeric("Purples", domain = carreaux_sf$revenu_moyen, na.color = "transparent")
  pal_jeune <- colorNumeric("Oranges", domain = carreaux_sf$ind_0_3 + carreaux_sf$ind_4_5 + carreaux_sf$ind_6_10, na.color = "transparent")
  pal_aines <- colorNumeric("Reds", domain = carreaux_sf$ind_65_79 + carreaux_sf$ind_80p, na.color = "transparent")
  pal_maison<- colorNumeric("Greens", domain = carreaux_sf$part_maison, na.color = "transparent")
  
  url_contours <- paste0("https://geo.api.gouv.fr/departements/", code_dept, "/communes?format=geojson&geometry=contour")
  communes_wgs84 <- tryCatch(st_read(url_contours, quiet = TRUE), error = function(e) NULL)
  if(is.null(communes_wgs84) || nrow(communes_wgs84) == 0) next
  
  dept_centroides_wgs84 <- centroides_l93_trans %>%
    filter(INSEE_DEP == code_dept) %>%
    st_transform(4326)
  
  dept_lycees_publics_wgs84 <- lycees_clean %>%
    filter(statut_public_prive == "Public", type_etablissement == "Lycée", 
           code_departement %in% c(code_dept, paste0("0", code_dept), str_remove(code_dept, "^0"))) %>%
    st_transform(4326)
  
  dept_lycees_prives_wgs84 <- lycees_clean %>%
    filter(statut_public_prive == "Privé", type_etablissement == "Lycée", 
           code_departement %in% c(code_dept, paste0("0", code_dept), str_remove(code_dept, "^0"))) %>%
    st_transform(4326)
  
  # Carte Leaflet avec toutes les variables organisées
  m <- leaflet(width = "100%", height = "750px", options = leafletOptions(zoomControl = TRUE)) %>%
    addProviderTiles(providers$Esri.WorldGrayCanvas, group = "Fond de carte") %>%
    
    # 1. Démographie & Population
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.2,
      fillColor = ~pal_pop(ind), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau 200m :</b> ", idcar_200m, "<br>",
                      "<b>👥 Population totale :</b> ", ind, "<br>",
                      "<b>🏠 Ménages :</b> ", men, "<br>",
                      "<b>👶 0-10 ans :</b> ", ind_0_3 + ind_4_5 + ind_6_10, "<br>",
                      "<b>🧑 11-24 ans :</b> ", ind_11_17 + ind_18_24, "<br>",
                      "<b>👴 65 ans et + :</b> ", ind_65_79 + ind_80p),
      group = "1. Population totale"
    ) %>%
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.2,
      fillColor = ~pal_jeune(ind_0_3 + ind_4_5 + ind_6_10), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau 200m :</b> ", idcar_200m, "<br><b>👶 Enfants (0-10 ans) :</b> ", ind_0_3 + ind_4_5 + ind_6_10),
      group = "1b. Population jeune (0-10 ans)"
    ) %>%
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.2,
      fillColor = ~pal_aines(ind_65_79 + ind_80p), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau 200m :</b> ", idcar_200m, "<br><b>👴 Seniors (65 ans et +) :</b> ", ind_65_79 + ind_80p),
      group = "1c. Population senior (65+)"
    ) %>%
    
    # 2. Revenus & Pauvreté
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.2,
      fillColor = ~pal_pauv(taux_pauvrete), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau 200m :</b> ", idcar_200m, "<br>",
                      "<b>📉 Taux de pauvreté :</b> ", round(taux_pauvrete, 1), "%<br>",
                      "<b>Ménages pauvres :</b> ", men_pauv, " / ", men),
      group = "2. Taux de pauvreté"
    ) %>%
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.2,
      fillColor = ~pal_rev(revenu_moyen), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau 200m :</b> ", idcar_200m, "<br><b>💶 Salaire / Revenu moyen :</b> ", round(revenu_moyen, 0), " €"),
      group = "2b. Salaire & Revenu moyen"
    ) %>%
    
    # 3. Logement & Propriété
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.2,
      fillColor = ~pal_prop(part_proprio), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau 200m :</b> ", idcar_200m, "<br>",
                      "<b>🏢 Part propriétaires :</b> ", round(part_proprio, 1), "%<br>",
                      "<b>Ménages propriétaires :</b> ", men_prop, " / ", men),
      group = "3. Part de propriétaires"
    ) %>%
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.2,
      fillColor = ~pal_maison(part_maison), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau 200m :</b> ", idcar_200m, "<br>",
                      "<b>🏡 Part maisons :</b> ", round(part_maison, 1), "%<br>",
                      "<b>🏢 Part logements collectifs :</b> ", round(part_collectif, 1), "%<br>",
                      "<b>Surface moyenne :</b> ", round(men_surf, 1), " m²"),
      group = "3b. Type d'habitat (Maisons vs Collectif)"
    ) %>%
    
    # 4. Éléments géographiques et scolaires de superposition
    addPolygons(
      data = communes_wgs84, color = "#1e293b", weight = 1, fillColor = "transparent",
      label = ~paste0(nom, " (", code, ")"),
      group = "Communes"
    ) %>%
    addCircleMarkers(
      data = dept_centroides_wgs84, radius = 3, color = "#f43f5e", stroke = FALSE, fillOpacity = 0.8,
      popup = ~paste0("<b>Centroïde :</b> ", NOM, "<br><b>Dist. Lycée Public :</b> ", round(dist_lycee_public_km, 2), " km"),
      group = "Centroïdes population"
    ) %>%
    addCircleMarkers(
      data = dept_lycees_publics_wgs84, radius = 6, color = "black", fillColor = "black", stroke = TRUE, weight = 1, fillOpacity = 0.9,
      popup = ~paste0("<b>🏫 Lycée Public :</b> ", nom_etablissement, "<br><b>Commune :</b> ", nom_commune),
      group = "Lycées Publics"
    ) %>%
    addCircleMarkers(
      data = dept_lycees_prives_wgs84, radius = 6, color = "#475569", fillColor = "#888888", stroke = TRUE, weight = 1, fillOpacity = 0.9,
      popup = ~paste0("<b>🏫 Lycée Privé :</b> ", nom_etablissement, "<br><b>Commune :</b> ", nom_commune, "<br><b>Contrat :</b> ", type_contrat_prive),
      group = "Lycées Privés"
    ) %>%
    addLayersControl(
      baseGroups = c(
        "1. Population totale", 
        "1b. Population jeune (0-10 ans)", 
        "1c. Population senior (65+)", 
        "2. Taux de pauvreté", 
        "2b. Salaire & Revenu moyen", 
        "3. Part de propriétaires", 
        "3b. Type d'habitat (Maisons vs Collectif)"
      ),
      overlayGroups = c("Communes", "Centroïdes population", "Lycées Publics", "Lycées Privés"),
      options = layersControlOptions(collapsed = FALSE)
    )
  
  # Structure HTML
  dept_html <- tags$html(
    tags$head(
      tags$style(HTML("
        body { font-family: 'Segoe UI', sans-serif; background: #f8fafc; margin: 0; padding: 20px; }
        .container { max-width: 1400px; margin: 0 auto; background: white; padding: 25px; border-radius: 12px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); }
        .header-bar { display: flex; align-items: center; justify-content: space-between; margin-bottom: 20px; }
        h1 { font-size: 22px; color: #0f172a; margin: 0; }
        .btn-back { background: #38bdf8; color: #0f172a; padding: 8px 16px; border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 13px; transition: background 0.2s; }
        .btn-back:hover { background: #0ea5e9; }
      "))
    ),
    tags$body(
      tags$div(class = "container",
               tags$div(class = "header-bar",
                        tags$h1(paste("Département :", nom_dept, "(", code_dept, ") - Atlas 200m Intégral")),
                        tags$a(href = "index.html", class = "btn-back", "← Retour au menu national")
               ),
               m
      )
    )
  )
  
  filename <- paste0("carte_dept_", code_dept, ".html")
  filepath <- paste0(DOSSIER, filename)
  save_html(dept_html, file = filepath)
  
  dep_files <- rbind(dep_files, data.frame(code = code_dept, nom = nom_dept, filename = filename, stringsAsFactors = FALSE))
}

# 4. Page d'accueil index.html
list_items <- apply(dep_files, 1, function(row) {
  tags$li(
    tags$a(href = row["filename"], 
           tags$span(class = "code-badge", row["code"]),
           tags$span(class = "dep-name", row["nom"])
    )
  )
})

index_html <- tags$html(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background: #f1f5f9; margin: 0; padding: 40px; color: #1e293b; }
      .container { max-width: 950px; margin: 0 auto; background: white; padding: 40px; border-radius: 12px; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.05); }
      h1 { margin-top: 0; color: #0f172a; font-size: 26px; font-weight: 600; }
      p.subtitle { color: #64748b; font-size: 15px; margin-bottom: 25px; }
      .search-box { width: 100%; padding: 14px 20px; font-size: 16px; border: 1px solid #cbd5e1; border-radius: 8px; margin-bottom: 30px; box-sizing: border-box; outline: none; }
      .search-box:focus { border-color: #38bdf8; box-shadow: 0 0 0 3px rgba(56, 189, 248, 0.2); }
      ul.dep-grid { list-style: none; padding: 0; display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 12px; max-height: 60vh; overflow-y: auto; }
      ul.dep-grid li a { display: flex; align-items: center; background: #f8fafc; border: 1px solid #e2e8f0; padding: 14px 18px; border-radius: 8px; text-decoration: none; color: #334155; font-weight: 500; transition: all 0.2s; }
      ul.dep-grid li a:hover { background: #f0fdf4; border-color: #22c55e; transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.03); }
      .code-badge { background: #0f172a; color: white; padding: 4px 8px; border-radius: 4px; font-size: 12px; font-weight: 600; margin-right: 12px; }
      .dep-name { font-size: 15px; }
    "))
  ),
  tags$body(
    tags$div(class = "container",
             tags$h1("Atlas National des Lycées & Données INSEE (Maillage 200m Intégral)"),
             tags$p(class = "subtitle", "Toutes les variables socio-démographiques, de pauvreté, de salaires et de logements sont intégrées et basculables sur les cartes départementales."),
             tags$input(type = "text", id = "search", class = "search-box", placeholder = "🔍 Filtrer par nom ou numéro de département (ex: 35, Paris, Rhône)...", onkeyup = "filterDeps()"),
             tags$ul(id = "depList", class = "dep-grid", list_items),
             tags$script(HTML("
        function filterDeps() {
          let input = document.getElementById('search').value.toLowerCase();
          let items = document.getElementById('depList').getElementsByTagName('li');
          for (let i = 0; i < items.length; i++) {
            let text = items[i].textContent || items[i].innerText;
            items[i].style.display = text.toLowerCase().indexOf(input) > -1 ? '' : 'none';
          }
        }
      "))
    )
  )
)

index_path <- paste0(DOSSIER, "index.html")
save_html(index_html, file = index_path)

message("✅ Génération complète de l'atlas intégral avec toutes les variables réussie !")