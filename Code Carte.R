library(leaflet)
library(sf)
library(dplyr)
library(htmltools)
library(data.table)
library(stringr)
library(jsonlite)

# 1. Récupérer la liste des départements
deps_df <- fromJSON("https://geo.api.gouv.fr/departements") %>%
  select(code, nom) %>%
  arrange(code)

# 2. EXTRACTION GLOBALE DES COORDONNÉES DES CARREAUX (Correction ici)
message("Extraction globale des coordonnées des carreaux INSEE...")
donnees_carreaux$i_est_1km  <- as.numeric(str_extract(donnees_carreaux$idcar_1km, "(?<=E)\\d+"))
donnees_carreaux$i_nord_1km <- as.numeric(str_extract(donnees_carreaux$idcar_1km, "(?<=N)\\d+(?=E)"))

# 3. Préparation des centroïdes
centroides_l93_trans <- st_transform(centroides_l93, 2154)

# Dataframe pour stocker les liens du menu index
dep_files <- data.frame(code = character(), nom = character(), filename = character(), stringsAsFactors = FALSE)

# 4. Boucle de génération : un fichier HTML par département
for(i in 1:nrow(deps_df)) {
  code_dept <- deps_df$code[i]
  nom_dept <- deps_df$nom[i]
  
  message(sprintf("Génération du département %s - %s (%d/%d)...", code_dept, nom_dept, i, nrow(deps_df)))
  
  # Filtrer les carreaux du département
  carreaux_dept <- donnees_carreaux %>%
    filter(substr(as.character(lcog_geo), 1, 2) == code_dept) %>%
    filter(!is.na(i_est_1km) & !is.na(i_nord_1km))
  
  if(nrow(carreaux_dept) == 0) next 
  
  carreaux_sf <- carreaux_dept %>%
    st_as_sf(coords = c("i_est_1km", "i_nord_1km"), crs = 3035, remove = FALSE) %>%
    st_buffer(dist = 500, endCapStyle = "SQUARE") %>%
    st_transform(4326) %>%
    mutate(
      taux_pauvrete = ifelse(men > 0, (men_pauv / men) * 100, 0),
      part_proprio  = ifelse(men > 0, (men_prop / men) * 100, 0)
    )
  
  pal_pop  <- colorNumeric("Blues", domain = carreaux_sf$ind, na.color = "transparent")
  pal_pauv <- colorNumeric(c("#22c55e", "#fef08a", "#ef4444"), domain = carreaux_sf$taux_pauvrete, na.color = "transparent")
  pal_prop <- colorNumeric("YlGn", domain = carreaux_sf$part_proprio, na.color = "transparent")
  
  # Récupérer les contours des communes du département
  url_contours <- paste0("https://geo.api.gouv.fr/departements/", code_dept, "/communes?format=geojson&geometry=contour")
  communes_wgs84 <- tryCatch(st_read(url_contours, quiet = TRUE), error = function(e) NULL)
  if(is.null(communes_wgs84) || nrow(communes_wgs84) == 0) next
  
  # Filtrer centroïdes et lycées
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
  
  # Construction de la carte Leaflet
  m <- leaflet(width = "100%", height = "100%", options = leafletOptions(zoomControl = TRUE)) %>%
    addProviderTiles(providers$Esri.WorldGrayCanvas, group = "Fond de carte") %>%
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.3,
      fillColor = ~pal_pop(ind), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau :</b> ", idcar_1km, "<br><b>👥 Population :</b> ", ind),
      group = "1. Population (ind)"
    ) %>%
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.3,
      fillColor = ~pal_pauv(taux_pauvrete), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau :</b> ", idcar_1km, "<br><b>📉 Taux pauvreté :</b> ", round(taux_pauvrete, 1), "%"),
      group = "2. Taux de pauvreté (Vert -> Rouge)"
    ) %>%
    addPolygons(
      data = carreaux_sf, color = "#cbd5e1", weight = 0.3,
      fillColor = ~pal_prop(part_proprio), fillOpacity = 0.65,
      popup = ~paste0("<b>Carreau :</b> ", idcar_1km, "<br><b>🏢 Part propriétaires :</b> ", round(part_proprio, 1), "%"),
      group = "3. Part de propriétaires"
    ) %>%
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
      baseGroups = c("1. Population (ind)", "2. Taux de pauvreté (Vert -> Rouge)", "3. Part de propriétaires"),
      overlayGroups = c("Communes", "Centroïdes population", "Lycées Publics", "Lycées Privés"),
      options = layersControlOptions(collapsed = FALSE)
    )
  
  # Structure HTML de la page départementale
  dept_html <- tags$html(
    tags$head(
      tags$style(HTML("
        html, body { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; font-family: 'Segoe UI', sans-serif; background: #f8fafc; }
        
        /* En-tête fixe en haut */
        #header { position: absolute; top: 0; left: 0; right: 0; height: 50px; background: #0f172a; color: white; padding: 0 25px; display: flex; align-items: center; justify-content: space-between; z-index: 1000; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        #header h1 { margin: 0; font-size: 18px; font-weight: 600; }
        .btn-back { background: #38bdf8; color: #0f172a; padding: 6px 14px; border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 13px; transition: background 0.2s; }
        .btn-back:hover { background: #0ea5e9; }
        
        /* Positionnement absolu de la carte pour éviter le bug de hauteur à 0px */
        #map-container { position: absolute; top: 50px; bottom: 0; left: 0; right: 0; width: 100%; }
        .html-widget, .leaflet { width: 100% !important; height: 100% !important; }
      "))
    ),
    tags$body(
      tags$div(
        tags$div(id = "header",
                 tags$h1(paste("Département :", nom_dept, "(", code_dept, ")")),
                 tags$a(href = "index.html", class = "btn-back", "← Retour au menu national")
        ),
        tags$div(id = "map-container", m)
      )
    )
  )
  
  filename <- paste0("carte_dept_", code_dept, ".html")
  filepath <- paste0(DOSSIER, filename)
  save_html(dept_html, file = filepath)
  
  dep_files <- rbind(dep_files, data.frame(code = code_dept, nom = nom_dept, filename = filename, stringsAsFactors = FALSE))
}


# 5. GÉNÉRATION DE LA PAGE D'ACCUEIL (index.html pour GitHub Pages)
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
             tags$h1("Atlas des Lycées & Données INSEE - Sélection du Département"),
             tags$p(class = "subtitle", "Recherchez ou cliquez sur un département pour ouvrir sa carte interactive dédiée."),
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

# Nommé obligatoirement index.html pour GitHub Pages
index_path <- paste0(DOSSIER, "index.html")
save_html(index_html, file = index_path)

rstudioapi::viewer(normalizePath(index_path))
message("Génération complète terminée ! Fichiers prêts dans : ", DOSSIER)