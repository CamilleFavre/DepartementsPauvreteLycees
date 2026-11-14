library(sf); library(dplyr); library(ggplot2)

# ---------- 0) Contour France métropolitaine (Lambert-93) ----------
cache <- file.path(DOSSIER_OUT, "contour_france_2154.gpkg")
unlink(cache)                                    # supprime l'ancien cache vide

pts <- st_as_sf(d_sp, coords = c("x_pop", "y_pop"), crs = 2154)

# emprise Lambert-93 de la métropole (Corse incluse, DOM exclus)
bbox_met <- st_as_sfc(st_bbox(c(xmin = 80000,  xmax = 1300000,
                                ymin = 6000000, ymax = 7150000), crs = 2154))

# Option A : giscoR
fr <- tryCatch({
  if (!requireNamespace("giscoR", quietly = TRUE)) install.packages("giscoR")
  giscoR::gisco_get_countries(country = "FR", resolution = "10", year = "2020") %>%
    st_transform(2154) %>%
    st_make_valid() %>%
    st_intersection(bbox_met) %>%
    st_union() %>%
    st_simplify(dTolerance = 1500)
}, error = function(e) { message("giscoR indisponible : ", e$message); NULL })

# Option B (repli) : enveloppe des points de villes bufferisés
if (is.null(fr)) {
  fr <- tryCatch({
    u <- pts %>% st_buffer(20000) %>% st_union()
    tryCatch(st_concave_hull(u, ratio = 0.05), error = function(e) st_convex_hull(u)) %>%
      st_simplify(dTolerance = 1500)
  }, error = function(e) NULL)
}

# contrôle de santé : métropole ≈ 5.5e11 m² (~550 000 km²)
ok <- !is.null(fr) && !any(st_is_empty(fr)) && as.numeric(st_area(fr)) > 4e11
if (ok) {
  st_write(st_sf(geometry = st_sfc(st_geometry(fr)[[1]], crs = 2154)),
           cache, delete_dsn = TRUE, quiet = TRUE)
  cat("contour OK, aire km2 :", round(as.numeric(st_area(fr)) / 1e6), "- cache écrit\n")
} else {
  fr <- NULL
  cat("contour non produit -> carte sans outline\n")
}

# ---------- 1) Étiquettes et palette ----------
d_sp$bloc <- factor(d_sp$bloc)
lab <- d_sp %>%
  group_by(bloc) %>%
  summarise(x = mean(x_pop), y = mean(y_pop), n = n(), .groups = "drop")

if (!exists("pal10")) pal10 <- scales::hue_pal()(nlevels(d_sp$bloc))
pts$bloc <- d_sp$bloc

# ---------- 2) Carte ----------
p <- ggplot()
if (!is.null(fr)) {
  fr  <- st_transform(fr, st_crs(pts))
  bb  <- st_bbox(fr)
  p <- p + geom_sf(data = fr, fill = "grey93", color = "grey35", linewidth = 0.5)
}
p <- p +
  geom_sf(data = pts, aes(color = bloc), size = 1.4) +
  geom_text(data = lab, aes(x = x, y = y, label = paste0("B", bloc, " (n = ", n, ")")),
            size = 3.2, fontface = "bold", check_overlap = TRUE) +
  scale_color_manual(values = pal10, name = "k-means block (k = 10)") +
  theme_minimal(base_size = 12) +
  labs(title = "The 10 spatial CV blocks (k-means on population-weighted centroids, k = 10)",
       subtitle = "Each block serves once as the test fold: no neighbour of a test city is in the training set",
       x = NULL, y = NULL)

if (!is.null(fr)) {
  p <- p + coord_sf(xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]))
} else {
  p <- p + coord_sf()
}

print(p)
ggsave(file.path(DOSSIER_OUT, "carte_blocs_k10.png"), p, width = 10, height = 9, dpi = 300)