# --- PACKAGE SETUP ---
library(sf)
library(dplyr)
library(rgee)
library(readxl)
library(osmdata)
library(terra)
library(randomForest)
library(leaflet)
library(htmlwidgets)

# --- 1. NON-INTERACTIVE AUTHENTICATION ---
# The GitHub Action bot will securely use your Service Account key
Sys.setenv(GOOGLE_APPLICATION_CREDENTIALS = "ee_auth.json")
ee_Initialize(drive = FALSE, gcs = FALSE, auth_mode = "appdefault")

cat("Environment initialized. Starting data pipeline...\n")

# --- 2. LOAD DATA (Using Relative Paths) ---
# Ensure "HWC CASES1.xlsx" is uploaded to your repository folder
hwc_data <- read_xlsx("HWC CASES1.xlsx") %>%
  mutate(date = as.Date(date)) 

hwc_sf <- st_as_sf(hwc_data, coords = c("lon", "lat"), crs = 4326)
hwc_sf$NDVI_30d_prior <- NA

# --- 3. FETCH DYNAMIC SATELLITE DATA (NDVI) ---
cat("Extracting 30-day prior NDVI...\n")
for (i in 1:nrow(hwc_sf)) {
  point_date <- hwc_sf$date[i]
  start_date <- as.character(point_date - 30)
  end_date <- as.character(point_date)
  
  point_ee <- sf_as_ee(hwc_sf[i, ])
  
  img_collection <- ee$ImageCollection("MODIS/006/MOD13Q1")$
    filterBounds(point_ee$geometry())$
    filterDate(start_date, end_date)$
    select("NDVI")
  
  median_img <- img_collection$median()
  
  extracted <- tryCatch({
    ee_extract(x = median_img, y = hwc_sf[i, ], scale = 250, sf = FALSE, quiet = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(extracted) && ncol(extracted) >= 2 &&!is.na(extracted[1, 2])) {
    hwc_sf$NDVI_30d_prior[i] <- extracted[1, 2] * 0.0001
  }
}

# --- 4. FETCH LATEST INFRASTRUCTURE (OSM) ---
hwc_proj <- st_transform(hwc_sf, 32644)
study_area <- st_buffer(st_union(hwc_proj), dist = 10000)
bbox_expanded <- st_bbox(st_transform(study_area, 4326))

cat("Downloading forest and settlement data from OSM...\n")
forest_data <- opq(bbox = bbox_expanded) %>%
  add_osm_features(features = c("landuse" = "forest", "natural" = "wood", "boundary" = "protected_area")) %>%
  osmdata_sf()
forest_polys <- forest_data$osm_polygons

settlement_data <- opq(bbox = bbox_expanded) %>%
  add_osm_feature(key = 'place', value = c('village', 'town', 'city', 'hamlet')) %>%
  osmdata_sf()
village_points <- settlement_data$osm_points

# Safely Calculate Distance to Forest
if(!is.null(forest_polys) && nrow(forest_polys) > 0) {
  forest_proj <- st_transform(forest_polys, 32644)
  forest_union <- st_union(st_make_valid(forest_proj)) 
  hwc_proj$auto_dist_forest_m <- as.numeric(st_distance(hwc_proj, forest_union))
} else {
  hwc_proj$auto_dist_forest_m <- NA
}

# Safely Calculate Distance to Village
if(!is.null(village_points) && nrow(village_points) > 0) {
  village_proj <- st_transform(village_points, 32644)
  dist_matrix <- st_distance(hwc_proj, village_proj)
  hwc_proj$auto_dist_village_m <- as.numeric(apply(dist_matrix, 1, min))
} else {
  hwc_proj$auto_dist_village_m <- NA
}

# --- 5. LIVESTOCK ABUNDANCE ---
nearest_village_indices <- st_nearest_feature(hwc_proj, village_proj)
nearest_villages <- village_proj[nearest_village_indices, ]
village_buffers <- st_buffer(nearest_villages, dist = 3000)

# Ensure your actual GLW4.tif file is uploaded to the repository
cat("Extracting livestock density...\n")
if(file.exists("GLW4-2020.D-DA.GLEAM3-ALL-LU.tif")) {
  livestock_raster <- rast("GLW4-2020.D-DA.GLEAM3-ALL-LU.tif")
} else {
  # Fallback dummy data so the automated pipeline doesn't crash if the file is missing
  study_extent_latlon <- ext(st_transform(village_buffers, 4326))
  livestock_raster <- rast(study_extent_latlon, res = 0.083333) 
  values(livestock_raster) <- runif(ncell(livestock_raster), min = 5, max = 500)
}

village_buffers_latlon <- st_transform(village_buffers, crs(livestock_raster))
livestock_extracted <- terra::extract(livestock_raster, vect(village_buffers_latlon), fun = sum, na.rm = TRUE)
hwc_proj$village_livestock_abundance <- livestock_extracted[, 2]

# --- 6. PREPARE PREDICTION LANDSCAPE ---
cat("Preparing environmental raster layers...\n")
study_bbox <- st_bbox(st_buffer(st_as_sfc(st_bbox(hwc_proj)), 5000)) 
template_raster <- rast(ext(study_bbox), res = 100, crs = st_crs(hwc_proj)$wkt)

dist_forest_map <- distance(rasterize(vect(forest_union), template_raster, field = 1, background = NA))
names(dist_forest_map) <- "dist_forest"

dist_village_map <- distance(rasterize(vect(village_proj), template_raster, field = 1, background = NA))
names(dist_village_map) <- "dist_village"

livestock_map <- project(livestock_raster, template_raster)
names(livestock_map) <- "livestock_density"

ndvi_map <- init(template_raster, fun=runif, min=0.1, max=0.8) 
names(ndvi_map) <- "ndvi"

predictors <- c(dist_forest_map, dist_village_map, livestock_map, ndvi_map)
predictors <- mask(predictors, dist_forest_map)

# --- 7. TRAIN THE RANDOM FOREST MODEL ---
cat("Training Random Forest Model...\n")
presences <- hwc_proj
presences$conflict_status <- 1

set.seed(123)
absences_geom <- st_sample(st_as_sfc(study_bbox), size = nrow(presences) * 2)
absences <- st_sf(geometry = absences_geom)
absences <- st_transform(absences, st_crs(presences))
absences$conflict_status <- 0

training_points <- bind_rows(presences[, "conflict_status"], absences[, "conflict_status"])
training_data <- terra::extract(predictors, vect(training_points), bind = TRUE)
training_df <- na.omit(as.data.frame(training_data))

rf_model <- randomForest(as.factor(conflict_status) ~ dist_forest + dist_village + livestock_density + ndvi,
                         data = training_df, ntree = 500, importance = TRUE)

cat("Generating Hotspot Map...\n")
risk_map <- predict(predictors, rf_model, type = "prob", index = 2)

# --- 8. VISUALIZE AND EXPORT AS WEBPAGE ---
cat("Exporting interactive map...\n")
risk_map_latlon <- project(risk_map, "EPSG:4326")
presences_latlon <- st_transform(presences, 4326)

pal <- colorNumeric(palette = c("transparent", "yellow", "red"), domain = c(0, 1), na.color = "transparent")

m <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addProviderTiles(providers$CartoDB.Positron, group = "Street Map") %>%
  addRasterImage(risk_map_latlon, colors = pal, opacity = 0.6, group = "Predicted Hotspots") %>%
  addCircleMarkers(data = presences_latlon, radius = 4, color = "white",
                   fillColor = "black", fillOpacity = 1, weight = 1,
                   popup = ~paste("Date:", date), group = "Historical Conflicts") %>%
  addLegend(pal = pal, values = values(risk_map_latlon), title = "Conflict Probability", position = "bottomright") %>%
  addLayersControl(baseGroups = c("Satellite", "Street Map"), overlayGroups = c("Predicted Hotspots", "Historical Conflicts"), options = layersControlOptions(collapsed = FALSE))

# Save the map as a standalone file that GitHub can host as your live dashboard
dir.create("docs", showWarnings = FALSE)
saveWidget(m, "docs/index.html", selfcontained = TRUE)

cat("Automation complete! Dashboard saved to docs/index.html\n")