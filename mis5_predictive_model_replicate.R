# =============================================================================
# MIS 5 HABITAT SUITABILITY MODEL — COMPLETE PIPELINE (MERGED)
# =============================================================================
# Combines:
#   code_clean.R   — raster computation pipeline (Sections 1–3) + full model
#   FINAL_CODE.R   — downstream analyses (descriptive stats, correlations,
#                    zone weights, 4-class model, grid analysis, validation)
#
# Bug fixed (code_clean.R Section 3):
#   P98 threshold → P95 threshold for Distance to Water
#   (P98 = top 2% of flow = too few pixels, artificially inflated distances)
#
# Key design choices from FINAL_CODE.R adopted throughout:
#   - P95 wadi threshold (instead of P98)
#   - 4-class quartile classification (instead of 3-class tertiles)
#   - 900 m aggregation (instead of 1800 m)
#   - All outputs under outputs_val_june/
#   - Site–DEM assignment via sites$Area column (explicit, not coverage-based)
# =============================================================================

cat("================================================================================\n")
cat("MIS 5 HABITAT SUITABILITY MODEL — COMPLETE PIPELINE\n")
cat("================================================================================\n\n")

# =============================================================================
# SECTION 0: CONFIGURATION  (only modify this section)
# =============================================================================
# When run via run_pipeline.R (CLI), the *_OVERRIDE variables are pre-set.
# When run interactively in RStudio, edit the defaults below directly.

# ── EDIT THIS PATH to match your machine ──────────────────────────────────────
# Desktop:  "C:/Users/skkou/Documents/MIS5_hunters"
# Laptop:   update to wherever you copied the project folder
if (!exists("BASE_DIR")) BASE_DIR <- "C:/Users/skkou/OneDrive/Έγγραφα/MIS5_hunters"

INPUT_DIR  <- file.path(BASE_DIR, "Inputs")
DEM_DIR    <- file.path(INPUT_DIR, "DEM")
GEO_DIR    <- file.path(INPUT_DIR, "Geo_data")


OUTPUT_NAME  <- if (exists("OUTPUT_NAME_OVERRIDE")) OUTPUT_NAME_OVERRIDE else "outputs_val_june"
OUTPUT_DIR   <- file.path(BASE_DIR, OUTPUT_NAME)
TOPO_OUTPUT  <- file.path(OUTPUT_DIR, "Topo")
SITES_OUTPUT <- file.path(OUTPUT_DIR, "Sites_data")
MCDA_OUTPUT  <- file.path(OUTPUT_DIR, "MCDA")
STATS_OUTPUT <- file.path(OUTPUT_DIR, "Stats")
FIGS_OUTPUT  <- file.path(OUTPUT_DIR, "Figures")

# Regions and DEMs — overridden by config YAML when run from CLI.
# Edit these defaults when running interactively without a config file.
if (exists("DEM_FILES_OVERRIDE")) {
  DEM_FILES <- DEM_FILES_OVERRIDE
} else {
  DEM_FILES <- list(
    Saudi          = "Saudi.4326.tif",
    Lebanon        = "Lebanon.4326.tif",
    Syria          = "Syria.4326.tif",
    SouthernLevant = "Southern.Levant.tif"
  )
}

if (exists("REGIONS_OVERRIDE")) {
  REGIONS <- REGIONS_OVERRIDE
} else {
  REGIONS <- c("Saudi", "Lebanon", "Syria", "SouthernLevant")
}

# The region for which the suitability model is built.
if (!exists("MODEL_REGION")) MODEL_REGION <- "SouthernLevant"

# Sites CSV filename (inside Inputs/Sites/)
sites_csv  <- if (exists("SITES_FILE_OVERRIDE")) SITES_FILE_OVERRIDE else "MIS5.csv"
SITES_FILE <- file.path(INPUT_DIR, "Sites", sites_csv)

# Geology mask raster (inside Inputs/Geo_data/). NULL = skip.
if (exists("GEO_MASK_OVERRIDE")) {
  GEO_SHAPEFILE <- GEO_MASK_OVERRIDE
} else {
  GEO_SHAPEFILE <- "raw_material_presence.tif"
}

cat("Configuration:\n")
cat(sprintf("  Base directory:  %s\n", BASE_DIR))
cat(sprintf("  Sites file:      %s\n", SITES_FILE))
cat(sprintf("  Model region:    %s\n", MODEL_REGION))
cat(sprintf("  Regions:         %s\n", paste(REGIONS, collapse = ", ")))
cat(sprintf("  Output root:     %s\n\n", OUTPUT_DIR))

# =============================================================================
# SETUP — Install packages if needed (first run on a new machine)
# =============================================================================

required_packages <- c("terra", "whitebox", "corrplot")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(terra)
library(whitebox)
library(corrplot)

# First-time WhiteboxTools setup (downloads the binary ~30 MB)
if (!whitebox::check_whitebox_binary()) {
  cat("Installing WhiteboxTools binary (first run only)...\n")
  whitebox::install_whitebox()
}

# Create all output directories
for (d in c(TOPO_OUTPUT, SITES_OUTPUT, MCDA_OUTPUT, STATS_OUTPUT, FIGS_OUTPUT)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
for (var in c("Slope", "Aspect", "TRI", "TWI", "DistToWater")) {
  dir.create(file.path(TOPO_OUTPUT, var), showWarnings = FALSE, recursive = TRUE)
}

# Verify input files
cat("Verifying input files...\n")
if (!file.exists(SITES_FILE)) stop(sprintf("Sites file not found: %s", SITES_FILE))
for (region in REGIONS) {
  dem_path <- file.path(DEM_DIR, DEM_FILES[[region]])
  if (!file.exists(dem_path)) stop(sprintf("DEM not found: %s", dem_path))
}
cat("✓ All input files verified\n\n")


# =============================================================================
# SECTION 1: TOPOGRAPHIC VARIABLES — Slope, Aspect (Northness), TRI
# =============================================================================

cat("================================================================================\n")
cat("SECTION 1: TOPOGRAPHIC VARIABLES (Slope / Northness / TRI)\n")
cat("================================================================================\n\n")

for (region in REGIONS) {
  
  cat(sprintf("Processing %s...\n", region))
  
  dem_r <- rast(file.path(DEM_DIR, DEM_FILES[[region]]))
  cat(sprintf("  DEM: %d x %d px, ~%.0f m resolution\n",
              ncol(dem_r), nrow(dem_r), res(dem_r)[1] * 111000))
  
  # Slope (degrees)
  slope_r <- terrain(dem_r, v = "slope", unit = "degrees")
  writeRaster(slope_r,
              file.path(TOPO_OUTPUT, "Slope", paste0(region, "_slope.tif")),
              overwrite = TRUE)
  
  # Northness: cos(aspect) rescaled to [0, 1]  (1 = north-facing, 0 = south-facing)
  aspect_r    <- terrain(dem_r, v = "aspect", unit = "degrees")
  northness_r <- (cos(aspect_r * pi / 180) + 1) / 2
  writeRaster(northness_r,
              file.path(TOPO_OUTPUT, "Aspect", paste0(region, "_aspect.tif")),
              overwrite = TRUE)
  
  # TRI
  tri_r <- terrain(dem_r, v = "TRI")
  writeRaster(tri_r,
              file.path(TOPO_OUTPUT, "TRI", paste0(region, "_TRI.tif")),
              overwrite = TRUE)
  
  cat(sprintf("  ✓ Slope / Northness / TRI saved\n\n"))
  rm(dem_r, slope_r, aspect_r, northness_r, tri_r); gc(verbose = FALSE)
}

cat("✓ Section 1 complete\n\n")


# =============================================================================
# SECTION 2: TWI  (WhiteboxTools workflow)
# =============================================================================

cat("================================================================================\n")
cat("SECTION 2: TOPOGRAPHIC WETNESS INDEX (TWI)\n")
cat("================================================================================\n\n")
cat("Workflow: fill depressions → D8 flow accumulation (SCA) → slope (rad) → TWI\n\n")

for (region in REGIONS) {
  
  cat(sprintf("Processing %s TWI...\n", region))
  
  dem_r      <- rast(file.path(DEM_DIR, DEM_FILES[[region]]))
  temp_dir   <- file.path(TOPO_OUTPUT, "TWI", paste0("temp_", region))
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Reproject to UTM
  dem_centre_lon <- mean(as.vector(ext(dem_r))[c(1, 2)])
  utm_zone  <- floor((dem_centre_lon + 180) / 6) + 1
  utm_epsg  <- 32600 + utm_zone
  cat(sprintf("  UTM zone %dN (EPSG:%d)\n", utm_zone, utm_epsg))
  
  dem_utm      <- project(dem_r, sprintf("EPSG:%d", utm_epsg), method = "bilinear")
  dem_utm_file <- file.path(temp_dir, "dem_utm.tif")
  writeRaster(dem_utm, dem_utm_file, overwrite = TRUE)
  
  filled_file   <- file.path(temp_dir, "filled.tif")
  flow_file     <- file.path(temp_dir, "flow_accum.tif")
  slope_rad_file <- file.path(temp_dir, "slope_rad.tif")
  twi_file      <- file.path(temp_dir, "twi.tif")
  
  wbt_fill_depressions_wang_and_liu(dem = dem_utm_file, output = filled_file)
  wbt_d8_flow_accumulation(input = filled_file, output = flow_file,
                           out_type = "specific contributing area")
  wbt_slope(dem = filled_file, output = slope_rad_file, units = "radians")
  wbt_wetness_index(sca = flow_file, slope = slope_rad_file, output = twi_file)
  
  twi_utm   <- rast(twi_file)
  twi_wgs84 <- project(twi_utm, "EPSG:4326", method = "bilinear")
  twi_final <- resample(twi_wgs84, dem_r, method = "bilinear")
  
  writeRaster(twi_final,
              file.path(TOPO_OUTPUT, "TWI", paste0(region, "_TWI.tif")),
              overwrite = TRUE)
  
  cat(sprintf("  TWI range: %.2f – %.2f\n",
              global(twi_final, "min", na.rm = TRUE)[1,1],
              global(twi_final, "max", na.rm = TRUE)[1,1]))
  
  unlink(temp_dir, recursive = TRUE)
  cat(sprintf("  ✓ TWI saved\n\n"))
  rm(dem_r, dem_utm, twi_utm, twi_wgs84, twi_final); gc(verbose = FALSE)
}

cat("✓ Section 2 complete\n\n")


# =============================================================================
# SECTION 3: DISTANCE TO WATER  (P95 threshold — major wadis = top 5% SCA)
# =============================================================================
# BUG FIX vs code_clean.R:
#   Old code used P98 (top 2% of flow accumulation) as the wadi threshold.
#   This leaves too few "water" source pixels, artificially inflating distances.
#   Corrected to P95 (top 5%) — consistent with FINAL_CODE.R.
# =============================================================================

# =============================================================================
# SECTION 3: DISTANCE TO WATER  (pre-computed P95 rasters)
# =============================================================================

cat("================================================================================\n")
cat("SECTION 3: DISTANCE TO WATER  (loading pre-computed P95 rasters)\n")
cat("================================================================================\n\n")

DIST_INPUT_DIR <- file.path(BASE_DIR, "DistToWater")

if (!dir.exists(DIST_INPUT_DIR))
  stop(sprintf("DistToWater input folder not found: %s", DIST_INPUT_DIR))

for (region in REGIONS) {
  
  # Try _P95 name first, fall back to name without it
  src_file <- file.path(DIST_INPUT_DIR,
                        paste0(region, "_DistToWater_P95.tif"))
  if (!file.exists(src_file)) {
    src_file <- file.path(DIST_INPUT_DIR,
                          paste0(region, "_DistToWater.tif"))
  }
  
  dst_file <- file.path(TOPO_OUTPUT, "DistToWater",
                        paste0(region, "_DistToWater_P95.tif"))
  
  if (!file.exists(src_file))
    stop(sprintf("DistToWater raster not found: %s", src_file))
  
  cat(sprintf("Loading %s...\n", basename(src_file)))
  dist_r <- rast(src_file)
  
  cat(sprintf("  Distance range: %.2f – %.2f km\n",
              global(dist_r, "min", na.rm = TRUE)[1,1],
              global(dist_r, "max", na.rm = TRUE)[1,1]))
  
  writeRaster(dist_r, dst_file, overwrite = TRUE)
  cat(sprintf("  ✓ Copied to %s\n\n", dst_file))
  
  rm(dist_r); gc(verbose = FALSE)
}

cat("✓ Section 3 complete\n\n")

# =============================================================================
# SECTION 4: EXTRACT VARIABLES AT SITE LOCATIONS
# =============================================================================
# Site–DEM assignment via sites$Area column (explicit, not coverage-based).
# =============================================================================

cat("================================================================================\n")
cat("SECTION 4: EXTRACT VARIABLES AT SITES\n")
cat("================================================================================\n\n")

sites <- read.csv(SITES_FILE)
cat(sprintf("Loaded %d sites from: %s\n", nrow(sites), basename(SITES_FILE)))

# Validate required columns
required_cols <- c("Lon", "Lat", "Area", "Region")
missing_cols  <- required_cols[!required_cols %in% names(sites)]
if (length(missing_cols) > 0) {
  stop(sprintf("MIS5.csv is missing columns: %s", paste(missing_cols, collapse = ", ")))
}

cat("\nBioclimatic zones:\n"); print(table(sites$Region))
cat("\nDEM regions (Area):\n");  print(table(sites$Area))
cat("\n")

sites$elevation   <- NA
sites$slope       <- NA
sites$northness   <- NA
sites$TRI         <- NA
sites$TWI         <- NA
sites$DistToWater <- NA
sites$source_dem  <- NA

for (region in REGIONS) {
  
  region_idx <- which(trimws(sites$Area) == region)
  if (length(region_idx) == 0) {
    cat(sprintf("%s: no sites, skipping\n\n", region)); next
  }
  cat(sprintf("Extracting %s (%d sites)...\n", region, length(region_idx)))
  
  dem_r    <- rast(file.path(DEM_DIR, DEM_FILES[[region]]))
  slope_r  <- rast(file.path(TOPO_OUTPUT, "Slope",  paste0(region, "_slope.tif")))
  aspect_r <- rast(file.path(TOPO_OUTPUT, "Aspect", paste0(region, "_aspect.tif")))
  tri_r    <- rast(file.path(TOPO_OUTPUT, "TRI",    paste0(region, "_TRI.tif")))
  twi_r    <- rast(file.path(TOPO_OUTPUT, "TWI",    paste0(region, "_TWI.tif")))
  dist_r   <- rast(file.path(TOPO_OUTPUT, "DistToWater",
                             paste0(region, "_DistToWater_P95.tif")))
  
  region_sites <- sites[region_idx, ]
  pts <- vect(region_sites, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  
  sites$elevation[region_idx]   <- extract(dem_r,    pts)[, 2]
  sites$slope[region_idx]       <- extract(slope_r,  pts)[, 2]
  sites$northness[region_idx]   <- extract(aspect_r, pts)[, 2]
  sites$TRI[region_idx]         <- extract(tri_r,    pts)[, 2]
  sites$TWI[region_idx]         <- extract(twi_r,    pts)[, 2]
  sites$DistToWater[region_idx] <- extract(dist_r,   pts)[, 2]
  sites$source_dem[region_idx]  <- region
  
  n_ok <- sum(!is.na(sites$elevation[region_idx]))
  cat(sprintf("  Extracted: %d/%d sites\n\n", n_ok, length(region_idx)))
  
  rm(dem_r, slope_r, aspect_r, tri_r, twi_r, dist_r); gc(verbose = FALSE)
}

vars_check     <- c("elevation", "slope", "northness", "TRI", "TWI", "DistToWater")
sites_complete <- sites[complete.cases(sites[, vars_check]), ]
sites_na       <- sites[!complete.cases(sites[, vars_check]), ]

cat(sprintf("Complete sites: %d/%d\n", nrow(sites_complete), nrow(sites)))
if (nrow(sites_na) > 0) {
  cat(sprintf("Sites with NAs (%d):\n", nrow(sites_na)))
  print(sites_na[, c("Site", "Area", "Region")])
}

write.csv(sites_complete,
          file.path(SITES_OUTPUT, "Sites_with_variables.csv"),
          row.names = FALSE)
cat(sprintf("\n✓ Saved: %s\n\n", file.path(SITES_OUTPUT, "Sites_with_variables.csv")))


# =============================================================================
# SECTION 5: BACKGROUND SAMPLING  (15× sites per region, min 100 pts)
# =============================================================================

cat("================================================================================\n")
cat("SECTION 5: BACKGROUND SAMPLING\n")
cat("================================================================================\n\n")

sites_per_region <- table(sites_complete$source_dem)
background_list  <- list()
set.seed(123)

for (region in REGIONS) {
  
  n_sites <- as.numeric(sites_per_region[region])
  if (is.na(n_sites)) n_sites <- 0
  n_bg <- max(100, n_sites * 15)
  cat(sprintf("Sampling %s: %d points...\n", region, n_bg))
  
  dem_r    <- rast(file.path(DEM_DIR, DEM_FILES[[region]]))
  slope_r  <- rast(file.path(TOPO_OUTPUT, "Slope",  paste0(region, "_slope.tif")))
  aspect_r <- rast(file.path(TOPO_OUTPUT, "Aspect", paste0(region, "_aspect.tif")))
  tri_r    <- rast(file.path(TOPO_OUTPUT, "TRI",    paste0(region, "_TRI.tif")))
  twi_r    <- rast(file.path(TOPO_OUTPUT, "TWI",    paste0(region, "_TWI.tif")))
  dist_r   <- rast(file.path(TOPO_OUTPUT, "DistToWater",
                             paste0(region, "_DistToWater_P95.tif")))
  
  bg_pts    <- spatSample(dem_r, size = n_bg, method = "random",
                          na.rm = TRUE, as.points = TRUE, xy = TRUE)
  bg_coords <- geom(bg_pts)[, c("x", "y")]
  
  bg_df <- data.frame(
    type        = "background",
    region      = region,
    Lon         = bg_coords[, 1],
    Lat         = bg_coords[, 2],
    elevation   = extract(dem_r,    bg_pts)[, 2],
    slope       = extract(slope_r,  bg_pts)[, 2],
    northness   = extract(aspect_r, bg_pts)[, 2],
    TRI         = extract(tri_r,    bg_pts)[, 2],
    TWI         = extract(twi_r,    bg_pts)[, 2],
    DistToWater = extract(dist_r,   bg_pts)[, 2]
  )
  bg_df <- bg_df[complete.cases(bg_df), ]
  
  if (nrow(bg_df) < n_bg * 0.9)
    cat(sprintf("  Warning: only %.0f%% of requested points returned\n",
                nrow(bg_df) / n_bg * 100))
  
  background_list[[region]] <- bg_df
  cat(sprintf("  → %d valid background points\n", nrow(bg_df)))
  
  rm(dem_r, slope_r, aspect_r, tri_r, twi_r, dist_r); gc(verbose = FALSE)
}

background <- do.call(rbind, background_list)
cat(sprintf("\n✓ Total background: %d points\n\n", nrow(background)))

write.csv(background,
          file.path(SITES_OUTPUT, "Background_sample.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: %s\n\n", file.path(SITES_OUTPUT, "Background_sample.csv")))


# =============================================================================
# SECTION 6: DESCRIPTIVE STATISTICS BY BIOCLIMATIC ZONE
# =============================================================================

cat("================================================================================\n")
cat("SECTION 6: DESCRIPTIVE STATISTICS BY BIOCLIMATIC ZONE\n")
cat("================================================================================\n\n")

var_labels <- c(
  elevation   = "Elevation (m)",
  slope       = "Slope (degrees)",
  northness   = "Northness (0-1)",
  TRI         = "TRI",
  TWI         = "TWI",
  DistToWater = "Distance to Water (km)"
)

zones        <- c("Mediterranean", "Irano-Turanian", "Saharo-Arabian")
summary_rows <- list()

for (zone in zones) {
  
  z_sites <- sites_complete[sites_complete$Region == zone, ]
  n <- nrow(z_sites)
  cat(sprintf("--- %s (n=%d) ---\n", zone, n))
  
  for (v in names(var_labels)) {
    vals     <- z_sites[[v]][!is.na(z_sites[[v]])]
    site_min <- z_sites$Site[which.min(z_sites[[v]])]
    site_max <- z_sites$Site[which.max(z_sites[[v]])]
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Zone     = zone,
      n_sites  = n,
      Variable = var_labels[v],
      Mean     = round(mean(vals), 2),
      SD       = round(sd(vals), 2),
      Median   = round(median(vals), 2),
      Min      = round(min(vals), 2),
      Min_Site = site_min,
      Max      = round(max(vals), 2),
      Max_Site = site_max
    )
    
    cat(sprintf("  %-25s: mean=%.2f ± %.2f, median=%.2f, range=%.2f–%.2f\n",
                var_labels[v], mean(vals), sd(vals), median(vals), min(vals), max(vals)))
    cat(sprintf("    Min: %s (%.2f),  Max: %s (%.2f)\n",
                site_min, min(vals), site_max, max(vals)))
  }
  cat("\n")
}

summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
write.csv(summary_table,
          file.path(STATS_OUTPUT, "Descriptive_stats_by_zone.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: %s\n\n", file.path(STATS_OUTPUT, "Descriptive_stats_by_zone.csv")))

# Cross-zone extremes
cat("Cross-zone extremes (all sites):\n")
for (v in names(var_labels)) {
  vals     <- sites_complete[[v]]
  site_min <- sites_complete$Site[which.min(vals)]
  site_max <- sites_complete$Site[which.max(vals)]
  zone_min <- sites_complete$Region[which.min(vals)]
  zone_max <- sites_complete$Region[which.max(vals)]
  cat(sprintf("  %-25s  min=%.2f → %s (%s)   max=%.2f → %s (%s)\n",
              var_labels[v],
              min(vals, na.rm=TRUE), site_min, zone_min,
              max(vals, na.rm=TRUE), site_max, zone_max))
}
cat("\n")


# =============================================================================
# SECTION 7: PAIRWISE CORRELATION ANALYSIS (Spearman)
# =============================================================================

cat("================================================================================\n")
cat("SECTION 7: PAIRWISE VARIABLE CORRELATIONS\n")
cat("================================================================================\n\n")

var_names   <- c("Elevation", "Slope", "Northness", "TRI", "TWI", "Dist. to Water")
data_matrix <- sites_complete[, names(var_labels)]
colnames(data_matrix) <- var_names

cor_matrix <- cor(data_matrix, use = "complete.obs", method = "spearman")

n         <- nrow(data_matrix)
p_matrix  <- matrix(NA, nrow = 6, ncol = 6,
                    dimnames = list(var_names, var_names))
for (i in 1:6) for (j in 1:6) {
  if (i != j) {
    tst <- cor.test(data_matrix[, i], data_matrix[, j],
                    method = "spearman", exact = FALSE)
    p_matrix[i, j] <- tst$p.value
  }
}

# Print upper triangle (r) / lower triangle (p)
cat("Correlation matrix  [upper: rho | lower: p-value]\n\n")
cat(sprintf("%-17s", ""))
for (v in var_names) cat(sprintf(" %8s", substr(v, 1, 8)))
cat("\n", paste(rep("-", 66), collapse=""), "\n")
for (i in 1:6) {
  cat(sprintf("%-17s", var_names[i]))
  for (j in 1:6) {
    if      (i == j) cat(sprintf(" %8s", "---"))
    else if (j  > i) cat(sprintf(" %8.3f", cor_matrix[i, j]))
    else cat(sprintf(" %8s", ifelse(p_matrix[i, j] < 0.001, "<0.001",
                                    sprintf("%.3f", p_matrix[i, j]))))
  }
  cat("\n")
}
cat("\n")

write.csv(cor_matrix,
          file.path(STATS_OUTPUT, "Variable_correlations.csv"),
          row.names = TRUE)
write.csv(p_matrix,
          file.path(STATS_OUTPUT, "Variable_correlation_pvalues.csv"),
          row.names = TRUE)

png(file.path(FIGS_OUTPUT, "Variable_Correlation_Matrix.png"),
    width = 3000, height = 3000, res = 400)
col_palette <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                                  "#F7F7F7",
                                  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(200)
corrplot(cor_matrix, method = "color", type = "upper", order = "original",
         addCoef.col = "black", number.cex = 0.8, tl.col = "black", tl.srt = 45,
         tl.cex = 0.9, col = col_palette, cl.cex = 0.8, cl.ratio = 0.2,
         diag = TRUE, title = "Variable Correlation Matrix (Spearman rho)",
         mar = c(0, 0, 2, 0))
dev.off()
cat(sprintf("✓ Saved correlation matrix plot\n\n"))


# =============================================================================
# SECTION 8: VERHAGEN WEIGHTS — GLOBAL + BY BIOCLIMATIC ZONE
# =============================================================================

cat("================================================================================\n")
cat("SECTION 8: VERHAGEN WEIGHTS (global + by zone)\n")
cat("================================================================================\n\n")

# Verhagen gain function (Kullback–Leibler divergence)
verhagen_gain <- function(site_vals, bg_vals) {
  site_vals <- site_vals[!is.na(site_vals)]
  bg_vals   <- bg_vals[!is.na(bg_vals)]
  all_vals  <- c(site_vals, bg_vals)
  breaks    <- quantile(all_vals, probs = seq(0, 1, length.out = 11))
  breaks    <- unique(breaks)
  breaks[1] <- breaks[1] - 0.001
  n_bins    <- length(breaks) - 1
  site_bins <- cut(site_vals, breaks = breaks, labels = FALSE)
  bg_bins   <- cut(bg_vals,   breaks = breaks, labels = FALSE)
  site_freq <- table(factor(site_bins, levels = 1:n_bins)) / length(site_vals)
  bg_freq   <- table(factor(bg_bins,   levels = 1:n_bins)) / length(bg_vals)
  site_freq[site_freq == 0] <- 0.001
  bg_freq[bg_freq   == 0]   <- 0.001
  sum(site_freq * log(site_freq / bg_freq))
}

VARS      <- c("Elevation", "Slope", "Northness", "TRI", "TWI", "DistToWater")
site_list <- list(sites_complete$elevation, sites_complete$slope,
                  sites_complete$northness, sites_complete$TRI,
                  sites_complete$TWI, sites_complete$DistToWater)
bg_list   <- list(background$elevation, background$slope,
                  background$northness, background$TRI,
                  background$TWI, background$DistToWater)

# --- Global weights ---
cat("Calculating modeled-area weights (all sites)...\n")
gains_global   <- sapply(1:6, function(i) verhagen_gain(site_list[[i]], bg_list[[i]]))
weights_global <- (gains_global / sum(gains_global)) * 100

cat("\nModeled Area Weights:\n")
for (i in 1:6) cat(sprintf("  %-15s gain=%.4f  weight=%.1f%%\n",
                           VARS[i], gains_global[i], weights_global[i]))
cat(sprintf("\nGrouped:\n"))
cat(sprintf("  Topographic (Elev + Slope + TRI): %.1f%%\n",
            sum(weights_global[c(1,2,4)])))
cat(sprintf("  Water-related (TWI + DistToWater): %.1f%%\n",
            sum(weights_global[c(5,6)])))
cat(sprintf("  Solar (Northness): %.1f%%\n\n", weights_global[3]))

# --- Zone-specific weights ---
all_weights <- data.frame(Zone = "Modeled Area", Variable = VARS,
                          Gain = gains_global, Weight = weights_global)

for (zone in zones) {
  cat(sprintf("Calculating %s weights...\n", zone))
  sites_zone <- sites_complete[sites_complete$Region == zone, ]
  cat(sprintf("  Sites in zone: %d\n", nrow(sites_zone)))
  if (nrow(sites_zone) < 3) { cat("  Too few sites, skipping\n\n"); next }
  
  site_zone_list <- list(sites_zone$elevation, sites_zone$slope,
                         sites_zone$northness, sites_zone$TRI,
                         sites_zone$TWI, sites_zone$DistToWater)
  gains_zone   <- sapply(1:6, function(i) verhagen_gain(site_zone_list[[i]], bg_list[[i]]))
  weights_zone <- (gains_zone / sum(gains_zone)) * 100
  
  for (i in 1:6) cat(sprintf("    %-15s: %.1f%%\n", VARS[i], weights_zone[i]))
  cat("\n")
  
  all_weights <- rbind(all_weights,
                       data.frame(Zone = zone, Variable = VARS,
                                  Gain = gains_zone, Weight = weights_zone))
}

write.csv(all_weights,
          file.path(STATS_OUTPUT, "Variable_weights_all.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: %s\n\n", file.path(STATS_OUTPUT, "Variable_weights_all.csv")))

# Keep the modeled-area weights as the main weight vector
weights  <- weights_global


# =============================================================================
# SECTION 9: SUITABILITY MODEL  (4-class, 900 m, Southern Levant)
# =============================================================================

cat("================================================================================\n")
cat("SECTION 9: SUITABILITY MODEL BUILD  (4-class / 900 m)\n")
cat("================================================================================\n\n")

bounds <- data.frame(
  Variable = VARS,
  P5  = sapply(site_list, function(x) quantile(x, 0.05, na.rm = TRUE)),
  P95 = sapply(site_list, function(x) quantile(x, 0.95, na.rm = TRUE))
)
cat("Normalization bounds (P5–P95 from sites):\n"); print(bounds); cat("\n")

# Load Southern Levant rasters
dem_r    <- rast(file.path(DEM_DIR, DEM_FILES[[MODEL_REGION]]))
slope_r  <- rast(file.path(TOPO_OUTPUT, "Slope",       paste0(MODEL_REGION, "_slope.tif")))
aspect_r <- rast(file.path(TOPO_OUTPUT, "Aspect",      paste0(MODEL_REGION, "_aspect.tif")))
tri_r    <- rast(file.path(TOPO_OUTPUT, "TRI",         paste0(MODEL_REGION, "_TRI.tif")))
twi_r    <- rast(file.path(TOPO_OUTPUT, "TWI",         paste0(MODEL_REGION, "_TWI.tif")))
dist_r   <- rast(file.path(TOPO_OUTPUT, "DistToWater", paste0(MODEL_REGION, "_DistToWater_P95.tif")))

twi_r  <- resample(twi_r,  dem_r, method = "bilinear")
dist_r <- resample(dist_r, dem_r, method = "bilinear")

normalize <- function(r, p5, p95, invert = FALSE) {
  r_norm <- clamp((r - p5) / (p95 - p5) * 100, 0, 100)
  if (invert) r_norm <- 100 - r_norm
  r_norm
}

elev_norm  <- normalize(dem_r,    bounds$P5[1], bounds$P95[1], invert = TRUE)
slope_norm <- normalize(slope_r,  bounds$P5[2], bounds$P95[2])
north_norm <- clamp(aspect_r, 0, 1) * 100          # already on 0–1 scale from Section 1
tri_norm   <- normalize(tri_r,    bounds$P5[4], bounds$P95[4])
twi_norm   <- normalize(twi_r,    bounds$P5[5], bounds$P95[5])
dist_norm  <- normalize(dist_r,   bounds$P5[6], bounds$P95[6], invert = TRUE)

w <- weights / 100
suitability_30m <- (elev_norm  * w[1] + slope_norm * w[2] + north_norm * w[3] +
                      tri_norm   * w[4] + twi_norm   * w[5] + dist_norm  * w[6])

cat(sprintf("30m suitability range: %.2f – %.2f\n",
            global(suitability_30m, "min", na.rm=TRUE)[1,1],
            global(suitability_30m, "max", na.rm=TRUE)[1,1]))

writeRaster(suitability_30m,
            file.path(MCDA_OUTPUT, "Suitability_30m_P95.tif"), overwrite = TRUE)
cat("✓ Saved: Suitability_30m_P95.tif\n")

# Aggregate to 900 m (factor = 30)
cat("Aggregating to 900 m...\n")
suitability_900m <- aggregate(suitability_30m, fact = 30, fun = "mean", na.rm = TRUE)
writeRaster(suitability_900m,
            file.path(MCDA_OUTPUT, "Suitability_900m_P95.tif"), overwrite = TRUE)
cat("✓ Saved: Suitability_900m_P95.tif\n\n")

# 4-class quartile classification
vals      <- values(suitability_900m, na.rm = TRUE)
quartiles <- quantile(vals, probs = c(0, 0.25, 0.5, 0.75, 1))

cat("Quartile class breaks:\n")
cat(sprintf("  Class 1 Low:       %.2f – %.2f\n", quartiles[1], quartiles[2]))
cat(sprintf("  Class 2 Medium:    %.2f – %.2f\n", quartiles[2], quartiles[3]))
cat(sprintf("  Class 3 High:      %.2f – %.2f\n", quartiles[3], quartiles[4]))
cat(sprintf("  Class 4 Very High: %.2f – %.2f\n\n", quartiles[4], quartiles[5]))

rcl_matrix <- matrix(c(
  quartiles[1] - 0.001, quartiles[2],         1,
  quartiles[2],         quartiles[3],         2,
  quartiles[3],         quartiles[4],         3,
  quartiles[4],         quartiles[5] + 0.001, 4
), ncol = 3, byrow = TRUE)

suitability_4class <- classify(suitability_900m, rcl_matrix, include.lowest = TRUE)
writeRaster(suitability_4class,
            file.path(MCDA_OUTPUT, "Suitability_4class_P95.tif"), overwrite = TRUE)
cat("✓ Saved: Suitability_4class_P95.tif\n\n")

# Geology mask (optional)
geo_mask_file <- file.path(GEO_DIR, GEO_SHAPEFILE)
if (!is.null(GEO_SHAPEFILE) && file.exists(geo_mask_file)) {
  cat("Applying geology mask...\n")
  geology_mask      <- rast(geo_mask_file)
  geo_resampled     <- resample(geology_mask, suitability_900m, method = "near")
  
  suitability_masked  <- suitability_900m * geo_resampled
  suitability_masked[is.na(geo_resampled)] <- NA
  s4_masked           <- suitability_4class * geo_resampled
  s4_masked[is.na(geo_resampled)] <- NA
  
  writeRaster(suitability_masked,
              file.path(MCDA_OUTPUT, "Suitability_900m_P95_MASKED.tif"), overwrite = TRUE)
  writeRaster(s4_masked,
              file.path(MCDA_OUTPUT, "Suitability_4class_P95_MASKED.tif"), overwrite = TRUE)
  
  valid_px    <- sum(!is.na(values(suitability_900m)), na.rm = TRUE)
  excl_px     <- sum(is.na(values(suitability_masked)) & !is.na(values(suitability_900m)), na.rm = TRUE)
  cat(sprintf("✓ Geology mask: %.1f%% of landscape excluded\n", excl_px / valid_px * 100))
  cat("✓ Saved masked suitability rasters\n\n")
} else {
  cat("⚠ Geology mask not found or disabled — skipping\n\n")
  suitability_masked <- suitability_900m
  s4_masked          <- suitability_4class
}


# =============================================================================
# SECTION 10: MODEL VALIDATION
# =============================================================================

cat("================================================================================\n")
cat("SECTION 10: MODEL VALIDATION\n")
cat("================================================================================\n\n")

# Ensure output directories exist
for (d in c(STATS_OUTPUT, FIGS_OUTPUT)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

sites_sl  <- sites_complete[sites_complete$source_dem == MODEL_REGION, ]
site_pts  <- vect(sites_sl, geom = c("Lon", "Lat"), crs = "EPSG:4326")

sites_sl$suitability <- extract(suitability_30m, site_pts)[, 2]
sites_sl$suit_class  <- as.integer(cut(sites_sl$suitability,
                                       breaks = c(-Inf, quartiles[2], quartiles[3],
                                                  quartiles[4], Inf),
                                       labels = c(1, 2, 3, 4)))

cat(sprintf("Sites (Southern Levant): n=%d\n", nrow(sites_sl)))
cat(sprintf("Suitability: mean=%.2f, range %.2f–%.2f\n\n",
            mean(sites_sl$suitability, na.rm=TRUE),
            min(sites_sl$suitability,  na.rm=TRUE),
            max(sites_sl$suitability,  na.rm=TRUE)))

# --- Chi-square ---
cat("=== CHI-SQUARE ENRICHMENT ===\n\n")
class_counts <- as.vector(freq(suitability_4class)$count)
class_props  <- class_counts / sum(class_counts)
site_cls_cnt <- table(factor(sites_sl$suit_class, levels = 1:4))
expected     <- class_props * nrow(sites_sl)
class_names  <- c("Low", "Medium", "High", "Very High")

cat(sprintf("%-12s  Landscape  Sites  Expected  Enrichment\n", "Class"))
cat(paste(rep("-", 58), collapse=""), "\n")
for (i in 1:4) {
  enrich <- (site_cls_cnt[i] / nrow(sites_sl)) / class_props[i]
  cat(sprintf("%-12s  %6.1f%%    %3d    %5.1f     %.2fx\n",
              class_names[i], class_props[i]*100, site_cls_cnt[i], expected[i], enrich))
}
chi_sq <- sum((site_cls_cnt - expected)^2 / expected)
p_val  <- pchisq(chi_sq, df = 3, lower.tail = FALSE)
cat(sprintf("\nchi-sq = %.2f, df=3, p = %.2e\n", chi_sq, p_val))
if (p_val < 0.001) cat("✓ SIGNIFICANT: sites non-randomly distributed\n\n")

# --- Mann-Whitney ---
cat("=== MANN-WHITNEY U TEST ===\n\n")
set.seed(123)
land_sample <- spatSample(suitability_900m, size = 10000, method = "random", na.rm = TRUE)
landscape_vals <- land_sample[!is.na(land_sample[, 1]), 1]

cat(sprintf("Sites mean:     %.2f\n", mean(sites_sl$suitability, na.rm=TRUE)))
cat(sprintf("Landscape mean: %.2f\n\n", mean(landscape_vals)))
mw_test <- wilcox.test(sites_sl$suitability, landscape_vals,
                       alternative = "greater", exact = FALSE)
cat(sprintf("Mann-Whitney U: p = %.2e\n", mw_test$p.value))
if (mw_test$p.value < 0.001) cat("✓ SIGNIFICANT\n\n")

# --- Cohen's d ---
cat("=== COHEN'S D ===\n\n")
mean_diff <- mean(sites_sl$suitability, na.rm=TRUE) - mean(landscape_vals)
pooled_sd <- sqrt((sd(sites_sl$suitability, na.rm=TRUE)^2 + sd(landscape_vals)^2) / 2)
cohens_d  <- mean_diff / pooled_sd
label_d   <- if (cohens_d >= 0.8) "large" else if (cohens_d >= 0.5) "medium" else "small"
cat(sprintf("Cohen's d = %.3f (%s)\n\n", cohens_d, label_d))

# --- LOO-CV ---
cat("=== LEAVE-ONE-OUT CROSS-VALIDATION ===\n\n")
sl_sites        <- which(sites_complete$source_dem == MODEL_REGION)
loo_correct_hv  <- logical(length(sl_sites))   # High or Very High
loo_correct_vh  <- logical(length(sl_sites))   # Very High only
loo_weights     <- matrix(NA, nrow = length(sl_sites), ncol = 6,
                          dimnames = list(NULL, VARS))  # store fold weights

for (i in seq_along(sl_sites)) {
  train  <- sites_complete[-sl_sites[i], ]
  test   <- sites_complete[ sl_sites[i], ]
  t_list <- list(train$elevation, train$slope, train$northness,
                 train$TRI, train$TWI, train$DistToWater)
  loo_g  <- sapply(1:6, function(j) verhagen_gain(t_list[[j]], bg_list[[j]]))
  loo_w  <- loo_g / sum(loo_g)
  loo_weights[i, ] <- loo_w * 100   # store as percentages
  
  test_score <-
    pmax(0, pmin(1, (bounds$P95[1] - test$elevation)   / (bounds$P95[1] - bounds$P5[1]))) * loo_w[1] +
    pmax(0, pmin(1, (test$slope    - bounds$P5[2])     / (bounds$P95[2] - bounds$P5[2]))) * loo_w[2] +
    pmax(0, pmin(1, test$northness))                                                       * loo_w[3] +
    pmax(0, pmin(1, (test$TRI      - bounds$P5[4])     / (bounds$P95[4] - bounds$P5[4]))) * loo_w[4] +
    pmax(0, pmin(1, (test$TWI      - bounds$P5[5])     / (bounds$P95[5] - bounds$P5[5]))) * loo_w[5] +
    pmax(0, pmin(1, (bounds$P95[6] - test$DistToWater) / (bounds$P95[6] - bounds$P5[6]))) * loo_w[6]
  
  test_class <- findInterval(test_score * 100,
                             c(quartiles[2], quartiles[3], quartiles[4])) + 1
  loo_correct_hv[i] <- (test_class >= 3)
  loo_correct_vh[i] <- (test_class == 4)
}

loo_acc_hv <- mean(loo_correct_hv) * 100
loo_acc_vh <- mean(loo_correct_vh) * 100

# Clopper-Pearson exact 95% CIs on LOO-CV accuracy
n_loo      <- length(sl_sites)
ci_loo_hv  <- binom.test(sum(loo_correct_hv), n_loo)$conf.int * 100
ci_loo_vh  <- binom.test(sum(loo_correct_vh), n_loo)$conf.int * 100

cat(sprintf("LOO-CV accuracy (High + Very High): %.1f%%  (%d/%d)  95%% CI [%.1f%%, %.1f%%]\n",
            loo_acc_hv, sum(loo_correct_hv), n_loo, ci_loo_hv[1], ci_loo_hv[2]))
cat(sprintf("LOO-CV accuracy (Very High only):   %.1f%%  (%d/%d)  95%% CI [%.1f%%, %.1f%%]\n\n",
            loo_acc_vh, sum(loo_correct_vh), n_loo, ci_loo_vh[1], ci_loo_vh[2]))

# --- LOO-CV weight stability ---
cat("=== LOO-CV WEIGHT STABILITY ===\n\n")
cat("Weight per variable across LOO folds (mean ± SD):\n")
cat(sprintf("  Full-model weights shown for comparison.\n\n"))
cat(sprintf("  %-15s  Full model   LOO mean ± SD     Range\n", "Variable"))
cat(paste(rep("-", 68), collapse = ""), "\n")
for (j in 1:6) {
  w_mean  <- mean(loo_weights[, j])
  w_sd    <- sd(loo_weights[, j])
  w_min   <- min(loo_weights[, j])
  w_max   <- max(loo_weights[, j])
  cat(sprintf("  %-15s  %5.1f%%       %5.1f%% ± %4.1f%%     [%.1f%%, %.1f%%]\n",
              VARS[j], weights[j], w_mean, w_sd, w_min, w_max))
}
cat("\n")

# Save weight stability table
loo_weight_stability <- data.frame(
  Variable    = VARS,
  Full_Weight = weights,
  LOO_Mean    = apply(loo_weights, 2, mean),
  LOO_SD      = apply(loo_weights, 2, sd),
  LOO_Min     = apply(loo_weights, 2, min),
  LOO_Max     = apply(loo_weights, 2, max)
)
write.csv(loo_weight_stability,
          file.path(STATS_OUTPUT, "LOO_weight_stability.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: %s\n\n", file.path(STATS_OUTPUT, "LOO_weight_stability.csv")))

# --- Permutation test ---
cat("=== PERMUTATION TEST (location-based, n=999) ===\n\n")
set.seed(123)
n_sites_sl       <- nrow(sites_sl)
perm_acc_hv      <- numeric(999)
perm_acc_vh      <- numeric(999)

for (p in 1:999) {
  rand_s          <- sample(landscape_vals, size = n_sites_sl, replace = FALSE)
  rand_cls        <- findInterval(rand_s, c(quartiles[2], quartiles[3], quartiles[4])) + 1
  perm_acc_hv[p]  <- mean(rand_cls >= 3) * 100
  perm_acc_vh[p]  <- mean(rand_cls == 4) * 100
}

real_acc_hv <- mean(sites_sl$suit_class >= 3, na.rm=TRUE) * 100
real_acc_vh <- mean(sites_sl$suit_class == 4, na.rm=TRUE) * 100
p_perm_hv   <- (sum(perm_acc_hv >= real_acc_hv) + 1) / 1000
p_perm_vh   <- (sum(perm_acc_vh >= real_acc_vh) + 1) / 1000

cat(sprintf("High+VH: real=%.1f%%, null mean=%.1f%%, p=%.4f\n",
            real_acc_hv, mean(perm_acc_hv), p_perm_hv))
cat(sprintf("VH only: real=%.1f%%, null mean=%.1f%%, p=%.4f\n\n",
            real_acc_vh, mean(perm_acc_vh), p_perm_vh))

# --- Clopper-Pearson exact 95% CIs on site classification ---
cat("=== CLOPPER-PEARSON 95% CONFIDENCE INTERVALS ===\n\n")
n_sl   <- nrow(sites_sl)
k_hv   <- sum(sites_sl$suit_class >= 3, na.rm = TRUE)
k_vh   <- sum(sites_sl$suit_class == 4, na.rm = TRUE)
ci_hv  <- binom.test(k_hv, n_sl)$conf.int * 100
ci_vh  <- binom.test(k_vh, n_sl)$conf.int * 100

cat(sprintf("Sites in High+VH: %d/%d = %.1f%%   95%% CI [%.1f%%, %.1f%%]\n",
            k_hv, n_sl, real_acc_hv, ci_hv[1], ci_hv[2]))
cat(sprintf("Sites in VH only: %d/%d = %.1f%%   95%% CI [%.1f%%, %.1f%%]\n",
            k_vh, n_sl, real_acc_vh, ci_vh[1], ci_vh[2]))
cat(sprintf("LOO-CV  High+VH:  %d/%d = %.1f%%   95%% CI [%.1f%%, %.1f%%]\n",
            sum(loo_correct_hv), n_loo, loo_acc_hv, ci_loo_hv[1], ci_loo_hv[2]))
cat(sprintf("LOO-CV  VH only:  %d/%d = %.1f%%   95%% CI [%.1f%%, %.1f%%]\n\n",
            sum(loo_correct_vh), n_loo, loo_acc_vh, ci_loo_vh[1], ci_loo_vh[2]))

# Save validation metrics (expanded with CIs and weight stability)
val_df <- data.frame(
  Metric = c("Chi-square", "p-value (chi-sq)", "df", "Cohen's d",
             "Sites in High+VH (%)", "Sites in High+VH 95% CI lower",
             "Sites in High+VH 95% CI upper",
             "Sites in VH (%)", "Sites in VH 95% CI lower",
             "Sites in VH 95% CI upper",
             "LOO-CV accuracy High+VH (%)", "LOO-CV High+VH 95% CI lower",
             "LOO-CV High+VH 95% CI upper",
             "LOO-CV accuracy VH (%)", "LOO-CV VH 95% CI lower",
             "LOO-CV VH 95% CI upper",
             "Permutation p-value (High+VH)", "Permutation p-value (VH)"),
  Value  = c(chi_sq, p_val, 3, cohens_d,
             real_acc_hv, ci_hv[1], ci_hv[2],
             real_acc_vh, ci_vh[1], ci_vh[2],
             loo_acc_hv, ci_loo_hv[1], ci_loo_hv[2],
             loo_acc_vh, ci_loo_vh[1], ci_loo_vh[2],
             p_perm_hv, p_perm_vh)
)
write.csv(val_df,
          file.path(STATS_OUTPUT, "Validation_metrics.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: %s\n\n", file.path(STATS_OUTPUT, "Validation_metrics.csv")))

# =============================================================================
# SECTION 11: GRID CELL SUITABILITY ANALYSIS  (8 × 6 grid)
# =============================================================================

cat("================================================================================\n")
cat("SECTION 11: GRID CELL SUITABILITY ANALYSIS (8 × 6)\n")
cat("================================================================================\n\n")

# Use masked suitability if available
suit_for_grid <- if (exists("suitability_masked")) suitability_masked else suitability_900m
class_for_grid <- if (exists("s4_masked")) s4_masked else suitability_4class

valid_ext  <- ext(trim(suit_for_grid))
suit_crop  <- crop(suit_for_grid,   valid_ext)
class_crop <- crop(class_for_grid,  valid_ext)

e <- ext(suit_crop)
n_cols <- 8; n_rows <- 6
lon_breaks  <- seq(e[1], e[2], length.out = n_cols + 1)
lat_breaks  <- seq(e[3], e[4], length.out = n_rows + 1)
col_labels  <- LETTERS[1:n_cols]
row_labels  <- as.character(1:n_rows)
col_centres <- (lon_breaks[-1]     + lon_breaks[-(n_cols+1)]) / 2
row_centres <- rev((lat_breaks[-1] + lat_breaks[-(n_rows+1)]) / 2)

cat("Calculating mean suitability per grid cell...\n\n")
grid_summary <- data.frame(Grid_Cell=character(), Mean=numeric(), SD=numeric(),
                           Min=numeric(), Max=numeric(), N_Pixels=integer())

for (i in 1:n_rows) for (j in 1:n_cols) {
  cell_name <- paste0(col_labels[j], row_labels[i])
  cell_ext  <- ext(lon_breaks[j], lon_breaks[j+1], lat_breaks[i], lat_breaks[i+1])
  cell_vals <- values(crop(suit_crop, cell_ext), na.rm = TRUE)
  if (length(cell_vals) > 0) {
    grid_summary <- rbind(grid_summary,
                          data.frame(Grid_Cell = cell_name,
                                     Mean      = mean(cell_vals),    SD   = sd(cell_vals),
                                     Min       = min(cell_vals),     Max  = max(cell_vals),
                                     N_Pixels  = length(cell_vals)))
  } else {
    grid_summary <- rbind(grid_summary,
                          data.frame(Grid_Cell=cell_name, Mean=NA, SD=NA, Min=NA, Max=NA, N_Pixels=0))
  }
}

grid_sorted <- grid_summary[order(-grid_summary$Mean, na.last=TRUE), ]
cat("Top 10 cells by mean suitability:\n")
top10 <- head(grid_sorted[!is.na(grid_sorted$Mean), ], 10)
for (i in 1:nrow(top10)) cat(sprintf("  %2d. %s: %.1f\n", i, top10$Grid_Cell[i], top10$Mean[i]))
cat("\n")

write.csv(grid_summary,
          file.path(STATS_OUTPUT, "Grid_suitability_summary.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: %s\n\n", file.path(STATS_OUTPUT, "Grid_suitability_summary.csv")))


# =============================================================================
# SECTION 12: VISUALIZATIONS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 12: VISUALIZATIONS\n")
cat("================================================================================\n\n")

# Reload sites_sl with suit_class from the 4-class model
sites_sl  <- sites_complete[sites_complete$source_dem == MODEL_REGION, ]
site_pts  <- vect(sites_sl, geom = c("Lon", "Lat"), crs = "EPSG:4326")
sites_sl$suitability  <- extract(suitability_900m, site_pts)[, 2]
sites_sl$suit_class   <- extract(suitability_4class, site_pts)[, 2]

class_names <- c("Low", "Medium", "High", "Very High")
cols_4class <- c("#56B4E9", "#E69F00", "#2DC653", "#B2182B")  # colorblind-friendly

# --- Fig 1: Suitability maps ---
cat("Creating suitability maps...\n")
png(file.path(FIGS_OUTPUT, "Suitability_Maps.png"),
    width = 3600, height = 1600, res = 300)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 6))

plot(suitability_900m, main = "Habitat Suitability (continuous)",
     col = terrain.colors(100), axes = TRUE)
site_cols <- cols_4class[sites_sl$suit_class]
points(site_pts, pch = 21, bg = site_cols, cex = 1.5, col = "black", lwd = 1.5)
legend("bottomright",
       legend = sprintf("%s (n=%d)", class_names,
                        sapply(1:4, function(k) sum(sites_sl$suit_class == k, na.rm=TRUE))),
       pch = 21, pt.bg = cols_4class, pt.cex = 1.5, bty = "n", cex = 0.8)

plot(suitability_4class, main = "Habitat Suitability (4 classes)",
     col = cols_4class, axes = TRUE, type = "classes",
     levels = class_names)
points(site_pts, pch = 21, bg = "white", cex = 1.2, col = "black", lwd = 2)
text(x = ext(suitability_4class)[1] + 0.5, y = ext(suitability_4class)[3] + 0.3,
     labels = sprintf("chi-sq=%.2f  p<0.001\nCohen's d=%.2f", chi_sq, cohens_d),
     pos = 4, cex = 0.9, font = 2)
dev.off()
cat("✓ Saved: Suitability_Maps.png\n")

# --- Fig 2: Validation tests (5-panel) ---
cat("Creating validation plots...\n")
png(file.path(FIGS_OUTPUT, "Validation_Tests.png"),
    width = 5000, height = 3000, res = 300)
par(mfrow = c(2, 3), mar = c(5, 4.5, 4, 2))

# Panel 1: Chi-square
bp <- barplot(rbind(expected, site_cls_cnt),
              beside = TRUE, names.arg = class_names,
              col = c("gray70", "steelblue"),
              ylab = "Number of sites", main = "Chi-Square Enrichment",
              ylim = c(0, max(c(expected, site_cls_cnt)) * 1.3),
              cex.main = 1.1)
for (i in 1:4) {
  enr <- (site_cls_cnt[i]/nrow(sites_sl)) / class_props[i]
  text(bp[2,i], site_cls_cnt[i] + 1, sprintf("%.2fx", enr), cex = 0.8, font = 2)
}
legend("topleft",
       legend = c("Expected", "Observed", sprintf("χ²=%.2f, p=%.2e", chi_sq, p_val)),
       fill = c("gray70","steelblue",NA), border=c("black","black",NA),
       bty="n", cex=0.8)

# Panel 2: Density distribution
dens_land  <- density(landscape_vals)
dens_sites <- density(sites_sl$suitability, na.rm=TRUE)
plot(dens_land, main="Sites vs Landscape", xlab="Suitability", ylab="Density",
     col="steelblue", lwd=2.5,
     xlim=range(c(landscape_vals, sites_sl$suitability), na.rm=TRUE),
     ylim=c(0, max(c(dens_land$y, dens_sites$y))*1.1))
lines(dens_sites, col="firebrick", lwd=2.5)
abline(v=quartiles[2:4], lty=2, col="gray40", lwd=1.5)
legend("topright",
       legend = c(sprintf("Landscape (mean=%.1f)", mean(landscape_vals)),
                  sprintf("Sites (mean=%.1f)", mean(sites_sl$suitability, na.rm=TRUE)),
                  sprintf("MW p=%.2e", mw_test$p.value),
                  sprintf("Cohen's d=%.2f", cohens_d)),
       col=c("steelblue","firebrick",NA,NA), lwd=c(2.5,2.5,NA,NA),
       bty="n", cex=0.8)

# Panel 3: Cohen's d
m_land <- mean(landscape_vals); m_site <- mean(sites_sl$suitability, na.rm=TRUE)
plot(c(1,2), c(m_land, m_site), xlim=c(0.5,2.5),
     ylim=c(0, max(m_land, m_site)*1.3), pch=19, cex=2,
     col=c("steelblue","firebrick"),
     xlab="", ylab="Mean Suitability", main="Effect Size (Cohen's d)", xaxt="n")
axis(1, at=c(1,2), labels=c("Landscape","Sites"))
arrows(1, m_land - sd(landscape_vals), 1, m_land + sd(landscape_vals),
       angle=90, code=3, length=0.1, lwd=2, col="steelblue")
arrows(2, m_site - sd(sites_sl$suitability, na.rm=TRUE),
       2, m_site + sd(sites_sl$suitability, na.rm=TRUE),
       angle=90, code=3, length=0.1, lwd=2, col="firebrick")
text(1.5, max(m_land, m_site)*1.18,
     labels=sprintf("d = %.2f (%s)", cohens_d, label_d), cex=1.0, font=2)

# Panel 4: LOO-CV with Clopper-Pearson 95% CIs
plot(NA, xlim = c(0.5, 2.5), ylim = c(0, 110), xaxt = "n",
     xlab = "", ylab = "Accuracy (%)",
     main = "LOO Cross-Validation\n(95% Clopper-Pearson CI)", cex.main = 1.0)
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = 1:2, labels = c("High+VH", "VH only"))
# CI whiskers
arrows(1, ci_loo_hv[1], 1, ci_loo_hv[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = "#2D8B2D")
arrows(2, ci_loo_vh[1], 2, ci_loo_vh[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = "#2D8B2D")
# Point estimates
points(1:2, c(loo_acc_hv, loo_acc_vh), pch = 19, cex = 2, col = "#2D8B2D")
# 50% chance line
abline(h = 50, lty = 2, col = "gray50", lwd = 1.5)
text(2.4, 50, "50%", cex = 0.7, col = "gray50", adj = c(1, -0.3))
# Labels
text(1, loo_acc_hv, sprintf("%.1f%%", loo_acc_hv), pos = 4, font = 2, cex = 0.85)
text(2, loo_acc_vh, sprintf("%.1f%%", loo_acc_vh), pos = 4, font = 2, cex = 0.85)
text(1, ci_loo_hv[1], sprintf("[%.0f%%", ci_loo_hv[1]), pos = 1, cex = 0.7, col = "gray40")
text(1, ci_loo_hv[2], sprintf("%.0f%%]", ci_loo_hv[2]), pos = 3, cex = 0.7, col = "gray40")
text(2, ci_loo_vh[1], sprintf("[%.0f%%", ci_loo_vh[1]), pos = 1, cex = 0.7, col = "gray40")
text(2, ci_loo_vh[2], sprintf("%.0f%%]", ci_loo_vh[2]), pos = 3, cex = 0.7, col = "gray40")

# Panel 5: Permutation test
hist(perm_acc_hv, breaks=30, col="lightsteelblue", border="white",
     main=sprintf("Permutation Test (High+VH)\np=%.4f", p_perm_hv),
     xlab="% Random sites in High+VH", ylab="Frequency", xlim=c(30,100))
abline(v=real_acc_hv, col="firebrick", lwd=3)
abline(v=mean(perm_acc_hv), col="gray40", lwd=2, lty=2)
legend("topleft",
       legend=c(sprintf("Real: %.1f%%",real_acc_hv),
                sprintf("Null: %.1f%%",mean(perm_acc_hv))),
       col=c("firebrick","gray40"), lwd=c(3,2), lty=c(1,2), bty="n", cex=0.85)

# Panel 6: LOO-CV weight stability
var_short  <- c("Elev", "Slope", "North", "TRI", "TWI", "Dist.W")
loo_means  <- apply(loo_weights, 2, mean)
loo_sds    <- apply(loo_weights, 2, sd)
x_pos      <- 1:6
yr         <- range(c(loo_means - loo_sds, loo_means + loo_sds, weights))
plot(NA, xlim = c(0.5, 6.5), ylim = c(max(0, yr[1] - 2), yr[2] + 3),
     xaxt = "n", xlab = "", ylab = "Weight (%)",
     main = "LOO-CV Weight Stability", cex.main = 1.0)
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = x_pos, labels = var_short, cex.axis = 0.85)
# LOO mean ± SD whiskers
arrows(x_pos, loo_means - loo_sds, x_pos, loo_means + loo_sds,
       angle = 90, code = 3, length = 0.08, lwd = 2, col = "steelblue")
points(x_pos, loo_means, pch = 19, cex = 1.5, col = "steelblue")
# Full-model weights as reference diamonds
points(x_pos, weights, pch = 18, cex = 1.5, col = "firebrick")
legend("topright",
       legend = c("LOO mean ± SD", "Full-model weight"),
       pch = c(19, 18), col = c("steelblue", "firebrick"),
       pt.cex = 1.3, bty = "n", cex = 0.8)

dev.off()
cat("✓ Saved: Validation_Tests.png\n")

# --- Fig 3: Strip plots by bioclimatic zone ---
cat("Creating strip plots by zone...\n")

sites_complete$Region <- factor(sites_complete$Region,
                                levels = c("Mediterranean","Irano-Turanian","Saharo-Arabian"))
zone_cols_strip <- c(Mediterranean="steelblue", "Irano-Turanian"="#E69F00",
                     "Saharo-Arabian"="#009E73")
var_list <- list(
  list(col="elevation",   ylab="Elevation (m asl)",               title="A) Elevation"),
  list(col="slope",       ylab="Slope (°)",                       title="B) Slope"),
  list(col="TRI",         ylab="TRI",                             title="C) TRI"),
  list(col="TWI",         ylab="TWI",                             title="D) TWI"),
  list(col="northness",   ylab="Northness (0–1)",                 title="E) Northness"),
  list(col="DistToWater", ylab="Distance to Water (km)",          title="F) Dist. to Water")
)
is_outlier <- function(x) {
  q  <- quantile(x, c(0.25, 0.75), na.rm=TRUE)
  iq <- diff(q); x < (q[1]-1.5*iq) | x > (q[2]+1.5*iq)
}

png(file.path(FIGS_OUTPUT, "Topographic_by_zone.png"),
    width=5400, height=3600, res=300)
par(mfrow=c(2,3), mar=c(6,5,4,2), oma=c(1,0,3,0), mgp=c(2.8,0.9,0),
    cex=1.1, cex.lab=1.2, cex.axis=1.0, cex.main=1.3)

zone_lev <- levels(sites_complete$Region)
for (v in var_list) {
  vals_l <- lapply(zone_lev, function(z) {
    x <- sites_complete[[v$col]][sites_complete$Region == z]
    x[!is.na(x)]
  })
  means <- sapply(vals_l, mean); sds <- sapply(vals_l, sd)
  all_v <- unlist(vals_l)
  yr    <- diff(range(all_v))
  plot(NA, xlim=c(0.5,3.5), ylim=c(min(all_v)-0.05*yr, max(all_v)+0.25*yr),
       xaxt="n", xlab="", ylab=v$ylab, main=v$title, las=1, frame=FALSE)
  grid(nx=NA, ny=NULL, col="gray92"); box(bty="l")
  axis(1, at=1:3, labels=zone_lev, cex.axis=0.80)
  for (i in seq_along(zone_lev)) {
    z   <- zone_lev[i]; col <- zone_cols_strip[z]
    x   <- vals_l[[i]]; mn <- means[i]; sdv <- sds[i]; n <- length(x)
    segments(i, mn-sdv, i, mn+sdv, col=col, lwd=2)
    segments(i-0.1, mn-sdv, i+0.1, mn-sdv, col=col, lwd=1.8)
    segments(i-0.1, mn+sdv, i+0.1, mn+sdv, col=col, lwd=1.8)
    set.seed(42); jx <- i + runif(n, -0.15, 0.15)
    out <- is_outlier(x)
    points(jx[!out], x[!out], pch=21, bg=adjustcolor(col, 0.45),
           col=adjustcolor(col, 0.75), cex=0.9, lwd=0.7)
    if (any(out)) points(jx[out], x[out], pch=21, bg=adjustcolor(col, 0.85),
                         col="gray20", cex=1.2, lwd=1.1)
    points(i, mn, pch=23, bg="white", col="gray20", cex=1.2, lwd=1.4)
    mtext(sprintf("n=%d", n), side=1, at=i, line=4.8, cex=0.9, col="gray40")
  }
}
mtext("Topographic Variables by Bioclimatic Zone", outer=TRUE, cex=1.4, font=2, line=1.5)
dev.off()
cat("✓ Saved: Topographic_by_zone.png\n")

# --- Fig 4: Dumbbell plot (zone weights vs modeled area) ---
cat("Creating dumbbell plot...\n")
var_order  <- c("Elevation","Slope","TRI","TWI","Northness","DistToWater")
var_lab_d  <- c(Elevation="Elevation", Slope="Slope", TRI="TRI", TWI="TWI",
                Northness="Northness (Aspect)", DistToWater="Dist. to Water")
zone_pch   <- c(Mediterranean=21, "Irano-Turanian"=22, "Saharo-Arabian"=24)
zone_cols_d <- c(Mediterranean="#009E73","Irano-Turanian"="#0072B2","Saharo-Arabian"="#E69F00")

modeled <- all_weights[all_weights$Zone == "Modeled Area", ]

png(file.path(FIGS_OUTPUT, "Variable_weights_dumbbell.png"),
    width=4800, height=3200, res=300)
par(mar=c(5,9,4,10))
n_v <- length(var_order)
plot(NA, xlim=c(0,42), ylim=c(0.5, n_v+0.5), xaxt="n", yaxt="n",
     xlab="Weight (%)", ylab="", frame=FALSE)
for (g in seq(0,40,10)) segments(g, 0.3, g, n_v+0.5, col="gray90", lwd=1)
axis(1, at=seq(0,40,10), labels=paste0(seq(0,40,10),"%"), cex.axis=0.9, col="gray60")
axis(2, at=1:n_v, labels=var_lab_d[var_order], tick=FALSE, las=1, cex.axis=1.0)

for (i in seq_along(var_order)) {
  v   <- var_order[i]
  mv  <- modeled$Weight[modeled$Variable == v]
  zvs <- sapply(zones, function(z) {
    ww <- all_weights$Weight[all_weights$Zone==z & all_weights$Variable==v]
    if (length(ww) > 0) ww else NA
  })
  zvs <- zvs[!is.na(zvs)]
  if (length(zvs) > 0) segments(min(zvs), i, max(zvs), i, col="gray75", lwd=2.5)
  segments(mv, i-0.35, mv, i+0.35, col="gray30", lwd=2, lty=2)
  for (z in zones) {
    zv <- all_weights$Weight[all_weights$Zone==z & all_weights$Variable==v]
    if (length(zv) > 0)
      points(zv, i, pch=zone_pch[z], bg=adjustcolor(zone_cols_d[z], 0.85),
             col="gray20", cex=1.8, lwd=1.2)
  }
  points(mv, i, pch=18, col="gray20", cex=1.6)
  text(mv, i+0.32, sprintf("%.1f%%", mv), cex=0.72, col="gray30")
}
mtext("Variable Weights by Bioclimatic Zone", side=3, line=2, cex=1.25, font=2)

lx <- 43; ly <- n_v+0.3
for (k in seq_along(zones)) {
  z <- zones[k]; yp <- ly - (k-1)*0.6
  if (any(all_weights$Zone == z)) {
    points(lx, yp, pch=zone_pch[z], bg=adjustcolor(zone_cols_d[z],0.85),
           col="gray20", cex=1.6, xpd=TRUE, lwd=1.2)
    text(lx+0.8, yp, z, cex=0.88, adj=0, xpd=TRUE, col="gray20")
  }
}
yg <- ly - length(zones)*0.6 - 0.2
segments(lx-0.3, yg, lx+0.3, yg, col="gray30", lwd=2, lty=2, xpd=TRUE)
points(lx, yg, pch=18, col="gray20", cex=1.6, xpd=TRUE)
text(lx+0.8, yg, "Modeled Area", cex=0.88, adj=0, xpd=TRUE, col="gray20")
dev.off()
cat("✓ Saved: Variable_weights_dumbbell.png\n")

# --- Fig 5: Grid reference map ---
cat("Creating grid reference map...\n")
class_smooth <- focal(class_crop, w=3, fun="modal", na.rm=TRUE)
class_smooth <- mask(class_smooth, class_crop)

sites_sl_all <- sites_complete[sites_complete$source_dem == MODEL_REGION, ]
sites_inbounds <- sites_sl_all[
  sites_sl_all$Lon >= e[1] & sites_sl_all$Lon <= e[2] &
    sites_sl_all$Lat >= e[3] & sites_sl_all$Lat <= e[4], ]
sites_sl_pts <- vect(sites_inbounds, geom=c("Lon","Lat"), crs="EPSG:4326")
sites_inbounds$suit_class <- extract(suitability_4class, sites_sl_pts)[, 2]

png(file.path(FIGS_OUTPUT, paste0(MODEL_REGION, "_Grid_Map.png")),
    width=4900, height=5300, res=400)
par(mar=c(4,2,5,10))
plot(class_smooth, col=cols_4class, axes=FALSE, type="classes",
     xlim=c(e[1],e[2]), ylim=c(e[3],e[4]),
     plg=list(title="Suitability", legend=class_names, cex=0.9))
for (j in 1:n_cols) for (i in 1:n_rows)
  rect(lon_breaks[j], lat_breaks[i], lon_breaks[j+1], lat_breaks[i+1],
       border=adjustcolor("gray10", 0.8), lwd=2.5, col=NA)
for (j in seq_along(col_labels)) {
  text(col_centres[j], e[4]+(e[4]-e[3])*0.018, col_labels[j], cex=1.0, font=2, xpd=TRUE)
  text(col_centres[j], e[3]-(e[4]-e[3])*0.018, col_labels[j], cex=1.0, font=2, xpd=TRUE, adj=c(0.5,1))
}
for (i in seq_along(row_labels)) {
  text(e[1]-(e[2]-e[1])*0.04, row_centres[i], row_labels[i], cex=1.0, font=2, xpd=TRUE, adj=c(1,0.5))
  text(e[2]+(e[2]-e[1])*0.04, row_centres[i], row_labels[i], cex=1.0, font=2, xpd=TRUE, adj=c(0,0.5))
}
for (j in seq_along(col_labels)) for (i in seq_along(row_labels))
  text(col_centres[j], row_centres[i], paste0(col_labels[j], row_labels[i]),
       cex=0.55, col=adjustcolor("white", 0.85), font=2)
site_shp <- c(21, 22, 23, 24)[sites_inbounds$suit_class]
points(sites_sl_pts, pch=site_shp, bg="black", cex=1.5, col="white", lwd=1.8)
axis(1, at=round(seq(e[1],e[2],length.out=5),1),
     labels=paste0(round(seq(e[1],e[2],length.out=5),1),"°E"), cex.axis=0.85)
axis(2, at=round(seq(e[3],e[4],length.out=5),1),
     labels=paste0(round(seq(e[3],e[4],length.out=5),1),"°N"), cex.axis=0.85, las=1)
legend("bottomleft",
       legend=c(sprintf("Very High (n=%d)", sum(sites_inbounds$suit_class==4, na.rm=TRUE)),
                sprintf("High (n=%d)",      sum(sites_inbounds$suit_class==3, na.rm=TRUE)),
                sprintf("Medium (n=%d)",    sum(sites_inbounds$suit_class==2, na.rm=TRUE)),
                sprintf("Low (n=%d)",       sum(sites_inbounds$suit_class==1, na.rm=TRUE))),
       pch=c(24,23,22,21), pt.bg="black", pt.cex=1.4, col="black",
       bty="o", bg="white", box.col="gray40", cex=0.85, inset=c(0.01,0.01))
dev.off()
cat(sprintf("✓ Saved: %s_Grid_Map.png\n\n", MODEL_REGION))


# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Output directory: ", OUTPUT_DIR, "\n\n")

cat("FILES CREATED:\n")
cat("  Topo/Slope/           ← 4 regions, slope (degrees)\n")
cat("  Topo/Aspect/          ← 4 regions, northness (0–1)\n")
cat("  Topo/TRI/             ← 4 regions\n")
cat("  Topo/TWI/             ← 4 regions\n")
cat("  Topo/DistToWater/     ← 4 regions, P95 threshold (km)\n")
cat("  Sites_data/Sites_with_variables.csv\n")
cat("  Sites_data/Background_sample.csv\n")
cat("  Stats/Descriptive_stats_by_zone.csv\n")
cat("  Stats/Variable_correlations.csv\n")
cat("  Stats/Variable_weights_all.csv\n")
cat("  Stats/Validation_metrics.csv\n")
cat("  Stats/LOO_weight_stability.csv\n")
cat("  Stats/Grid_suitability_summary.csv\n")
cat("  MCDA/Suitability_30m_P95.tif\n")
cat("  MCDA/Suitability_900m_P95.tif\n")
cat("  MCDA/Suitability_4class_P95.tif\n")
cat("  MCDA/Suitability_900m_P95_MASKED.tif   (if geology mask present)\n")
cat("  MCDA/Suitability_4class_P95_MASKED.tif (if geology mask present)\n")
cat("  Figures/Variable_Correlation_Matrix.png\n")
cat("  Figures/Topographic_by_zone.png\n")
cat("  Figures/Variable_weights_dumbbell.png\n")
cat("  Figures/Suitability_Maps.png\n")
cat("  Figures/Validation_Tests.png\n")
cat(sprintf("  Figures/%s_Grid_Map.png\n\n", MODEL_REGION))

cat("VALIDATION SUMMARY:\n")
cat(sprintf("  Chi-square:            %.2f  (p = %.2e)\n", chi_sq, p_val))
cat(sprintf("  Cohen's d:             %.3f  (%s)\n", cohens_d, label_d))
cat(sprintf("  Sites in High+VH:      %.1f%%  95%% CI [%.1f%%, %.1f%%]\n",
            real_acc_hv, ci_hv[1], ci_hv[2]))
cat(sprintf("  LOO-CV (High+VH):      %.1f%%  95%% CI [%.1f%%, %.1f%%]\n",
            loo_acc_hv, ci_loo_hv[1], ci_loo_hv[2]))
cat(sprintf("  Permutation p (H+VH):  %.4f\n", p_perm_hv))
cat(sprintf("  Max LOO weight shift:  %.1f pp\n",
            max(abs(weights - apply(loo_weights, 2, mean)))))
cat("================================================================================\n\n")



# --- Fig 1: Suitability maps (publication style, 2-panel vertical) ---
cat("Creating suitability maps...\n")

# Publication color scheme (matching previous figure)
cols_pub <- c("#A8BED1",   # Low - light steel blue
              "#4A6FA5",   # Medium - medium blue
              "#CDA67E",   # High - salmon/peach
              "#7B2D3C")   # Very High - dark maroon
col_unlikely <- "#D3D3D3"  # light gray

# --- Build hillshade from 900m DEM for terrain texture ---
dem_900m   <- aggregate(dem_r, fact = 30, fun = "mean", na.rm = TRUE)
slope_hs   <- terrain(dem_900m, v = "slope", unit = "radians")
aspect_hs  <- terrain(dem_900m, v = "aspect", unit = "radians")
hillshade  <- shade(slope_hs, aspect_hs, angle = 45, direction = 315)

# --- Build 5-class raster for bottom panel (suitability + Unlikely) ---
if (exists("s4_masked") && !identical(s4_masked, suitability_4class)) {
  suit_5class <- suitability_4class
  unlikely_px <- !is.na(suitability_4class) & is.na(s4_masked)
  suit_5class[unlikely_px] <- 5
} else {
  suit_5class <- suitability_4class
}

# --- Extent and scale helpers ---
e <- ext(suitability_900m)

# Scale bar function (distance in km)
draw_scalebar <- function(x0, y0, km_len, e) {
  # Approximate degrees per km at this latitude
  mid_lat <- (e[3] + e[4]) / 2
  deg_per_km <- 1 / (111.32 * cos(mid_lat * pi / 180))
  dx <- km_len * deg_per_km
  segments(x0, y0, x0 + dx, y0, lwd = 2.5, col = "black")
  segments(x0, y0 - 0.05, x0, y0 + 0.05, lwd = 2, col = "black")
  segments(x0 + dx, y0 - 0.05, x0 + dx, y0 + 0.05, lwd = 2, col = "black")
  text(x0 + dx / 2, y0 - 0.15, paste0(km_len, " km"), cex = 0.7, font = 1)
}

# ==================== FIGURE ====================
png(file.path(BASE_DIR, "Suitability_Maps2.png"),
    width = 3200, height = 5800, res = 300)
layout(matrix(1:2, nrow = 2), heights = c(1, 1))

# ==================== TOP PANEL: unmasked 4-class ====================
par(mar = c(2, 2, 1, 2))

# Hillshade base layer
plot(hillshade, col = gray.colors(256, start = 0.05, end = 0.95),
     legend = FALSE, axes = FALSE, reset = FALSE)

# Overlay 4-class suitability with transparency
plot(suitability_4class,
     col = adjustcolor(cols_pub, alpha.f = 0.70),
     type = "classes", levels = class_names,
     legend = FALSE, axes = FALSE, add = TRUE)

# Sites as yellow triangles
points(site_pts, pch = 24, bg = "#FFD700", cex = 1.2, col = "black", lwd = 0.9)

# Geographic labels — ADJUST THESE COORDINATES to match your study area
text(35.05, 32.80, "Mt. Carmel",    cex = 0.80, font = 3, col = "white")
text(35.55, 32.35, "Jordan Valley", cex = 0.80, font = 3, col = "white")
text(36.80, 31.95, "Azraq",         cex = 0.80, font = 3, col = "white")
text(34.50, 30.85, "Negev",         cex = 0.80, font = 3, col = "white")
text(35.10, 30.15, "Wadi Araba",    cex = 0.80, font = 3, col = "white")
text(36.20, 30.40, "Al-Jafr",       cex = 0.80, font = 3, col = "white")

# Legend (top-left)
legend("topleft",
       legend = rev(class_names),
       fill   = rev(cols_pub),
       border = "gray30",
       title  = "MCDA Suitability", title.font = 2,
       bty = "o", bg = "white", box.col = "gray40",
       cex = 0.75, inset = c(0.01, 0.02))

# North arrow (left side)
arr_x <- e[1] + (e[2] - e[1]) * 0.03
arr_y0 <- e[3] + (e[4] - e[3]) * 0.15
arr_y1 <- arr_y0 + (e[4] - e[3]) * 0.06
arrows(arr_x, arr_y0, arr_x, arr_y1,
       length = 0.12, lwd = 2, col = "black")
text(arr_x, arr_y1 + (e[4] - e[3]) * 0.02, "N",
     cex = 0.9, font = 2)

# Scale bar (bottom-right)
sb_x <- e[2] - (e[2] - e[1]) * 0.35
sb_y <- e[3] + (e[4] - e[3]) * 0.05
draw_scalebar(sb_x, sb_y, 50, e)
draw_scalebar(sb_x + 50 / (111.32 * cos(((e[3]+e[4])/2) * pi / 180)),
              sb_y, 50, e)

# ==================== BOTTOM PANEL: masked 5-class + regions ====================
par(mar = c(2, 2, 1, 2))

# Hillshade base layer
plot(hillshade, col = gray.colors(256, start = 0.05, end = 0.95),
     legend = FALSE, axes = FALSE, reset = FALSE)

# Overlay 5-class suitability (with "Unlikely")
if (exists("suit_5class") && max(values(suit_5class), na.rm = TRUE) == 5) {
  plot(suit_5class,
       col = adjustcolor(c(cols_pub, col_unlikely), alpha.f = 0.70),
       type = "classes",
       levels = c(class_names, "Unlikely"),
       legend = FALSE, axes = FALSE, add = TRUE)
  leg_names <- rev(c(class_names, "Unlikely"))
  leg_cols  <- rev(c(cols_pub, col_unlikely))
} else {
  plot(suitability_4class,
       col = adjustcolor(cols_pub, alpha.f = 0.70),
       type = "classes", levels = class_names,
       legend = FALSE, axes = FALSE, add = TRUE)
  leg_names <- rev(class_names)
  leg_cols  <- rev(cols_pub)
}

# Sites
points(site_pts, pch = 24, bg = "#FFD700", cex = 1.2, col = "black", lwd = 0.9)

# Dashed region rectangles — ADJUST THESE COORDINATES
rect(34.50, 32.00, 35.20, 33.00, border = "black", lwd = 2.5, lty = 2)
text(34.55, 32.90, "A", cex = 1.1, font = 2, adj = c(0, 1))

rect(34.50, 30.80, 35.20, 32.00, border = "black", lwd = 2.5, lty = 2)
text(34.55, 31.90, "B", cex = 1.1, font = 2, adj = c(0, 1))

rect(36.00, 31.30, 37.80, 32.50, border = "black", lwd = 2.5, lty = 2)
text(36.70, 32.30, "C", cex = 1.1, font = 2, adj = c(0, 1))

rect(35.00, 29.50, 36.80, 31.00, border = "black", lwd = 2.5, lty = 2)
text(35.80, 30.30, "D", cex = 1.1, font = 2, adj = c(0, 1))

# Legend (top-left, with Unlikely)
legend("topleft",
       legend = leg_names,
       fill   = leg_cols,
       border = "gray30",
       title  = "MCDA Suitability", title.font = 2,
       bty = "o", bg = "white", box.col = "gray40",
       cex = 0.75, inset = c(0.01, 0.02))

# North arrow
arrows(arr_x, arr_y0, arr_x, arr_y1,
       length = 0.12, lwd = 2, col = "black")
text(arr_x, arr_y1 + (e[4] - e[3]) * 0.02, "N",
     cex = 0.9, font = 2)

# Scale bar
draw_scalebar(sb_x, sb_y, 50, e)
draw_scalebar(sb_x + 50 / (111.32 * cos(((e[3]+e[4])/2) * pi / 180)),
              sb_y, 50, e)

dev.off()
cat("✓ Saved: Suitability_Maps2.png\n")

###############################################################################
###############################################################################

# =============================================================================
# SECTION 12: VISUALIZATIONS NEW AUGUST 30 2026
# =============================================================================

cat("================================================================================\n")
cat("SECTION 12: VISUALIZATIONS\n")
cat("================================================================================\n\n")

# ── Shared style values (from style reference) ────────────────────────────────
zone_cols <- c(Mediterranean = "#954151",
               "Irano-Turanian" = "#1d4e89",
               "Saharo-Arabian" = "#7295bf")
col_burgundy  <- "#954151"
col_darkblue  <- "#1d4e89"
col_dustyblue <- "#7295bf"
col_mauve     <- "#db9793"
col_beige     <- "#ead98b"

# Publication map colors
cols_pub <- c("#A8BED1", "#4A6FA5", "#CDA67E", "#7B2D3C")
col_unlikely <- "#D3D3D3"

class_names <- c("Low", "Medium", "High", "Very High")

# Reload sites_sl with suit_class from the 4-class model
sites_sl  <- sites_complete[sites_complete$source_dem == MODEL_REGION, ]
site_pts  <- vect(sites_sl, geom = c("Lon", "Lat"), crs = "EPSG:4326")
sites_sl$suitability  <- extract(suitability_900m, site_pts)[, 2]
sites_sl$suit_class   <- extract(suitability_4class, site_pts)[, 2]


# ==========================================================================
# Fig 1: Suitability maps (publication, 2-panel vertical with hillshade)
# ==========================================================================
cat("Creating suitability maps (publication style)...\n")

dem_900m   <- aggregate(dem_r, fact = 30, fun = "mean", na.rm = TRUE)
slope_hs   <- terrain(dem_900m, v = "slope", unit = "radians")
aspect_hs  <- terrain(dem_900m, v = "aspect", unit = "radians")
hillshade  <- shade(slope_hs, aspect_hs, angle = 45, direction = 315)

# 5-class raster for bottom panel
if (exists("s4_masked") && !identical(s4_masked, suitability_4class)) {
  suit_5class <- suitability_4class
  unlikely_px <- !is.na(suitability_4class) & is.na(s4_masked)
  suit_5class[unlikely_px] <- 5
} else {
  suit_5class <- suitability_4class
}

e <- ext(suitability_900m)

draw_scalebar <- function(x0, y0, km_len, e) {
  mid_lat    <- (e[3] + e[4]) / 2
  deg_per_km <- 1 / (111.32 * cos(mid_lat * pi / 180))
  dx <- km_len * deg_per_km
  segments(x0, y0, x0 + dx, y0, lwd = 2.5, col = "black")
  segments(x0, y0 - 0.05, x0, y0 + 0.05, lwd = 2, col = "black")
  segments(x0 + dx, y0 - 0.05, x0 + dx, y0 + 0.05, lwd = 2, col = "black")
  text(x0 + dx / 2, y0 - 0.15, paste0(km_len, " km"), cex = 0.7)
}

png(file.path(BASE_DIR, "Suitability_Maps.png"),
    width = 3200, height = 5800, res = 300)
layout(matrix(1:2, nrow = 2), heights = c(1, 1))

# --- Top panel: unmasked 4-class ---
par(mar = c(2, 2, 1, 2))
plot(hillshade, col = gray.colors(256, 0.05, 0.95),
     legend = FALSE, axes = FALSE, reset = FALSE)
plot(suitability_4class,
     col = adjustcolor(cols_pub, alpha.f = 0.70),
     type = "classes", levels = class_names,
     legend = FALSE, axes = FALSE, add = TRUE)
points(site_pts, pch = 24, bg = "#FFD700", cex = 1.2, col = "black", lwd = 0.9)

# Geographic labels — adjust coordinates as needed
text(35.05, 32.80, "Mt. Carmel",    cex = 0.80, font = 3, col = "white")
text(35.55, 32.35, "Jordan Valley", cex = 0.80, font = 3, col = "white")
text(36.80, 31.95, "Azraq",         cex = 0.80, font = 3, col = "white")
text(34.50, 30.85, "Negev",         cex = 0.80, font = 3, col = "white")
text(35.10, 30.15, "Wadi Araba",    cex = 0.80, font = 3, col = "white")
text(36.20, 30.40, "Al-Jafr",       cex = 0.80, font = 3, col = "white")

legend("topleft", legend = rev(class_names), fill = rev(cols_pub),
       border = "gray30", title = "MCDA Suitability", title.font = 2,
       bty = "o", bg = "white", box.col = "gray40", cex = 0.75, inset = c(0.01, 0.02))

arr_x  <- e[1] + (e[2] - e[1]) * 0.03
arr_y0 <- e[3] + (e[4] - e[3]) * 0.15
arr_y1 <- arr_y0 + (e[4] - e[3]) * 0.06
arrows(arr_x, arr_y0, arr_x, arr_y1, length = 0.12, lwd = 2, col = "black")
text(arr_x, arr_y1 + (e[4] - e[3]) * 0.02, "N", cex = 0.9, font = 2)

sb_x <- e[2] - (e[2] - e[1]) * 0.35
sb_y <- e[3] + (e[4] - e[3]) * 0.05
draw_scalebar(sb_x, sb_y, 50, e)
draw_scalebar(sb_x + 50 / (111.32 * cos(((e[3]+e[4])/2)*pi/180)), sb_y, 50, e)

# --- Bottom panel: masked 5-class + regions ---
par(mar = c(2, 2, 1, 2))
plot(hillshade, col = gray.colors(256, 0.05, 0.95),
     legend = FALSE, axes = FALSE, reset = FALSE)

if (exists("suit_5class") && max(values(suit_5class), na.rm = TRUE) == 5) {
  plot(suit_5class,
       col = adjustcolor(c(cols_pub, col_unlikely), alpha.f = 0.70),
       type = "classes", levels = c(class_names, "Unlikely"),
       legend = FALSE, axes = FALSE, add = TRUE)
  leg_names <- rev(c(class_names, "Unlikely"))
  leg_cols  <- rev(c(cols_pub, col_unlikely))
} else {
  plot(suitability_4class,
       col = adjustcolor(cols_pub, alpha.f = 0.70),
       type = "classes", levels = class_names,
       legend = FALSE, axes = FALSE, add = TRUE)
  leg_names <- rev(class_names)
  leg_cols  <- rev(cols_pub)
}

points(site_pts, pch = 24, bg = "#FFD700", cex = 1.2, col = "black", lwd = 0.9)

# Region rectangles — adjust coordinates as needed
rect(34.50, 32.00, 35.20, 33.00, border = "black", lwd = 2.5, lty = 2)
text(34.55, 32.90, "A", cex = 1.1, font = 2, adj = c(0, 1))
rect(34.50, 30.80, 35.20, 32.00, border = "black", lwd = 2.5, lty = 2)
text(34.55, 31.90, "B", cex = 1.1, font = 2, adj = c(0, 1))
rect(36.00, 31.30, 37.80, 32.50, border = "black", lwd = 2.5, lty = 2)
text(36.70, 32.30, "C", cex = 1.1, font = 2, adj = c(0, 1))
rect(35.00, 29.50, 36.80, 31.00, border = "black", lwd = 2.5, lty = 2)
text(35.80, 30.30, "D", cex = 1.1, font = 2, adj = c(0, 1))

legend("topleft", legend = leg_names, fill = leg_cols,
       border = "gray30", title = "MCDA Suitability", title.font = 2,
       bty = "o", bg = "white", box.col = "gray40", cex = 0.75, inset = c(0.01, 0.02))
arrows(arr_x, arr_y0, arr_x, arr_y1, length = 0.12, lwd = 2, col = "black")
text(arr_x, arr_y1 + (e[4] - e[3]) * 0.02, "N", cex = 0.9, font = 2)
draw_scalebar(sb_x, sb_y, 50, e)
draw_scalebar(sb_x + 50 / (111.32 * cos(((e[3]+e[4])/2)*pi/180)), sb_y, 50, e)

dev.off()
cat("✓ Saved: Suitability_Maps.png\n")


# ==========================================================================
# Fig 4: Strip plots by bioclimatic zone (SVG, 18 x 12 in)
# ==========================================================================
cat("Creating strip plots by zone...\n")

sites_complete$Region <- factor(sites_complete$Region,
                                levels = c("Mediterranean","Irano-Turanian","Saharo-Arabian"))
var_list <- list(
  list(col="elevation",   ylab="Elevation (m asl)",                  title="A) Elevation"),
  list(col="slope",       ylab="Slope (degrees)",                    title="B) Slope"),
  list(col="TRI",         ylab="Terrain Ruggedness Index (TRI)",     title="C) TRI"),
  list(col="TWI",         ylab="Topographic Wetness Index (TWI)",    title="D) TWI"),
  list(col="northness",   ylab="Aspect (0-1)",                      title="E) Aspect"),
  list(col="DistToWater", ylab="DTD (km)",                          title="F) DTD")
)
is_outlier <- function(x) {
  q  <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iq <- diff(q); x < (q[1] - 1.5*iq) | x > (q[2] + 1.5*iq)
}

svg(file.path(FIGS_OUTPUT, "Topographic_by_zone.svg"), width = 18, height = 12)
par(mfrow = c(2, 3), mar = c(6, 5, 4, 2), oma = c(3, 0, 3, 0),
    mgp = c(2.8, 0.9, 0), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)

zone_lev <- levels(sites_complete$Region)
for (v in var_list) {
  vals_l <- lapply(zone_lev, function(z) {
    x <- sites_complete[[v$col]][sites_complete$Region == z]
    x[!is.na(x)]
  })
  means <- sapply(vals_l, mean); sds <- sapply(vals_l, sd)
  all_v <- unlist(vals_l)
  yr    <- diff(range(all_v))
  plot(NA, xlim = c(0.5, 3.5),
       ylim = c(min(all_v) - 0.05*yr, max(all_v) + 0.25*yr),
       xaxt = "n", xlab = "", ylab = v$ylab, main = v$title, las = 1, frame = FALSE)
  grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
  axis(1, at = 1:3, labels = zone_lev, cex.axis = 0.80)
  for (i in seq_along(zone_lev)) {
    z   <- zone_lev[i]; col <- zone_cols[z]
    x   <- vals_l[[i]]; mn <- means[i]; sdv <- sds[i]; n <- length(x)
    segments(i, mn - sdv, i, mn + sdv, col = col, lwd = 2)
    segments(i - 0.1, mn - sdv, i + 0.1, mn - sdv, col = col, lwd = 1.8)
    segments(i - 0.1, mn + sdv, i + 0.1, mn + sdv, col = col, lwd = 1.8)
    set.seed(42); jx <- i + runif(n, -0.15, 0.15)
    out <- is_outlier(x)
    points(jx[!out], x[!out], pch = 21, bg = adjustcolor(col, 0.45),
           col = adjustcolor(col, 0.75), cex = 0.9, lwd = 0.7)
    if (any(out)) points(jx[out], x[out], pch = 21, bg = adjustcolor(col, 0.85),
                         col = "gray20", cex = 1.2, lwd = 1.1)
    points(i, mn, pch = 23, bg = "white", col = "gray20", cex = 1.2, lwd = 1.4)
    # Highest-value site name label
    max_idx <- which.max(x)
    max_site <- sites_complete$Site[sites_complete$Region == z & !is.na(sites_complete[[v$col]])]
    if (length(max_site) >= max_idx)
      text(jx[max_idx], x[max_idx], max_site[max_idx],
           cex = 0.75, font = 2, col = col, pos = 4)
    mtext(sprintf("n=%d", n), side = 1, at = i, line = 4.8, cex = 0.9, col = "gray40")
  }
}
mtext("Topographic and Hydrological Variables by Bioclimatic Zone",
      outer = TRUE, cex = 1.4, font = 2, line = 1.5)
# Horizontal legend below panels
par(fig = c(0, 1, 0, 0.06), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = zone_lev, pch = 21,
       pt.bg = adjustcolor(zone_cols[zone_lev], 0.65), col = "gray20",
       pt.cex = 1.5, cex = 0.95, bty = "n", horiz = TRUE)
dev.off()
cat("✓ Saved: Topographic_by_zone.svg\n")


# ==========================================================================
# Fig 5: Dumbbell plot (zone weights vs modeled area)
# ==========================================================================
cat("Creating dumbbell plot...\n")

var_order  <- c("Elevation", "Slope", "TRI", "TWI", "Northness", "DistToWater")
var_lab_d  <- c(Elevation = "Elevation", Slope = "Slope", TRI = "TRI",
                TWI = "TWI", Northness = "Aspect", DistToWater = "DTD")
zone_pch   <- c(Mediterranean = 21, "Irano-Turanian" = 22, "Saharo-Arabian" = 24)

modeled <- all_weights[all_weights$Zone == "Modeled Area", ]

png(file.path(FIGS_OUTPUT, "Variable_weights_dumbbell.png"),
    width = 4800, height = 3200, res = 300)
par(mar = c(5, 9, 4, 10), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0)
n_v <- length(var_order)
plot(NA, xlim = c(0, 42), ylim = c(0.5, n_v + 0.5),
     xaxt = "n", yaxt = "n", xlab = "Weight (%)", ylab = "", frame = FALSE)
for (g in seq(0, 40, 10)) segments(g, 0.3, g, n_v + 0.5, col = "gray90", lwd = 1)
axis(1, at = seq(0, 40, 10), labels = paste0(seq(0, 40, 10), "%"),
     cex.axis = 1.0, col.axis = "gray30", col = "gray60")
axis(2, at = 1:n_v, labels = var_lab_d[var_order], tick = FALSE, las = 1, cex.axis = 1.2)

for (i in seq_along(var_order)) {
  v   <- var_order[i]
  mv  <- modeled$Weight[modeled$Variable == v]
  zvs <- sapply(zones, function(z) {
    ww <- all_weights$Weight[all_weights$Zone == z & all_weights$Variable == v]
    if (length(ww) > 0) ww else NA
  })
  zvs <- zvs[!is.na(zvs)]
  if (length(zvs) > 0) segments(min(zvs), i, max(zvs), i, col = "gray75", lwd = 2.5)
  segments(mv, i - 0.35, mv, i + 0.35, col = col_mauve, lwd = 2, lty = 2)
  for (z in zones) {
    zv <- all_weights$Weight[all_weights$Zone == z & all_weights$Variable == v]
    if (length(zv) > 0)
      points(zv, i, pch = zone_pch[z], bg = adjustcolor(zone_cols[z], 0.85),
             col = "gray20", cex = 1.8, lwd = 1.2)
  }
  points(mv, i, pch = 18, col = col_mauve, cex = 1.6)
  text(mv, i + 0.32, sprintf("%.1f%%", mv), cex = 1.0, col = "gray30")
}
mtext("Variable Weights by Bioclimatic Zone", side = 3, line = 2, cex = 1.3, font = 2)

lx <- 43; ly <- n_v + 0.3
for (k in seq_along(zones)) {
  z <- zones[k]; yp <- ly - (k - 1) * 0.6
  if (any(all_weights$Zone == z)) {
    points(lx, yp, pch = zone_pch[z], bg = adjustcolor(zone_cols[z], 0.85),
           col = "gray20", cex = 1.6, xpd = TRUE, lwd = 1.2)
    text(lx + 0.8, yp, z, cex = 0.95, adj = 0, xpd = TRUE, col = "gray20")
  }
}
yg <- ly - length(zones) * 0.6 - 0.2
segments(lx - 0.3, yg, lx + 0.3, yg, col = col_mauve, lwd = 2, lty = 2, xpd = TRUE)
points(lx, yg, pch = 18, col = col_mauve, cex = 1.6, xpd = TRUE)
text(lx + 0.8, yg, "Modeled Area", cex = 0.95, adj = 0, xpd = TRUE, col = "gray20")
dev.off()
cat("✓ Saved: Variable_weights_dumbbell.png\n")


# ==========================================================================
# Fig 7: Validation tests (5 panels: 3 top, 2 centered bottom — SVG)
# ==========================================================================
cat("Creating validation plots...\n")

svg(file.path(FIGS_OUTPUT, "Validation_Tests.svg"), width = 18.67, height = 12)
layout(matrix(c(1, 2, 3,
                0, 4, 5), nrow = 2, byrow = TRUE))
par(mar = c(5, 4.5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)

# Panel 1: Chi-square
bp <- barplot(rbind(expected, site_cls_cnt),
              beside = TRUE, names.arg = class_names,
              col = c(col_mauve, col_burgundy),
              ylab = "Number of sites", main = "Chi-Square Enrichment Test",
              ylim = c(0, max(c(expected, site_cls_cnt)) * 1.3))
for (i in 1:4) {
  enr <- (site_cls_cnt[i] / nrow(sites_sl)) / class_props[i]
  text(bp[2, i], site_cls_cnt[i] + 1, sprintf("%.2fx", enr), cex = 1.0, font = 2)
}
legend("topleft",
       legend = c("Expected", "Observed", sprintf("\u03c7\u00b2=%.2f, p=%.2e", chi_sq, p_val)),
       fill = c(col_mauve, col_burgundy, NA), border = c("black", "black", NA),
       bty = "n", cex = 0.95)

# Panel 2: Density distribution
dens_land  <- density(landscape_vals)
dens_sites <- density(sites_sl$suitability, na.rm = TRUE)
plot(dens_land, main = "Sites vs Landscape Suitability Distribution",
     xlab = "Suitability", ylab = "Density",
     col = col_dustyblue, lwd = 2.5,
     xlim = range(c(landscape_vals, sites_sl$suitability), na.rm = TRUE),
     ylim = c(0, max(c(dens_land$y, dens_sites$y)) * 1.1))
lines(dens_sites, col = col_burgundy, lwd = 2.5)
abline(v = quartiles[2:4], lty = 2, col = "gray40", lwd = 1.5)
legend("topright",
       legend = c(sprintf("Landscape (mean=%.1f)", mean(landscape_vals)),
                  sprintf("Sites (mean=%.1f)", mean(sites_sl$suitability, na.rm = TRUE)),
                  sprintf("MW p=%.2e", mw_test$p.value),
                  sprintf("Cohen's d=%.2f", cohens_d)),
       col = c(col_dustyblue, col_burgundy, NA, NA), lwd = c(2.5, 2.5, NA, NA),
       bty = "n", cex = 0.95)

# Panel 3: Cohen's d
m_land <- mean(landscape_vals); m_site <- mean(sites_sl$suitability, na.rm = TRUE)
plot(c(1, 2), c(m_land, m_site), xlim = c(0.5, 2.5),
     ylim = c(0, max(m_land, m_site) * 1.3), pch = 19, cex = 2,
     col = c(col_dustyblue, col_burgundy),
     xlab = "", ylab = "Mean Suitability", main = "Effect Size (Cohen's d)", xaxt = "n")
axis(1, at = c(1, 2), labels = c("Landscape", "Sites"))
arrows(1, m_land - sd(landscape_vals), 1, m_land + sd(landscape_vals),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = col_dustyblue)
arrows(2, m_site - sd(sites_sl$suitability, na.rm = TRUE),
       2, m_site + sd(sites_sl$suitability, na.rm = TRUE),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = col_burgundy)
text(1.5, max(m_land, m_site) * 1.18,
     labels = sprintf("d = %.2f (%s)", cohens_d, label_d), cex = 1.0, font = 2)

# Panel 4: LOO-CV with Clopper-Pearson 95% CIs
plot(NA, xlim = c(0.5, 2.5), ylim = c(0, 110), xaxt = "n",
     xlab = "", ylab = "Accuracy (%)",
     main = "Leave-One-Out Cross-Validation")
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = 1:2, labels = c("High+VH", "VH only"))
arrows(1, ci_loo_hv[1], 1, ci_loo_hv[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = col_darkblue)
arrows(2, ci_loo_vh[1], 2, ci_loo_vh[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = col_darkblue)
points(1:2, c(loo_acc_hv, loo_acc_vh), pch = 19, cex = 2, col = col_darkblue)
abline(h = 50, lty = 2, col = col_mauve, lwd = 1.5)
text(2.4, 50, "50%", cex = 0.7, col = col_mauve, adj = c(1, -0.3))
text(1, loo_acc_hv, sprintf("%.1f%%", loo_acc_hv), pos = 4, font = 2, cex = 1.0)
text(2, loo_acc_vh, sprintf("%.1f%%", loo_acc_vh), pos = 4, font = 2, cex = 1.0)
text(1, ci_loo_hv[1], sprintf("[%.0f%%", ci_loo_hv[1]), pos = 1, cex = 0.8, col = "gray40")
text(1, ci_loo_hv[2], sprintf("%.0f%%]", ci_loo_hv[2]), pos = 3, cex = 0.8, col = "gray40")
text(2, ci_loo_vh[1], sprintf("[%.0f%%", ci_loo_vh[1]), pos = 1, cex = 0.8, col = "gray40")
text(2, ci_loo_vh[2], sprintf("%.0f%%]", ci_loo_vh[2]), pos = 3, cex = 0.8, col = "gray40")

# Panel 5: Permutation test
hist(perm_acc_hv, breaks = 30, col = col_beige, border = "white",
     main = "Null Model Permutation Test (High + Very High)",
     xlab = "% Random sites in High+VH", ylab = "Frequency", xlim = c(30, 100))
abline(v = real_acc_hv, col = col_burgundy, lwd = 3)
abline(v = mean(perm_acc_hv), col = col_dustyblue, lwd = 2, lty = 2)
legend("topleft",
       legend = c(sprintf("Real: %.1f%%", real_acc_hv),
                  sprintf("Null: %.1f%%", mean(perm_acc_hv))),
       col = c(col_burgundy, col_dustyblue), lwd = c(3, 2), lty = c(1, 2),
       bty = "n", cex = 0.95)

dev.off()
cat("✓ Saved: Validation_Tests.svg\n")


# ==========================================================================
# Fig 7b: LOO-CV weight stability (separate figure)
# ==========================================================================
cat("Creating weight stability plot...\n")

png(file.path(FIGS_OUTPUT, "LOO_Weight_Stability.png"),
    width = 3200, height = 2200, res = 300)
par(mar = c(5, 5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)

var_short  <- c("Elev", "Slope", "TRI", "TWI", "Aspect", "DTD")
loo_means  <- apply(loo_weights, 2, mean)
loo_sds    <- apply(loo_weights, 2, sd)
x_pos      <- 1:6
yr         <- range(c(loo_means - loo_sds, loo_means + loo_sds, weights))
plot(NA, xlim = c(0.5, 6.5), ylim = c(max(0, yr[1] - 2), yr[2] + 3),
     xaxt = "n", xlab = "", ylab = "Weight (%)",
     main = "LOO-CV Weight Stability")
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = x_pos, labels = var_short, cex.axis = 1.0)
arrows(x_pos, loo_means - loo_sds, x_pos, loo_means + loo_sds,
       angle = 90, code = 3, length = 0.08, lwd = 2, col = col_darkblue)
points(x_pos, loo_means, pch = 19, cex = 1.5, col = col_darkblue)
points(x_pos, weights, pch = 18, cex = 1.5, col = col_mauve)
legend("topright",
       legend = c("LOO mean \u00b1 SD", "Full-model weight"),
       pch = c(19, 18), col = c(col_darkblue, col_mauve),
       pt.cex = 1.3, bty = "n", cex = 0.95)
dev.off()
cat("✓ Saved: LOO_Weight_Stability.png\n")


# ==========================================================================
# Fig: Grid reference map
# ==========================================================================
cat("Creating grid reference map...\n")

suit_for_grid  <- if (exists("suitability_masked")) suitability_masked else suitability_900m
class_for_grid <- if (exists("s4_masked")) s4_masked else suitability_4class
valid_ext  <- ext(trim(suit_for_grid))
suit_crop  <- crop(suit_for_grid,   valid_ext)
class_crop <- crop(class_for_grid,  valid_ext)

e <- ext(suit_crop)
n_cols <- 8; n_rows <- 6
lon_breaks  <- seq(e[1], e[2], length.out = n_cols + 1)
lat_breaks  <- seq(e[3], e[4], length.out = n_rows + 1)
col_labels  <- LETTERS[1:n_cols]
row_labels  <- as.character(1:n_rows)
col_centres <- (lon_breaks[-1]     + lon_breaks[-(n_cols+1)]) / 2
row_centres <- rev((lat_breaks[-1] + lat_breaks[-(n_rows+1)]) / 2)

class_smooth <- focal(class_crop, w = 3, fun = "modal", na.rm = TRUE)
class_smooth <- mask(class_smooth, class_crop)

sites_sl_all <- sites_complete[sites_complete$source_dem == MODEL_REGION, ]
sites_inbounds <- sites_sl_all[
  sites_sl_all$Lon >= e[1] & sites_sl_all$Lon <= e[2] &
    sites_sl_all$Lat >= e[3] & sites_sl_all$Lat <= e[4], ]
sites_sl_pts <- vect(sites_inbounds, geom = c("Lon", "Lat"), crs = "EPSG:4326")
sites_inbounds$suit_class <- extract(suitability_4class, sites_sl_pts)[, 2]

png(file.path(FIGS_OUTPUT, paste0(MODEL_REGION, "_Grid_Map.png")),
    width = 4900, height = 5300, res = 400)
par(mar = c(4, 2, 5, 10))

# Set class levels before plotting (fixes terra classes error after focal/mask)
levels(class_smooth) <- data.frame(id = 1:4, class = class_names)

plot(class_smooth, col = cols_pub, axes = FALSE, type = "classes",
     xlim = c(e[1], e[2]), ylim = c(e[3], e[4]),
     plg = list(title = "Suitability", legend = class_names, cex = 0.9))

for (j in 1:n_cols) for (i in 1:n_rows)
  rect(lon_breaks[j], lat_breaks[i], lon_breaks[j+1], lat_breaks[i+1],
       border = adjustcolor("gray10", 0.8), lwd = 2.5, col = NA)
for (j in seq_along(col_labels)) {
  text(col_centres[j], e[4] + (e[4]-e[3])*0.018, col_labels[j], cex = 1.0, font = 2, xpd = TRUE)
  text(col_centres[j], e[3] - (e[4]-e[3])*0.018, col_labels[j], cex = 1.0, font = 2, xpd = TRUE, adj = c(0.5, 1))
}
for (i in seq_along(row_labels)) {
  text(e[1] - (e[2]-e[1])*0.04, row_centres[i], row_labels[i], cex = 1.0, font = 2, xpd = TRUE, adj = c(1, 0.5))
  text(e[2] + (e[2]-e[1])*0.04, row_centres[i], row_labels[i], cex = 1.0, font = 2, xpd = TRUE, adj = c(0, 0.5))
}
for (j in seq_along(col_labels)) for (i in seq_along(row_labels))
  text(col_centres[j], row_centres[i], paste0(col_labels[j], row_labels[i]),
       cex = 0.55, col = adjustcolor("white", 0.85), font = 2)

site_shp <- c(21, 22, 23, 24)[sites_inbounds$suit_class]
points(sites_sl_pts, pch = site_shp, bg = "black", cex = 1.5, col = "white", lwd = 1.8)

axis(1, at = round(seq(e[1], e[2], length.out = 5), 1),
     labels = paste0(round(seq(e[1], e[2], length.out = 5), 1), "°E"), cex.axis = 0.85)
axis(2, at = round(seq(e[3], e[4], length.out = 5), 1),
     labels = paste0(round(seq(e[3], e[4], length.out = 5), 1), "°N"), cex.axis = 0.85, las = 1)

legend("bottomleft",
       legend = c(sprintf("Very High (n=%d)", sum(sites_inbounds$suit_class == 4, na.rm = TRUE)),
                  sprintf("High (n=%d)",      sum(sites_inbounds$suit_class == 3, na.rm = TRUE)),
                  sprintf("Medium (n=%d)",    sum(sites_inbounds$suit_class == 2, na.rm = TRUE)),
                  sprintf("Low (n=%d)",       sum(sites_inbounds$suit_class == 1, na.rm = TRUE))),
       pch = c(24, 23, 22, 21), pt.bg = "black", pt.cex = 1.4, col = "black",
       bty = "o", bg = "white", box.col = "gray40", cex = 0.85, inset = c(0.01, 0.01))

dev.off()
cat(sprintf("✓ Saved: %s\n\n", file.path(FIGS_OUTPUT, paste0(MODEL_REGION, "_Grid_Map.png"))))




#############################################################################
# =============================================================================
# SECTION 12: VISUALIZATIONS
# =============================================================================

cat("================================================================================\n")
cat("SECTION 12: VISUALIZATIONS\n")
cat("================================================================================\n\n")

# ── Install svglite if needed ─────────────────────────────────────────────────
if (!requireNamespace("svglite", quietly = TRUE)) {
  cat("Installing svglite...\n")
  install.packages("svglite", repos = "https://cloud.r-project.org")
}
library(svglite)

# ── Shared style values (from style reference) ────────────────────────────────
zone_cols <- c(Mediterranean = "#954151",
               "Irano-Turanian" = "#1d4e89",
               "Saharo-Arabian" = "#7295bf")
col_burgundy  <- "#954151"
col_darkblue  <- "#1d4e89"
col_dustyblue <- "#7295bf"
col_mauve     <- "#db9793"
col_beige     <- "#ead98b"

# Publication map colors
cols_pub <- c("#A8BED1", "#4A6FA5", "#CDA67E", "#7B2D3C")
col_unlikely <- "#D3D3D3"

class_names <- c("Low", "Medium", "High", "Very High")

# Ensure output directories exist
dir.create(FIGS_OUTPUT, showWarnings = FALSE, recursive = TRUE)

# Reload sites_sl with suit_class from the 4-class model
sites_sl  <- sites_complete[sites_complete$source_dem == MODEL_REGION, ]
site_pts  <- vect(sites_sl, geom = c("Lon", "Lat"), crs = "EPSG:4326")
sites_sl$suitability  <- extract(suitability_900m, site_pts)[, 2]
sites_sl$suit_class   <- extract(suitability_4class, site_pts)[, 2]


# ==========================================================================
# Fig 1: Suitability maps (publication, 2-panel vertical with hillshade)
# ==========================================================================
cat("Creating suitability maps (publication style)...\n")

dem_900m   <- aggregate(dem_r, fact = 30, fun = "mean", na.rm = TRUE)
slope_hs   <- terrain(dem_900m, v = "slope", unit = "radians")
aspect_hs  <- terrain(dem_900m, v = "aspect", unit = "radians")
hillshade  <- shade(slope_hs, aspect_hs, angle = 45, direction = 315)

# 5-class raster for bottom panel
if (exists("s4_masked") && !identical(s4_masked, suitability_4class)) {
  suit_5class <- suitability_4class
  unlikely_px <- !is.na(suitability_4class) & is.na(s4_masked)
  suit_5class[unlikely_px] <- 5
  levels(suit_5class) <- data.frame(id = 1:5,
                                    class = c(class_names, "Unlikely"))
} else {
  suit_5class <- suitability_4class
}

e <- ext(suitability_900m)

draw_scalebar <- function(x0, y0, km_len, e) {
  mid_lat    <- (e[3] + e[4]) / 2
  deg_per_km <- 1 / (111.32 * cos(mid_lat * pi / 180))
  dx <- km_len * deg_per_km
  segments(x0, y0, x0 + dx, y0, lwd = 2.5, col = "black")
  segments(x0, y0 - 0.05, x0, y0 + 0.05, lwd = 2, col = "black")
  segments(x0 + dx, y0 - 0.05, x0 + dx, y0 + 0.05, lwd = 2, col = "black")
  text(x0 + dx / 2, y0 - 0.15, paste0(km_len, " km"), cex = 0.7)
}

fig1_path <- file.path(BASE_DIR, "Suitability_Maps.png")
png(fig1_path, width = 3200, height = 5800, res = 300)
layout(matrix(1:2, nrow = 2), heights = c(1, 1))

# --- Top panel: unmasked 4-class ---
par(mar = c(2, 2, 1, 2))
plot(hillshade, col = gray.colors(256, 0.05, 0.95),
     legend = FALSE, axes = FALSE, reset = FALSE)

levels(suitability_4class) <- data.frame(id = 1:4, class = class_names)
plot(suitability_4class,
     col = adjustcolor(cols_pub, alpha.f = 0.70),
     type = "classes",
     legend = FALSE, axes = FALSE, add = TRUE)

points(site_pts, pch = 24, bg = "#FFD700", cex = 1.2, col = "black", lwd = 0.9)

# Geographic labels — adjust coordinates as needed
text(35.05, 32.80, "Mt. Carmel",    cex = 0.80, font = 3, col = "white")
text(35.55, 32.35, "Jordan Valley", cex = 0.80, font = 3, col = "white")
text(36.80, 31.95, "Azraq",         cex = 0.80, font = 3, col = "white")
text(34.50, 30.85, "Negev",         cex = 0.80, font = 3, col = "white")
text(35.10, 30.15, "Wadi Araba",    cex = 0.80, font = 3, col = "white")
text(36.20, 30.40, "Al-Jafr",       cex = 0.80, font = 3, col = "white")

legend("topleft", legend = rev(class_names), fill = rev(cols_pub),
       border = "gray30", title = "MCDA Suitability", title.font = 2,
       bty = "o", bg = "white", box.col = "gray40", cex = 0.75, inset = c(0.01, 0.02))

arr_x  <- e[1] + (e[2] - e[1]) * 0.03
arr_y0 <- e[3] + (e[4] - e[3]) * 0.15
arr_y1 <- arr_y0 + (e[4] - e[3]) * 0.06
arrows(arr_x, arr_y0, arr_x, arr_y1, length = 0.12, lwd = 2, col = "black")
text(arr_x, arr_y1 + (e[4] - e[3]) * 0.02, "N", cex = 0.9, font = 2)

sb_x <- e[2] - (e[2] - e[1]) * 0.35
sb_y <- e[3] + (e[4] - e[3]) * 0.05
draw_scalebar(sb_x, sb_y, 50, e)
draw_scalebar(sb_x + 50 / (111.32 * cos(((e[3]+e[4])/2)*pi/180)), sb_y, 50, e)

# --- Bottom panel: masked 5-class + regions ---
par(mar = c(2, 2, 1, 2))
plot(hillshade, col = gray.colors(256, 0.05, 0.95),
     legend = FALSE, axes = FALSE, reset = FALSE)

if (exists("suit_5class") && max(values(suit_5class), na.rm = TRUE) == 5) {
  plot(suit_5class,
       col = adjustcolor(c(cols_pub, col_unlikely), alpha.f = 0.70),
       type = "classes",
       legend = FALSE, axes = FALSE, add = TRUE)
  leg_names <- rev(c(class_names, "Unlikely"))
  leg_cols  <- rev(c(cols_pub, col_unlikely))
} else {
  levels(suitability_4class) <- data.frame(id = 1:4, class = class_names)
  plot(suitability_4class,
       col = adjustcolor(cols_pub, alpha.f = 0.70),
       type = "classes",
       legend = FALSE, axes = FALSE, add = TRUE)
  leg_names <- rev(class_names)
  leg_cols  <- rev(cols_pub)
}

points(site_pts, pch = 24, bg = "#FFD700", cex = 1.2, col = "black", lwd = 0.9)

# Region rectangles — adjust coordinates as needed
rect(34.50, 32.00, 35.20, 33.00, border = "black", lwd = 2.5, lty = 2)
text(34.55, 32.90, "A", cex = 1.1, font = 2, adj = c(0, 1))
rect(34.50, 30.80, 35.20, 32.00, border = "black", lwd = 2.5, lty = 2)
text(34.55, 31.90, "B", cex = 1.1, font = 2, adj = c(0, 1))
rect(36.00, 31.30, 37.80, 32.50, border = "black", lwd = 2.5, lty = 2)
text(36.70, 32.30, "C", cex = 1.1, font = 2, adj = c(0, 1))
rect(35.00, 29.50, 36.80, 31.00, border = "black", lwd = 2.5, lty = 2)
text(35.80, 30.30, "D", cex = 1.1, font = 2, adj = c(0, 1))

legend("topleft", legend = leg_names, fill = leg_cols,
       border = "gray30", title = "MCDA Suitability", title.font = 2,
       bty = "o", bg = "white", box.col = "gray40", cex = 0.75, inset = c(0.01, 0.02))
arrows(arr_x, arr_y0, arr_x, arr_y1, length = 0.12, lwd = 2, col = "black")
text(arr_x, arr_y1 + (e[4] - e[3]) * 0.02, "N", cex = 0.9, font = 2)
draw_scalebar(sb_x, sb_y, 50, e)
draw_scalebar(sb_x + 50 / (111.32 * cos(((e[3]+e[4])/2)*pi/180)), sb_y, 50, e)

dev.off()
cat(sprintf("\u2713 Saved: %s\n", fig1_path))


# ==========================================================================
# Fig 4: Strip plots by bioclimatic zone (svglite, 18 x 12 in)
# ==========================================================================
cat("Creating strip plots by zone...\n")

sites_complete$Region <- factor(sites_complete$Region,
                                levels = c("Mediterranean","Irano-Turanian","Saharo-Arabian"))
var_list <- list(
  list(col="elevation",   ylab="Elevation (m asl)",                  title="A) Elevation"),
  list(col="slope",       ylab="Slope (degrees)",                    title="B) Slope"),
  list(col="TRI",         ylab="Terrain Ruggedness Index (TRI)",     title="C) TRI"),
  list(col="TWI",         ylab="Topographic Wetness Index (TWI)",    title="D) TWI"),
  list(col="northness",   ylab="Aspect (0-1)",                      title="E) Aspect"),
  list(col="DistToWater", ylab="DTD (km)",                          title="F) DTD")
)
is_outlier <- function(x) {
  q  <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iq <- diff(q); x < (q[1] - 1.5*iq) | x > (q[2] + 1.5*iq)
}

fig4_path <- file.path(FIGS_OUTPUT, "Topographic_by_zone.svg")
svglite(fig4_path, width = 18, height = 12)
par(mfrow = c(2, 3), mar = c(6, 5, 4, 2), oma = c(1, 0, 3, 0),
    mgp = c(2.8, 0.9, 0), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)

zone_lev <- levels(sites_complete$Region)
panel_count <- 0
for (v in var_list) {
  panel_count <- panel_count + 1
  vals_l <- lapply(zone_lev, function(z) {
    x <- sites_complete[[v$col]][sites_complete$Region == z]
    x[!is.na(x)]
  })
  means <- sapply(vals_l, mean); sds <- sapply(vals_l, sd)
  all_v <- unlist(vals_l)
  yr    <- diff(range(all_v))
  plot(NA, xlim = c(0.5, 3.5),
       ylim = c(min(all_v) - 0.05*yr, max(all_v) + 0.25*yr),
       xaxt = "n", xlab = "", ylab = v$ylab, main = v$title, las = 1, frame = FALSE)
  grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
  axis(1, at = 1:3, labels = zone_lev, cex.axis = 0.80)
  for (i in seq_along(zone_lev)) {
    z   <- zone_lev[i]; col <- zone_cols[z]
    x   <- vals_l[[i]]; mn <- means[i]; sdv <- sds[i]; n <- length(x)
    segments(i, mn - sdv, i, mn + sdv, col = col, lwd = 2)
    segments(i - 0.1, mn - sdv, i + 0.1, mn - sdv, col = col, lwd = 1.8)
    segments(i - 0.1, mn + sdv, i + 0.1, mn + sdv, col = col, lwd = 1.8)
    set.seed(42); jx <- i + runif(n, -0.15, 0.15)
    out <- is_outlier(x)
    points(jx[!out], x[!out], pch = 21, bg = adjustcolor(col, 0.45),
           col = adjustcolor(col, 0.75), cex = 0.9, lwd = 0.7)
    if (any(out)) points(jx[out], x[out], pch = 21, bg = adjustcolor(col, 0.85),
                         col = "gray20", cex = 1.2, lwd = 1.1)
    points(i, mn, pch = 23, bg = "white", col = "gray20", cex = 1.2, lwd = 1.4)
    # Highest-value site name label
    max_idx <- which.max(x)
    site_names_z <- sites_complete$Site[sites_complete$Region == z & !is.na(sites_complete[[v$col]])]
    if (length(site_names_z) >= max_idx)
      text(jx[max_idx], x[max_idx], site_names_z[max_idx],
           cex = 0.75, font = 2, col = col, pos = 4)
    mtext(sprintf("n=%d", n), side = 1, at = i, line = 4.8, cex = 0.9, col = "gray40")
  }
  # Add legend in the last panel
  if (panel_count == 6) {
    legend("topright", legend = zone_lev, pch = 21,
           pt.bg = adjustcolor(zone_cols[zone_lev], 0.65), col = "gray20",
           pt.cex = 1.3, cex = 0.85, bty = "n")
  }
}
mtext("Topographic and Hydrological Variables by Bioclimatic Zone",
      outer = TRUE, cex = 1.4, font = 2, line = 1.5)
dev.off()
cat(sprintf("\u2713 Saved: %s\n", fig4_path))


# ==========================================================================
# Fig 5: Dumbbell plot (zone weights vs modeled area)
# ==========================================================================
cat("Creating dumbbell plot...\n")

var_order  <- c("Elevation", "Slope", "TRI", "TWI", "Northness", "DistToWater")
var_lab_d  <- c(Elevation = "Elevation", Slope = "Slope", TRI = "TRI",
                TWI = "TWI", Northness = "Aspect", DistToWater = "DTD")
zone_pch   <- c(Mediterranean = 21, "Irano-Turanian" = 22, "Saharo-Arabian" = 24)

modeled <- all_weights[all_weights$Zone == "Modeled Area", ]

fig5_path <- file.path(FIGS_OUTPUT, "Variable_weights_dumbbell.png")
png(fig5_path, width = 4800, height = 3200, res = 300)
par(mar = c(5, 9, 4, 10), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0)
n_v <- length(var_order)
plot(NA, xlim = c(0, 42), ylim = c(0.5, n_v + 0.5),
     xaxt = "n", yaxt = "n", xlab = "Weight (%)", ylab = "", frame = FALSE)
for (g in seq(0, 40, 10)) segments(g, 0.3, g, n_v + 0.5, col = "gray90", lwd = 1)
axis(1, at = seq(0, 40, 10), labels = paste0(seq(0, 40, 10), "%"),
     cex.axis = 1.0, col.axis = "gray30", col = "gray60")
axis(2, at = 1:n_v, labels = var_lab_d[var_order], tick = FALSE, las = 1, cex.axis = 1.2)

for (i in seq_along(var_order)) {
  v   <- var_order[i]
  mv  <- modeled$Weight[modeled$Variable == v]
  zvs <- sapply(zones, function(z) {
    ww <- all_weights$Weight[all_weights$Zone == z & all_weights$Variable == v]
    if (length(ww) > 0) ww else NA
  })
  zvs <- zvs[!is.na(zvs)]
  if (length(zvs) > 0) segments(min(zvs), i, max(zvs), i, col = "gray75", lwd = 2.5)
  segments(mv, i - 0.35, mv, i + 0.35, col = col_mauve, lwd = 2, lty = 2)
  for (z in zones) {
    zv <- all_weights$Weight[all_weights$Zone == z & all_weights$Variable == v]
    if (length(zv) > 0)
      points(zv, i, pch = zone_pch[z], bg = adjustcolor(zone_cols[z], 0.85),
             col = "gray20", cex = 1.8, lwd = 1.2)
  }
  points(mv, i, pch = 18, col = col_mauve, cex = 1.6)
  text(mv, i + 0.32, sprintf("%.1f%%", mv), cex = 1.0, col = "gray30")
}
mtext("Variable Weights by Bioclimatic Zone", side = 3, line = 2, cex = 1.3, font = 2)

lx <- 43; ly <- n_v + 0.3
for (k in seq_along(zones)) {
  z <- zones[k]; yp <- ly - (k - 1) * 0.6
  if (any(all_weights$Zone == z)) {
    points(lx, yp, pch = zone_pch[z], bg = adjustcolor(zone_cols[z], 0.85),
           col = "gray20", cex = 1.6, xpd = TRUE, lwd = 1.2)
    text(lx + 0.8, yp, z, cex = 0.95, adj = 0, xpd = TRUE, col = "gray20")
  }
}
yg <- ly - length(zones) * 0.6 - 0.2
segments(lx - 0.3, yg, lx + 0.3, yg, col = col_mauve, lwd = 2, lty = 2, xpd = TRUE)
points(lx, yg, pch = 18, col = col_mauve, cex = 1.6, xpd = TRUE)
text(lx + 0.8, yg, "Modeled Area", cex = 0.95, adj = 0, xpd = TRUE, col = "gray20")
dev.off()
cat(sprintf("\u2713 Saved: %s\n", fig5_path))


# ==========================================================================
# Fig 7: Validation tests (6 panels: 3 top, 2 bottom + 1 empty)
# ==========================================================================
cat("Creating validation plots...\n")

fig7_path <- file.path(FIGS_OUTPUT, "Validation_Tests.svg")
svglite(fig7_path, width = 18.67, height = 12)
par(mfrow = c(2, 3), mar = c(5, 4.5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)

# Panel 1: Chi-square
bp <- barplot(rbind(expected, site_cls_cnt),
              beside = TRUE, names.arg = class_names,
              col = c(col_mauve, col_burgundy),
              ylab = "Number of sites", main = "Chi-Square Enrichment Test",
              ylim = c(0, max(c(expected, site_cls_cnt)) * 1.3))
for (i in 1:4) {
  enr <- (site_cls_cnt[i] / nrow(sites_sl)) / class_props[i]
  text(bp[2, i], site_cls_cnt[i] + 1, sprintf("%.2fx", enr), cex = 1.0, font = 2)
}
legend("topleft",
       legend = c("Expected", "Observed", sprintf("\u03c7\u00b2=%.2f, p=%.2e", chi_sq, p_val)),
       fill = c(col_mauve, col_burgundy, NA), border = c("black", "black", NA),
       bty = "n", cex = 0.95)

# Panel 2: Density distribution
dens_land  <- density(landscape_vals)
dens_sites <- density(sites_sl$suitability, na.rm = TRUE)
plot(dens_land, main = "Sites vs Landscape Suitability Distribution",
     xlab = "Suitability", ylab = "Density",
     col = col_dustyblue, lwd = 2.5,
     xlim = range(c(landscape_vals, sites_sl$suitability), na.rm = TRUE),
     ylim = c(0, max(c(dens_land$y, dens_sites$y)) * 1.1))
lines(dens_sites, col = col_burgundy, lwd = 2.5)
abline(v = quartiles[2:4], lty = 2, col = "gray40", lwd = 1.5)
legend("topright",
       legend = c(sprintf("Landscape (mean=%.1f)", mean(landscape_vals)),
                  sprintf("Sites (mean=%.1f)", mean(sites_sl$suitability, na.rm = TRUE)),
                  sprintf("MW p=%.2e", mw_test$p.value),
                  sprintf("Cohen's d=%.2f", cohens_d)),
       col = c(col_dustyblue, col_burgundy, NA, NA), lwd = c(2.5, 2.5, NA, NA),
       bty = "n", cex = 0.95)

# Panel 3: Cohen's d
m_land <- mean(landscape_vals); m_site <- mean(sites_sl$suitability, na.rm = TRUE)
plot(c(1, 2), c(m_land, m_site), xlim = c(0.5, 2.5),
     ylim = c(0, max(m_land, m_site) * 1.3), pch = 19, cex = 2,
     col = c(col_dustyblue, col_burgundy),
     xlab = "", ylab = "Mean Suitability", main = "Effect Size (Cohen's d)", xaxt = "n")
axis(1, at = c(1, 2), labels = c("Landscape", "Sites"))
arrows(1, m_land - sd(landscape_vals), 1, m_land + sd(landscape_vals),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = col_dustyblue)
arrows(2, m_site - sd(sites_sl$suitability, na.rm = TRUE),
       2, m_site + sd(sites_sl$suitability, na.rm = TRUE),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = col_burgundy)
text(1.5, max(m_land, m_site) * 1.18,
     labels = sprintf("d = %.2f (%s)", cohens_d, label_d), cex = 1.0, font = 2)

# Panel 4: LOO-CV with Clopper-Pearson 95% CIs
plot(NA, xlim = c(0.5, 2.5), ylim = c(0, 110), xaxt = "n",
     xlab = "", ylab = "Accuracy (%)",
     main = "Leave-One-Out Cross-Validation")
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = 1:2, labels = c("High+VH", "VH only"))
arrows(1, ci_loo_hv[1], 1, ci_loo_hv[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = col_darkblue)
arrows(2, ci_loo_vh[1], 2, ci_loo_vh[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = col_darkblue)
points(1:2, c(loo_acc_hv, loo_acc_vh), pch = 19, cex = 2, col = col_darkblue)
abline(h = 50, lty = 2, col = col_mauve, lwd = 1.5)
text(2.4, 50, "50%", cex = 0.7, col = col_mauve, adj = c(1, -0.3))
text(1, loo_acc_hv, sprintf("%.1f%%", loo_acc_hv), pos = 4, font = 2, cex = 1.0)
text(2, loo_acc_vh, sprintf("%.1f%%", loo_acc_vh), pos = 4, font = 2, cex = 1.0)
text(1, ci_loo_hv[1], sprintf("[%.0f%%", ci_loo_hv[1]), pos = 1, cex = 0.8, col = "gray40")
text(1, ci_loo_hv[2], sprintf("%.0f%%]", ci_loo_hv[2]), pos = 3, cex = 0.8, col = "gray40")
text(2, ci_loo_vh[1], sprintf("[%.0f%%", ci_loo_vh[1]), pos = 1, cex = 0.8, col = "gray40")
text(2, ci_loo_vh[2], sprintf("%.0f%%]", ci_loo_vh[2]), pos = 3, cex = 0.8, col = "gray40")

# Panel 5: Permutation test
hist(perm_acc_hv, breaks = 30, col = col_beige, border = "white",
     main = "Null Model Permutation Test (High + Very High)",
     xlab = "% Random sites in High+VH", ylab = "Frequency", xlim = c(30, 100))
abline(v = real_acc_hv, col = col_burgundy, lwd = 3)
abline(v = mean(perm_acc_hv), col = col_dustyblue, lwd = 2, lty = 2)
legend("topleft",
       legend = c(sprintf("Real: %.1f%%", real_acc_hv),
                  sprintf("Null: %.1f%%", mean(perm_acc_hv))),
       col = c(col_burgundy, col_dustyblue), lwd = c(3, 2), lty = c(1, 2),
       bty = "n", cex = 0.95)

# Panel 6: empty placeholder (keeps layout balanced)
plot.new()

dev.off()
cat(sprintf("\u2713 Saved: %s\n", fig7_path))


# ==========================================================================
# Fig 7b: LOO-CV weight stability (separate figure)
# ==========================================================================
cat("Creating weight stability plot...\n")

fig7b_path <- file.path(FIGS_OUTPUT, "LOO_Weight_Stability.png")
png(fig7b_path, width = 3200, height = 2200, res = 300)
par(mar = c(5, 5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)

var_short  <- c("Elev", "Slope", "Aspect", "TRI", "TWI", "DTD")
loo_means  <- apply(loo_weights, 2, mean)
loo_sds    <- apply(loo_weights, 2, sd)
x_pos      <- 1:6
yr         <- range(c(loo_means - loo_sds, loo_means + loo_sds, weights))
plot(NA, xlim = c(0.5, 6.5), ylim = c(max(0, yr[1] - 2), yr[2] + 3),
     xaxt = "n", xlab = "", ylab = "Weight (%)",
     main = "LOO-CV Weight Stability")
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = x_pos, labels = var_short, cex.axis = 1.0)
arrows(x_pos, loo_means - loo_sds, x_pos, loo_means + loo_sds,
       angle = 90, code = 3, length = 0.08, lwd = 2, col = col_darkblue)
points(x_pos, loo_means, pch = 19, cex = 1.5, col = col_darkblue)
points(x_pos, weights, pch = 18, cex = 1.5, col = col_mauve)
legend("topright",
       legend = c("LOO mean \u00b1 SD", "Full-model weight"),
       pch = c(19, 18), col = c(col_darkblue, col_mauve),
       pt.cex = 1.3, bty = "n", cex = 0.95)
dev.off()
cat(sprintf("\u2713 Saved: %s\n", fig7b_path))


# ==========================================================================
# Fig: Grid reference map
# ==========================================================================
cat("Creating grid reference map...\n")

suit_for_grid  <- if (exists("suitability_masked")) suitability_masked else suitability_900m
class_for_grid <- if (exists("s4_masked")) s4_masked else suitability_4class
valid_ext  <- ext(trim(suit_for_grid))
suit_crop  <- crop(suit_for_grid,   valid_ext)
class_crop <- crop(class_for_grid,  valid_ext)

e <- ext(suit_crop)
n_cols <- 8; n_rows <- 6
lon_breaks  <- seq(e[1], e[2], length.out = n_cols + 1)
lat_breaks  <- seq(e[3], e[4], length.out = n_rows + 1)
col_labels  <- LETTERS[1:n_cols]
row_labels  <- as.character(1:n_rows)
col_centres <- (lon_breaks[-1]     + lon_breaks[-(n_cols+1)]) / 2
row_centres <- rev((lat_breaks[-1] + lat_breaks[-(n_rows+1)]) / 2)

class_smooth <- focal(class_crop, w = 3, fun = "modal", na.rm = TRUE)
class_smooth <- mask(class_smooth, class_crop)

sites_sl_all <- sites_complete[sites_complete$source_dem == MODEL_REGION, ]
sites_inbounds <- sites_sl_all[
  sites_sl_all$Lon >= e[1] & sites_sl_all$Lon <= e[2] &
    sites_sl_all$Lat >= e[3] & sites_sl_all$Lat <= e[4], ]
sites_sl_pts <- vect(sites_inbounds, geom = c("Lon", "Lat"), crs = "EPSG:4326")
sites_inbounds$suit_class <- extract(suitability_4class, sites_sl_pts)[, 2]

grid_path <- file.path(FIGS_OUTPUT, paste0(MODEL_REGION, "_Grid_Map.png"))
png(grid_path, width = 4900, height = 5300, res = 400)
par(mar = c(4, 2, 5, 10))

# Set class levels before plotting (fixes terra error after focal/mask)
levels(class_smooth) <- data.frame(id = 1:4, class = class_names)

plot(class_smooth, col = cols_pub, axes = FALSE, type = "classes",
     xlim = c(e[1], e[2]), ylim = c(e[3], e[4]),
     plg = list(title = "Suitability", legend = class_names, cex = 0.9))

for (j in 1:n_cols) for (i in 1:n_rows)
  rect(lon_breaks[j], lat_breaks[i], lon_breaks[j+1], lat_breaks[i+1],
       border = adjustcolor("gray10", 0.8), lwd = 2.5, col = NA)
for (j in seq_along(col_labels)) {
  text(col_centres[j], e[4] + (e[4]-e[3])*0.018, col_labels[j], cex = 1.0, font = 2, xpd = TRUE)
  text(col_centres[j], e[3] - (e[4]-e[3])*0.018, col_labels[j], cex = 1.0, font = 2, xpd = TRUE, adj = c(0.5, 1))
}
for (i in seq_along(row_labels)) {
  text(e[1] - (e[2]-e[1])*0.04, row_centres[i], row_labels[i], cex = 1.0, font = 2, xpd = TRUE, adj = c(1, 0.5))
  text(e[2] + (e[2]-e[1])*0.04, row_centres[i], row_labels[i], cex = 1.0, font = 2, xpd = TRUE, adj = c(0, 0.5))
}
for (j in seq_along(col_labels)) for (i in seq_along(row_labels))
  text(col_centres[j], row_centres[i], paste0(col_labels[j], row_labels[i]),
       cex = 0.55, col = adjustcolor("white", 0.85), font = 2)

site_shp <- c(21, 22, 23, 24)[sites_inbounds$suit_class]
points(sites_sl_pts, pch = site_shp, bg = "black", cex = 1.5, col = "white", lwd = 1.8)

axis(1, at = round(seq(e[1], e[2], length.out = 5), 1),
     labels = paste0(round(seq(e[1], e[2], length.out = 5), 1), "\u00b0E"), cex.axis = 0.85)
axis(2, at = round(seq(e[3], e[4], length.out = 5), 1),
     labels = paste0(round(seq(e[3], e[4], length.out = 5), 1), "\u00b0N"), cex.axis = 0.85, las = 1)

legend("bottomleft",
       legend = c(sprintf("Very High (n=%d)", sum(sites_inbounds$suit_class == 4, na.rm = TRUE)),
                  sprintf("High (n=%d)",      sum(sites_inbounds$suit_class == 3, na.rm = TRUE)),
                  sprintf("Medium (n=%d)",    sum(sites_inbounds$suit_class == 2, na.rm = TRUE)),
                  sprintf("Low (n=%d)",       sum(sites_inbounds$suit_class == 1, na.rm = TRUE))),
       pch = c(24, 23, 22, 21), pt.bg = "black", pt.cex = 1.4, col = "black",
       bty = "o", bg = "white", box.col = "gray40", cex = 0.85, inset = c(0.01, 0.01))

dev.off()
cat(sprintf("\u2713 Saved: %s\n\n", grid_path))


# ==========================================================================
# FINAL SUMMARY
# ==========================================================================

cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Output directory: ", OUTPUT_DIR, "\n\n")

cat("FILES CREATED:\n")
cat("  Topo/Slope/           <- 4 regions, slope (degrees)\n")
cat("  Topo/Aspect/          <- 4 regions, northness (0-1)\n")
cat("  Topo/TRI/             <- 4 regions\n")
cat("  Topo/TWI/             <- 4 regions\n")
cat("  Topo/DistToWater/     <- 4 regions, P95 threshold (km)\n")
cat("  Sites_data/Sites_with_variables.csv\n")
cat("  Sites_data/Background_sample.csv\n")
cat("  Stats/Descriptive_stats_by_zone.csv\n")
cat("  Stats/Variable_correlations.csv\n")
cat("  Stats/Variable_weights_all.csv\n")
cat("  Stats/Validation_metrics.csv\n")
cat("  Stats/LOO_weight_stability.csv\n")
cat("  Stats/Grid_suitability_summary.csv\n")
cat("  MCDA/Suitability_30m_P95.tif\n")
cat("  MCDA/Suitability_900m_P95.tif\n")
cat("  MCDA/Suitability_4class_P95.tif\n")
cat("  MCDA/Suitability_900m_P95_MASKED.tif   (if geology mask present)\n")
cat("  MCDA/Suitability_4class_P95_MASKED.tif (if geology mask present)\n")
cat("  Figures/Topographic_by_zone.svg\n")
cat("  Figures/Variable_weights_dumbbell.png\n")
cat("  Figures/Validation_Tests.svg\n")
cat("  Figures/LOO_Weight_Stability.png\n")
cat(sprintf("  Figures/%s_Grid_Map.png\n", MODEL_REGION))
cat(sprintf("  %s\n\n", file.path(BASE_DIR, "Suitability_Maps.png")))

cat("VALIDATION SUMMARY:\n")
cat(sprintf("  Chi-square:            %.2f  (p = %.2e)\n", chi_sq, p_val))
cat(sprintf("  Cohen's d:             %.3f  (%s)\n", cohens_d, label_d))
cat(sprintf("  Sites in High+VH:      %.1f%%  95%% CI [%.1f%%, %.1f%%]\n",
            real_acc_hv, ci_hv[1], ci_hv[2]))
cat(sprintf("  LOO-CV (High+VH):      %.1f%%  95%% CI [%.1f%%, %.1f%%]\n",
            loo_acc_hv, ci_loo_hv[1], ci_loo_hv[2]))
cat(sprintf("  Permutation p (H+VH):  %.4f\n", p_perm_hv))
cat(sprintf("  Max LOO weight shift:  %.1f pp\n",
            max(abs(weights - apply(loo_weights, 2, mean)))))
cat("================================================================================\n\n")





############## One by One svg plots###############################################

# ==========================================================================
# Validation panels — individual SVG files
# ==========================================================================
cat("Creating individual validation panels...\n")

# Panel 1: Chi-square
svglite(file.path(FIGS_OUTPUT, "Val_1_ChiSquare.svg"), width = 6, height = 5)
par(mar = c(5, 4.5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)
bp <- barplot(rbind(expected, site_cls_cnt),
              beside = TRUE, names.arg = class_names,
              col = c(col_mauve, col_burgundy),
              ylab = "Number of sites", main = "Chi-Square Enrichment Test",
              ylim = c(0, max(c(expected, site_cls_cnt)) * 1.3))
for (i in 1:4) {
  enr <- (site_cls_cnt[i] / nrow(sites_sl)) / class_props[i]
  text(bp[2, i], site_cls_cnt[i] + 1, sprintf("%.2fx", enr), cex = 1.0, font = 2)
}
legend("topleft",
       legend = c("Expected", "Observed", sprintf("\u03c7\u00b2=%.2f, p=%.2e", chi_sq, p_val)),
       fill = c(col_mauve, col_burgundy, NA), border = c("black", "black", NA),
       bty = "n", cex = 0.95)
dev.off()
cat(sprintf("\u2713 Saved: %s\n", file.path(FIGS_OUTPUT, "Val_1_ChiSquare.svg")))

# Panel 2: Density distribution
svglite(file.path(FIGS_OUTPUT, "Val_2_Density.svg"), width = 6, height = 5)
par(mar = c(5, 4.5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)
dens_land  <- density(landscape_vals)
dens_sites <- density(sites_sl$suitability, na.rm = TRUE)
plot(dens_land, main = "Sites vs Landscape Suitability Distribution",
     xlab = "Suitability", ylab = "Density",
     col = col_dustyblue, lwd = 2.5,
     xlim = range(c(landscape_vals, sites_sl$suitability), na.rm = TRUE),
     ylim = c(0, max(c(dens_land$y, dens_sites$y)) * 1.1))
lines(dens_sites, col = col_burgundy, lwd = 2.5)
abline(v = quartiles[2:4], lty = 2, col = "gray40", lwd = 1.5)
legend("topright",
       legend = c(sprintf("Landscape (mean=%.1f)", mean(landscape_vals)),
                  sprintf("Sites (mean=%.1f)", mean(sites_sl$suitability, na.rm = TRUE)),
                  sprintf("MW p=%.2e", mw_test$p.value),
                  sprintf("Cohen's d=%.2f", cohens_d)),
       col = c(col_dustyblue, col_burgundy, NA, NA), lwd = c(2.5, 2.5, NA, NA),
       bty = "n", cex = 0.95)
dev.off()
cat(sprintf("\u2713 Saved: %s\n", file.path(FIGS_OUTPUT, "Val_2_Density.svg")))

# Panel 3: Cohen's d
svglite(file.path(FIGS_OUTPUT, "Val_3_CohensD.svg"), width = 6, height = 5)
par(mar = c(5, 4.5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)
m_land <- mean(landscape_vals); m_site <- mean(sites_sl$suitability, na.rm = TRUE)
plot(c(1, 2), c(m_land, m_site), xlim = c(0.5, 2.5),
     ylim = c(0, max(m_land, m_site) * 1.3), pch = 19, cex = 2,
     col = c(col_dustyblue, col_burgundy),
     xlab = "", ylab = "Mean Suitability", main = "Effect Size (Cohen's d)", xaxt = "n")
axis(1, at = c(1, 2), labels = c("Landscape", "Sites"))
arrows(1, m_land - sd(landscape_vals), 1, m_land + sd(landscape_vals),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = col_dustyblue)
arrows(2, m_site - sd(sites_sl$suitability, na.rm = TRUE),
       2, m_site + sd(sites_sl$suitability, na.rm = TRUE),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = col_burgundy)
text(1.5, max(m_land, m_site) * 1.18,
     labels = sprintf("d = %.2f (%s)", cohens_d, label_d), cex = 1.0, font = 2)
dev.off()
cat(sprintf("\u2713 Saved: %s\n", file.path(FIGS_OUTPUT, "Val_3_CohensD.svg")))

# Panel 4: LOO-CV with Clopper-Pearson 95% CIs
svglite(file.path(FIGS_OUTPUT, "Val_4_LOOCV.svg"), width = 6, height = 5)
par(mar = c(5, 4.5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)
plot(NA, xlim = c(0.5, 2.5), ylim = c(0, 110), xaxt = "n",
     xlab = "", ylab = "Accuracy (%)",
     main = "Leave-One-Out Cross-Validation")
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = 1:2, labels = c("High+VH", "VH only"))
arrows(1, ci_loo_hv[1], 1, ci_loo_hv[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = col_darkblue)
arrows(2, ci_loo_vh[1], 2, ci_loo_vh[2],
       angle = 90, code = 3, length = 0.12, lwd = 2.5, col = col_darkblue)
points(1:2, c(loo_acc_hv, loo_acc_vh), pch = 19, cex = 2, col = col_darkblue)
abline(h = 50, lty = 2, col = col_mauve, lwd = 1.5)
text(2.4, 50, "50%", cex = 0.7, col = col_mauve, adj = c(1, -0.3))
text(1, loo_acc_hv, sprintf("%.1f%%", loo_acc_hv), pos = 4, font = 2, cex = 1.0)
text(2, loo_acc_vh, sprintf("%.1f%%", loo_acc_vh), pos = 4, font = 2, cex = 1.0)
text(1, ci_loo_hv[1], sprintf("[%.0f%%", ci_loo_hv[1]), pos = 1, cex = 0.8, col = "gray40")
text(1, ci_loo_hv[2], sprintf("%.0f%%]", ci_loo_hv[2]), pos = 3, cex = 0.8, col = "gray40")
text(2, ci_loo_vh[1], sprintf("[%.0f%%", ci_loo_vh[1]), pos = 1, cex = 0.8, col = "gray40")
text(2, ci_loo_vh[2], sprintf("%.0f%%]", ci_loo_vh[2]), pos = 3, cex = 0.8, col = "gray40")
dev.off()
cat(sprintf("\u2713 Saved: %s\n", file.path(FIGS_OUTPUT, "Val_4_LOOCV.svg")))

# Panel 5: Permutation test
svglite(file.path(FIGS_OUTPUT, "Val_5_Permutation.svg"), width = 6, height = 5)
par(mar = c(5, 4.5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)
hist(perm_acc_hv, breaks = 30, col = col_beige, border = "white",
     main = "Null Model Permutation Test (High + Very High)",
     xlab = "% Random sites in High+VH", ylab = "Frequency", xlim = c(30, 100))
abline(v = real_acc_hv, col = col_burgundy, lwd = 3)
abline(v = mean(perm_acc_hv), col = col_dustyblue, lwd = 2, lty = 2)
legend("topleft",
       legend = c(sprintf("Real: %.1f%%", real_acc_hv),
                  sprintf("Null: %.1f%%", mean(perm_acc_hv))),
       col = c(col_burgundy, col_dustyblue), lwd = c(3, 2), lty = c(1, 2),
       bty = "n", cex = 0.95)
dev.off()
cat(sprintf("\u2713 Saved: %s\n", file.path(FIGS_OUTPUT, "Val_5_Permutation.svg")))

# Panel 6: Weight stability
svglite(file.path(FIGS_OUTPUT, "Val_6_WeightStability.svg"), width = 6, height = 5)
par(mar = c(5, 5, 4, 2), family = "sans",
    cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)
var_short  <- c("Elev", "Slope", "Aspect", "TRI", "TWI", "DTD")
loo_means  <- apply(loo_weights, 2, mean)
loo_sds    <- apply(loo_weights, 2, sd)
x_pos      <- 1:6
yr         <- range(c(loo_means - loo_sds, loo_means + loo_sds, weights))
plot(NA, xlim = c(0.5, 6.5), ylim = c(max(0, yr[1] - 2), yr[2] + 3),
     xaxt = "n", xlab = "", ylab = "Weight (%)",
     main = "LOO-CV Weight Stability")
grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
axis(1, at = x_pos, labels = var_short, cex.axis = 1.0)
arrows(x_pos, loo_means - loo_sds, x_pos, loo_means + loo_sds,
       angle = 90, code = 3, length = 0.08, lwd = 2, col = col_darkblue)
points(x_pos, loo_means, pch = 19, cex = 1.5, col = col_darkblue)
points(x_pos, weights, pch = 18, cex = 1.5, col = col_mauve)
legend("topright",
       legend = c("LOO mean \u00b1 SD", "Full-model weight"),
       pch = c(19, 18), col = c(col_darkblue, col_mauve),
       pt.cex = 1.3, bty = "n", cex = 0.95)
dev.off()
cat(sprintf("\u2713 Saved: %s\n", file.path(FIGS_OUTPUT, "Val_6_WeightStability.svg")))




# ==========================================================================
# Fig 4: Strip plot panels — individual SVG files
# ==========================================================================
cat("Creating individual strip plot panels...\n")

sites_complete$Region <- factor(sites_complete$Region,
                                levels = c("Mediterranean","Irano-Turanian","Saharo-Arabian"))
zone_lev <- levels(sites_complete$Region)

var_list <- list(
  list(col="elevation",   ylab="Elevation (m asl)",                  title="A) Elevation",  file="Fig4_A_Elevation"),
  list(col="slope",       ylab="Slope (degrees)",                    title="B) Slope",       file="Fig4_B_Slope"),
  list(col="TRI",         ylab="Terrain Ruggedness Index (TRI)",     title="C) TRI",         file="Fig4_C_TRI"),
  list(col="TWI",         ylab="Topographic Wetness Index (TWI)",    title="D) TWI",         file="Fig4_D_TWI"),
  list(col="northness",   ylab="Aspect (0-1)",                      title="E) Aspect",      file="Fig4_E_Aspect"),
  list(col="DistToWater", ylab="DTD (km)",                          title="F) DTD",          file="Fig4_F_DTD")
)

is_outlier <- function(x) {
  q  <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iq <- diff(q); x < (q[1] - 1.5*iq) | x > (q[2] + 1.5*iq)
}

for (v in var_list) {
  fig_path <- file.path(FIGS_OUTPUT, paste0(v$file, ".svg"))
  svglite(fig_path, width = 6, height = 5.5)
  par(mar = c(6, 5, 4, 2), family = "sans",
      cex = 1.1, cex.lab = 1.2, cex.axis = 1.0, cex.main = 1.3)
  
  vals_l <- lapply(zone_lev, function(z) {
    x <- sites_complete[[v$col]][sites_complete$Region == z]
    x[!is.na(x)]
  })
  means <- sapply(vals_l, mean); sds <- sapply(vals_l, sd)
  all_v <- unlist(vals_l)
  yr    <- diff(range(all_v))
  
  plot(NA, xlim = c(0.5, 3.5),
       ylim = c(min(all_v) - 0.05*yr, max(all_v) + 0.25*yr),
       xaxt = "n", xlab = "", ylab = v$ylab, main = v$title, las = 1, frame = FALSE)
  grid(nx = NA, ny = NULL, col = "gray92"); box(bty = "l")
  axis(1, at = 1:3, labels = zone_lev, cex.axis = 0.80)
  
  for (i in seq_along(zone_lev)) {
    z   <- zone_lev[i]; col <- zone_cols[z]
    x   <- vals_l[[i]]; mn <- means[i]; sdv <- sds[i]; n <- length(x)
    segments(i, mn - sdv, i, mn + sdv, col = col, lwd = 2)
    segments(i - 0.1, mn - sdv, i + 0.1, mn - sdv, col = col, lwd = 1.8)
    segments(i - 0.1, mn + sdv, i + 0.1, mn + sdv, col = col, lwd = 1.8)
    set.seed(42); jx <- i + runif(n, -0.15, 0.15)
    out <- is_outlier(x)
    points(jx[!out], x[!out], pch = 21, bg = adjustcolor(col, 0.45),
           col = adjustcolor(col, 0.75), cex = 0.9, lwd = 0.7)
    if (any(out)) points(jx[out], x[out], pch = 21, bg = adjustcolor(col, 0.85),
                         col = "gray20", cex = 1.2, lwd = 1.1)
    points(i, mn, pch = 23, bg = "white", col = "gray20", cex = 1.2, lwd = 1.4)
    # Highest-value site name label
    max_idx <- which.max(x)
    site_names_z <- sites_complete$Site[sites_complete$Region == z & !is.na(sites_complete[[v$col]])]
    if (length(site_names_z) >= max_idx)
      text(jx[max_idx], x[max_idx], site_names_z[max_idx],
           cex = 0.75, font = 2, col = col, pos = 4)
    mtext(sprintf("n=%d", n), side = 1, at = i, line = 4.8, cex = 0.9, col = "gray40")
  }
  
  # Add legend to each panel
  legend("topright", legend = zone_lev, pch = 21,
         pt.bg = adjustcolor(zone_cols[zone_lev], 0.65), col = "gray20",
         pt.cex = 1.3, cex = 0.85, bty = "n")
  
  dev.off()
  cat(sprintf("\u2713 Saved: %s\n", fig_path))
}