# --- SETUP AND LIBRARIES ---
rm(list=ls())
gc()

# Define your base directory (MUST be the Linux path)
dir <- "C:/Users/steph/Documents/SardineTMB"

load.lib <- c("ncdf4", "tidyverse", "data.table", "sf", "tools")
install.lib <- load.lib[!load.lib %in% installed.packages()]

for(lib in install.lib) install.packages(lib, dependencies = TRUE)
sapply(load.lib, require, character = TRUE)

# Turn on s2 and use geos since you are computing distance and area
sf_use_s2(TRUE)
# --- END OF SETUP ---

# --- Create CRS object ---
aea_custom <- "+proj=aea +lat_1=30 +lat_2=50 +lat_0=40 +lon_0=-124.5 +datum=WGS84 +units=m +no_defs"
target_crs <- st_crs(aea_custom)

proj_string <- target_crs$proj4string  # The PROJ version
wkt_string  <- target_crs$wkt          # The WKT version
# --- End create CRS object ---

# Create directory to save downloaded ROMS dataset
ROMSpath <- "Data/ROMS2/Regridded"
if (!dir.exists(ROMSpath)) {
  dir.create(ROMSpath, recursive = TRUE)
}

# Input ROMS files
fls <- list.files(file.path(dir, ROMSpath), pattern = "\\.nc$", full.names = TRUE)

# --- Check structure ---
res.check <- fls[1]
res.check <- nc_open(res.check)

lon.ck <- ncvar_get(res.check, "x")
lat.ck <- ncvar_get(res.check, "y")

lon.ck <- t(lon.ck)
lat.ck <- t(lat.ck)

dx <- abs(lon.ck[, -1] - lon.ck[, -ncol(lon.ck)])
dy <- abs(lat.ck[-1, ] - lat.ck[-nrow(lat.ck), ])

mean_dx <- mean(dx, na.rm = TRUE)
mean_dy <- mean(dy, na.rm = TRUE)

cat("Estimated average resolution:\n")
cat("Longitude:", round(mean_dx, 5), "°\n")
cat("Latitude :", round(mean_dy, 5), "°\n")

lon_range <- c(min(lon.ck), max(lon.ck))
lat_range <- c(min(lat.ck), max(lon.ck))

bbox_aea <- st_bbox(c(xmin = lon_range[1], xmax = lon_range[2],
                      ymin = lat_range[1], ymax = lat_range[2]), 
                      crs = aea_custom)

bbox_proj <- st_as_sfc(bbox_aea) |>
  st_transform(crs = st_crs(4326)) |>
  st_bbox()

min_x <- bbox_proj["xmin"]
max_x <- bbox_proj["xmax"]
min_y <- bbox_proj["ymin"]
max_y <- bbox_proj["ymax"]

nc_close(res.check)