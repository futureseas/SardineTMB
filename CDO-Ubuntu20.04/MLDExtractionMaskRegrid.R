# --- SETUP AND LIBRARIES ---
rm(list=ls())
gc()

# Base directory (Linux/WSL path)
dir <- "/mnt/c/Users/steph/Documents/SardineTMB"

load.lib <- c("ncdf4", "tidyverse", "data.table", "sf", "tools")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies = TRUE)
sapply(load.lib, require, character = TRUE)

sf_use_s2(TRUE)
# --- END OF SETUP ---

# --- CREATE CRS OBJECT ---
aea_custom <- "+proj=aea +lat_1=30 +lat_2=50 +lat_0=40 +lon_0=-124.5 +datum=WGS84 +units=m +no_defs"
target_crs <- st_crs(aea_custom)

proj_string <- target_crs$proj4string
wkt_string  <- target_crs$wkt
# --- END CREATE CRS OBJECT ---

# Input directory check
ROMSpath <- "Data/ROMS2.testMixedInput"
if (!dir.exists(file.path(dir, ROMSpath))) {
  dir.create(file.path(dir, ROMSpath), recursive = TRUE)
}

fls <- list.files(file.path(dir, ROMSpath), pattern = "\\.nc$", full.names = TRUE)

# --- CHECK STRUCTURE ---
if (length(fls) == 0) stop("No NetCDF files found in specified input directory.")

res.check <- nc_open(fls[1])
lon.ck <- t(ncvar_get(res.check, "longitude"))
lat.ck <- t(ncvar_get(res.check, "latitude"))

dx <- abs(lon.ck[, -1] - lon.ck[, -ncol(lon.ck)])
dy <- abs(lat.ck[-1, ] - lat.ck[-nrow(lat.ck), ])

mean_dx <- mean(dx, na.rm = TRUE)
mean_dy <- mean(dy, na.rm = TRUE)

cat("Estimated average resolution:\n")
cat("Longitude:", round(mean_dx, 5), "°\n")
cat("Latitude :", round(mean_dy, 5), "°\n")

lon_range  <- c(min(lon.ck), max(lon.ck))
lat_range  <- c(min(lat.ck), max(lat.ck))

bbox_wgs84 <- st_bbox(c(xmin = lon_range[1], xmax = lon_range[2],
                        ymin = lat_range[1], ymax = lat_range[2]), 
                      crs = st_crs(4326))

bbox_proj <- st_as_sfc(bbox_wgs84) |>
  st_transform(crs = aea_custom) |>
  st_bbox()

min_x <- bbox_proj["xmin"]
max_x <- bbox_proj["xmax"]
min_y <- bbox_proj["ymin"]
max_y <- bbox_proj["ymax"]

nc_close(res.check)
# --- END STRUCTURE CHECK ---

# --- APPLY MASK (MLD ONLY) ---
fls <- fls[!grepl("bathy", basename(fls), ignore.case = TRUE)]

masked_subfolder <- file.path(dir, "Data", "ROMS2", "Masked")
if (!dir.exists(masked_subfolder)) dir.create(masked_subfolder, recursive = TRUE)

for (fl in fls) {
  
  message("Processing file: ", basename(fl))
  
  # Extract 4-digit year from input filename
  yr <- str_extract(basename(fl), "\\d{4}")
  if (is.na(yr)) yr <- "unknown"
  
  nc_in <- nc_open(fl)
  
  # Focus strictly on mld
  var_names <- intersect(names(nc_in$var), "mld")
  if (length(var_names) == 0) {
    message("  Skipping: 'mld' variable not found in ", basename(fl))
    nc_close(nc_in)
    next
  }
  
  mask <- ncvar_get(nc_in, "mask")
  lon  <- ncvar_get(nc_in, "longitude")
  lat  <- ncvar_get(nc_in, "latitude")
  time <- ncvar_get(nc_in, "time")
  
  for (var_name in var_names) {
    message("  Masking variable: ", var_name)
    
    var_data <- ncvar_get(nc_in, var_name)
    var_data <- drop(var_data)
    
    for (t in seq_len(dim(var_data)[3])) {
      var_data[,,t][mask == 0] <- NA
    }
    
    dim_time <- ncdim_def("time", "seconds since 1900-01-01 00:00:00", vals = time, unlim = TRUE)
    dim_lat  <- ncdim_def("eta_rho", "degrees_north", vals = seq_len(dim(lat)[2]))
    dim_lon  <- ncdim_def("xi_rho", "degrees_east", vals = seq_len(dim(lon)[1]))
    
    var_def <- ncvar_def(
      var_name,
      units = tryCatch(nc_in$var[[var_name]]$units, error = function(e) "m"),
      dim = list(dim_lon, dim_lat, dim_time),
      missval = NA,
      prec = "float"
    )
    
    dim_nbnds    <- ncdim_def("nbnds", "", 1:4)
    lon_def      <- ncvar_def("lon", "degrees_east", list(dim_lon, dim_lat), prec = "float")
    lat_def      <- ncvar_def("lat", "degrees_north", list(dim_lon, dim_lat), prec = "float")
    lon_bnds_def <- ncvar_def("lon_bnds", "degrees_east", list(dim_nbnds, dim_lon, dim_lat), prec = "double")
    lat_bnds_def <- ncvar_def("lat_bnds", "degrees_north", list(dim_nbnds, dim_lon, dim_lat), prec = "double")
    
    # Rename file output structure to mld_hindcast_<year>_masked.nc
    masked_fname <- paste0("mld_hindcast_", yr, "_masked.nc")
    masked_full  <- file.path(masked_subfolder, masked_fname)
    
    nc_out <- nc_create(masked_full, list(var_def, lon_def, lat_def, lon_bnds_def, lat_bnds_def))
    
    lon_bnds <- array(NA, dim = c(4, dim(lon)[1], dim(lon)[2]))
    lat_bnds <- array(NA, dim = c(4, dim(lat)[1], dim(lat)[2]))
    
    calc_delta <- function(m) {
      d <- matrix(0, nrow=nrow(m), ncol=ncol(m))
      d[1:(nrow(m)-1), ] <- diff(m) 
      d[nrow(m), ] <- d[nrow(m)-1, ]
      return(d * 0.5)
    }
    
    d_lon <- calc_delta(lon)
    d_lat <- t(calc_delta(t(lat)))
    
    lon_bnds[1,,] <- lon - d_lon; lat_bnds[1,,] <- lat - d_lat
    lon_bnds[2,,] <- lon + d_lon; lat_bnds[2,,] <- lat - d_lat
    lon_bnds[3,,] <- lon + d_lon; lat_bnds[3,,] <- lat + d_lat
    lon_bnds[4,,] <- lon - d_lon; lat_bnds[4,,] <- lat + d_lat
    
    ncvar_put(nc_out, var_name, var_data)
    ncvar_put(nc_out, "lon", lon)
    ncvar_put(nc_out, "lat", lat)
    ncvar_put(nc_out, "lon_bnds", lon_bnds)
    ncvar_put(nc_out, "lat_bnds", lat_bnds)
    
    # Variable Attributes
    ncatt_put(nc_out, "mld", "long_name", "Mixed layer depth")
    ncatt_put(nc_out, "mld", "units", "m")
    ncatt_put(nc_out, "mld", "standard_name", "ocean_mixed_layer_thickness_defined_by_sigma_theta")
    ncatt_put(nc_out, "mld", "coordinates", "lon lat")
    
    # Coordinate Metadata
    ncatt_put(nc_out, "lon", "standard_name", "longitude")
    ncatt_put(nc_out, "lon", "long_name", "longitude of RHO points")
    ncatt_put(nc_out, "lon", "units", "degrees_east")
    ncatt_put(nc_out, "lon", "axis", "X")
    ncatt_put(nc_out, "lon", "bounds", "lon_bnds")
    
    ncatt_put(nc_out, "lat", "standard_name", "latitude")
    ncatt_put(nc_out, "lat", "long_name", "latitude of RHO points")
    ncatt_put(nc_out, "lat", "units", "degrees_north")
    ncatt_put(nc_out, "lat", "axis", "Y")
    ncatt_put(nc_out, "lat", "bounds", "lat_bnds")
    
    ncatt_put(nc_out, "lon_bnds", "long_name", "Longitude cell boundary corner coordinates")
    ncatt_put(nc_out, "lat_bnds", "long_name", "Latitude cell boundary corner coordinates")
    
    # Global Attributes
    ncatt_put(nc_out, 0, "Conventions", "CF-1.6, ACDD-1.3")
    ncatt_put(nc_out, 0, "title", "Subsetted and Masked ROMS-NEP Surface Data - MLD")
    ncatt_put(nc_out, 0, "history", paste("Created", Sys.time()))
    ncatt_put(nc_out, 0, "author", "Stephanie Hopkins, Mer Pozo Buil, Jerome Fiechter")
    ncatt_put(nc_out, 0, "creator_name", "Stephanie Hopkins")
    ncatt_put(nc_out, 0, "creator_email", "shopkin1@ucsc