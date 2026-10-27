# --- SETUP AND LIBRARIES ---
rm(list=ls())
gc()

# Define your base directory (MUST be the Linux path)
dir <- "/mnt/c/Users/steph/Documents/SardineTMB"

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
ROMSpath <- "Data/ROMS2"
if (!dir.exists(ROMSpath)) {
  dir.create(ROMSpath, recursive = TRUE)
}

# Input ROMS files
fls <- list.files(file.path(dir, ROMSpath), pattern = "\\.nc$", full.names = TRUE)

# --- Check structure ---
res.check <- fls[1]
res.check <- nc_open(res.check)

lon.ck <- ncvar_get(res.check, "longitude")
lat.ck <- ncvar_get(res.check, "latitude")

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
lat_range <- c(min(lat.ck), max(lat.ck))

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
# --- End structure check ---

# --- Apply mask to each variable and time step ---
masked.files.is.made = TRUE
if(masked.files.is.made == FALSE) {

  fls <- list.files(file.path(dir, ROMSpath), pattern = "\\.nc$", full.names = TRUE)
  fls <- fls[!grepl("bathy", basename(fls), ignore.case = TRUE)]

  if (length(fls) == 0) stop("No non-bathymetry NetCDF files found in specified directory.")

  masked_subfolder <- file.path(dir, "Data", "ROMS2", "Masked")
  if (!dir.exists(masked_subfolder)) dir.create(masked_subfolder, recursive = TRUE)

  for (fl in fls) {
    
    message("Processing file: ", basename(fl))

    # Extract 4-digit year from input filename
    yr_match <- regmatches(basename(fl), regexpr("(19|20)\\d{2}", basename(fl)))
    yr <- if (length(yr_match) > 0) yr_match else NA

    nc_temp <- nc_open(fl, write = TRUE)
    if ("u_vagv" %in% names(nc_temp$var)) {
        ncvar_rename(nc_temp, "u_vagv", "u_vavg")
        message("  Renamed 'u_vagv' to 'u_vavg'")
    }
    nc_close(nc_temp)

    nc_in <- nc_open(fl)
   
    mask <- ncvar_get(nc_in, "mask")
    lon  <- ncvar_get(nc_in, "longitude")
    lat  <- ncvar_get(nc_in, "latitude")
    time <- ncvar_get(nc_in, "time")
    
    var_names <- setdiff(names(nc_in$var), c("longitude", "latitude", "level", "time", "mask", "dz"))
    
    for (var_name in var_names) {
      message("  Masking variable: ", var_name)
      
      var_data <- ncvar_get(nc_in, var_name)
      var_data <- drop(var_data)
      
      for (t in seq_len(dim(var_data)[3])) {
        var_data[,,t][mask == 0] <- NA
      }
      
      # Define dimensions
      dim_time <- ncdim_def("time", "seconds since 1900-01-01 00:00:00", vals = time, unlim = TRUE)
      dim_lat  <- ncdim_def("eta_rho", "1", vals = seq_len(dim(lat)[2]))
      dim_lon  <- ncdim_def("xi_rho", "1", vals = seq_len(dim(lon)[1]))
      dim_nbnds <- ncdim_def("nbnds", "1", 1:4)

      var_def <- ncvar_def(
        var_name,
        units = nc_in$var[[var_name]]$units,
        dim = list(dim_lon, dim_lat, dim_time),
        missval = NA,
        prec = "float"
      )
      
      lon_def <- ncvar_def("lon", "degrees_east", list(dim_lon, dim_lat), prec = "float")
      lat_def <- ncvar_def("lat", "degrees_north", list(dim_lon, dim_lat), prec = "float")
      
      lon_bnds_def <- ncvar_def("lon_bnds", "degrees_east", list(dim_nbnds, dim_lon, dim_lat), prec = "double")
      lat_bnds_def <- ncvar_def("lat_bnds", "degrees_north", list(dim_nbnds, dim_lon, dim_lat), prec = "double")

      # Construct cleaned output filename: VAR_hindcast_YYYY_masked.nc
      if (!is.na(yr)) {
        masked_fname <- paste0(var_name, "_hindcast_", yr, "_masked.nc")
      } else {
        masked_fname <- paste0(var_name, "_masked.nc")
      }
      masked_full <- file.path(masked_subfolder, masked_fname)
      
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
      d_lat <- calc_delta(t(lat))
      d_lat <- t(d_lat)

      lon_bnds[1,,] <- lon - d_lon; lat_bnds[1,,] <- lat - d_lat
      lon_bnds[2,,] <- lon + d_lon; lat_bnds[2,,] <- lat - d_lat
      lon_bnds[3,,] <- lon + d_lon; lat_bnds[3,,] <- lat + d_lat
      lon_bnds[4,,] <- lon - d_lon; lat_bnds[4,,] <- lat + d_lat

      # Write data
      ncvar_put(nc_out, var_name, var_data)
      ncvar_put(nc_out, "lon", lon)
      ncvar_put(nc_out, "lat", lat)
      ncvar_put(nc_out, "lon_bnds", lon_bnds)
      ncvar_put(nc_out, "lat_bnds", lat_bnds)

      # --- Variable Specific Metadata ---
      if (var_name == "mld") {
        ncatt_put(nc_out, "mld", "long_name", "Mixed layer depth")
        ncatt_put(nc_out, "mld", "units", "m")
        ncatt_put(nc_out, "mld", "standard_name", "ocean_mixed_layer_thickness_defined_by_sigma_theta")
      } else {
        var_long <- tryCatch(nc_in$var[[var_name]]$longname, error = function(e) var_name)
        ncatt_put(nc_out, var_name, "long_name", var_long)
      }
      
      # --- Dimension Descriptive Metadata ---
      ncatt_put(nc_out, "xi_rho", "long_name", "x-dimension index on RHO grid")
      ncatt_put(nc_out, "eta_rho", "long_name", "y-dimension index on RHO grid")
      ncatt_put(nc_out, "nbnds", "long_name", "number of cell bounding box vertices")
      
      # --- Define Coordinate Metadata ---
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
      
      # --- Boundary Variables Metadata ---
      ncatt_put(nc_out, "lon_bnds", "long_name", "Longitude cell boundary corner coordinates")
      ncatt_put(nc_out, "lat_bnds", "long_name", "Latitude cell boundary corner coordinates")
     
      # --- Link Data Variable to Coordinates ---
      ncatt_put(nc_out, var_name, "coordinates", "lon lat")

      # --- Global Attributes (CF-1.6 & ACDD) ---
      ncatt_put(nc_out, 0, "Conventions", "CF-1.6, ACDD-1.3")
      ncatt_put(nc_out, 0, "title", "Subsetted and Masked ROMS-NEP Surface Data")
      ncatt_put(nc_out, 0, "history", paste("Created", Sys.time()))
      
      # Authorship & Attribution
      ncatt_put(nc_out, 0, "author", "Stephanie Hopkins, Mer Pozo Buil, Jerome Fiechter")
      ncatt_put(nc_out, 0, "creator_name", "Stephanie Hopkins")
      ncatt_put(nc_out, 0, "creator_email", "shopkin1@ucsc.edu")
      ncatt_put(nc_out, 0, "institution", "University of California Santa Cruz / NOAA SWFSC / NOAA NWFSC")
      
      # Provenance & Methods
      ncatt_put(nc_out, 0, "source", "ROMS-NEP (Northeast Pacific) regional ocean model runs provided by Mer Pozo Buil")
      ncatt_put(nc_out, 0, "references", "Lazier et al. (2002); Drenkard et al. (2025)")
      ncatt_put(nc_out, 0, "mld_computation", "Mixed Layer Depth computed using potential density with a fixed threshold criterion (0.03 kg/m3 difference from surface).")

      nc_close(nc_out)
    }
    
    nc_close(nc_in)
  }

  message("All files processed and masked!")
}
# --- End mask function ---

# --- Create empty grid ---
target_grid_nc  <- file.path(dir, ROMSpath, "target_10km_grid.nc")
target_grid_txt <- file.path(dir, ROMSpath, "target_10km_grid.txt")

x_vals <- seq(from = floor(min_x/10000)*10000, to = ceiling(max_x/10000)*10000, by = 10000)
y_vals <- seq(from = floor(min_y/10000)*10000, to = ceiling(max_y/10000)*10000, by = 10000)

grid_df     <- expand.grid(x = x_vals, y = y_vals)
grid_sf     <- st_as_sf(grid_df, coords = c("x", "y"), crs = target_crs)
grid_wgs84  <- st_transform(grid_sf, 4326)
coords      <- st_coordinates(grid_wgs84)

lon_2d <- matrix(coords[,1], nrow = length(x_vals), ncol = length(y_vals))
lat_2d <- matrix(coords[,2], nrow = length(x_vals), ncol = length(y_vals))

dx <- c(-5000, 5000, 5000, -5000)
dy <- c(-5000, -5000, 5000, 5000)

lon_corners <- array(NA, dim = c(4, length(x_vals), length(y_vals)))
lat_corners <- array(NA, dim = c(4, length(x_vals), length(y_vals)))

all_corners_df <- do.call(rbind, lapply(1:4, function(i) {
  expand.grid(x = x_vals + dx[i], y = y_vals + dy[i])
}))

all_corners_sf    <- st_as_sf(all_corners_df, coords = c("x", "y"), crs = target_crs)
all_corners_wgs84 <- st_transform(all_corners_sf, 4326)
all_coords        <- st_coordinates(all_corners_wgs84)

n_points <- length(x_vals) * length(y_vals)

for (i in 1:4) {
  start_idx <- ((i - 1) * n_points) + 1
  end_idx   <- i * n_points

  lon_corners[i,,] <- matrix(all_coords[start_idx:end_idx, 1], nrow = length(x_vals), ncol = length(y_vals))
  lat_corners[i,,] <- matrix(all_coords[start_idx:end_idx, 2], nrow = length(x_vals), ncol = length(y_vals))
}

dim_x  <- ncdim_def("x", "m", x_vals)
dim_y  <- ncdim_def("y", "m", y_vals)
dim_nv <- ncdim_def("nv", "", 1:4, create_dimvar = FALSE)

var_lon      <- ncvar_def("lon", "degrees_east", list(dim_x, dim_y), prec="double")
var_lat      <- ncvar_def("lat", "degrees_north", list(dim_x, dim_y), prec="double")
var_lon_bnds <- ncvar_def("lon_bnds", "degrees_east", list(dim_nv, dim_x, dim_y), prec="double")
var_lat_bnds <- ncvar_def("lat_bnds", "degrees_north", list(dim_nv, dim_x, dim_y), prec="double")
var_crs      <- ncvar_def("crs", "", list(), prec="integer")
var_dummy    <- ncvar_def("dummy", "1", list(dim_x, dim_y), prec="double")

nc <- nc_create(target_grid_nc, list(var_lon, var_lat, var_lon_bnds, var_lat_bnds, var_crs, var_dummy))

ncvar_put(nc, var_lon, lon_2d)
ncvar_put(nc, var_lat, lat_2d)
ncvar_put(nc, var_lon_bnds, lon_corners)
ncvar_put(nc, var_lat_bnds, lat_corners)
ncvar_put(nc, var_dummy, matrix(0, nrow=length(x_vals), ncol=length(y_vals)))

ncatt_put(nc, "lon", "standard_name", "longitude")
ncatt_put(nc, "lon", "units", "degrees_east")
ncatt_put(nc, "lon", "axis", "X")
ncatt_put(nc, "lon", "bounds", "lon_bnds")

ncatt_put(nc, "lat", "standard_name", "latitude")
ncatt_put(nc, "lat", "units", "degrees_north")
ncatt_put(nc, "lat", "axis", "Y")
ncatt_put(nc, "lat", "bounds", "lat_bnds")

ncatt_put(nc, "dummy", "coordinates", "lon lat")
ncatt_put(nc, "dummy", "grid_mapping", "crs")

ncatt_put(nc, "crs", "grid_mapping_name", "albers_conical_equal_area")
ncatt_put(nc, "crs", "standard_parallel", c(30, 50), prec = "double")
ncatt_put(nc, "crs", "longitude_of_central_meridian", -124.5)
ncatt_put(nc, "crs", "latitude_of_projection_origin", 40)
ncatt_put(nc, "crs", "false_easting", 0.0)
ncatt_put(nc, "crs", "false_northing", 0.0)
ncatt_put(nc, "crs", "crs_wkt", wkt_string)
ncatt_put(nc, "crs", "proj4text", proj_string)

ncatt_put(nc, "x", "standard_name", "projection_x_coordinate")
ncatt_put(nc, "x", "long_name", "Albers Equal Area Projected X coordinate")
ncatt_put(nc, "x", "axis", "X")
ncatt_put(nc, "y", "standard_name", "projection_y_coordinate")
ncatt_put(nc, "y", "long_name", "Albers Equal Area Projected Y coordinate")
ncatt_put(nc, "y", "axis", "Y")

ncatt_put(nc, 0, "crs_wkt", wkt_string)
ncatt_put(nc, 0, "Conventions", "CF-1.6")

nc_close(nc)

cdo_cmd <- sprintf("cdo griddes %s > %s", target_grid_nc, target_grid_txt)
if (.Platform$OS.type == "windows") {
  cdo_cmd <- paste("wsl", cdo_cmd)
}
system(cdo_cmd)
# --- End grid creation ---

# --- Regrid data ---
regrid_subfolder <- file.path(dir, ROMSpath, "Regridded")
if (!dir.exists(regrid_subfolder)) dir.create(regrid_subfolder, recursive = TRUE)

# Path helper for WSL compatibility
to_wsl_path <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  sub("^([A-Za-z]):", "/mnt/\\L\\1", p, perl = TRUE)
}

# Use the grid text description for CDO remapping
linux_target_grid <- file.path(dir, ROMSpath, "target_10km_grid.txt")

fls.msk <- list.files(file.path(dir, ROMSpath, "Masked"), pattern = "\\.nc$", full.names = TRUE)

for (f in fls.msk){
  
  var_file <- tools::file_path_sans_ext(basename(f))
  
  if (grepl("target_grid|target_10km_grid|dz_masked|bathy", var_file, ignore.case = TRUE)) {
    next
  }
  
  # 1. Remove leading 'allvars_' prefix if present
  cleaned <- sub("^allvars_", "", var_file)
  
  # 2. Remove trailing '_masked' suffix if present
  cleaned <- sub("_masked$", "", cleaned)

  # 3. Extract var_name by removing '_hindcast' and any following year/suffix
  internal_var <- sub("_(hindcast|\\d{4}).*", "", cleaned)

  if (is.na(internal_var) || internal_var == "") {
    warning(paste("Could not extract internal_var for file:", var_file, "- Skipping."))
    next
  }

  linux_out_file <- file.path(regrid_subfolder, paste0(var_file, "_regrid.nc"))
  
  method <- if (grepl("chla_vint|microzooplankton_vint|mesozooplankton_vint|Pzooplankton_vint|u_vavg|v_vavg|w_vavg",
                      internal_var, ignore.case = TRUE)) {
    "-remapcon"
  } else {
    "-remapbil"
  }
  
  target_grid_txt_cmd <- target_grid_txt
  f_cmd               <- f
  linux_out_file_cmd  <- linux_out_file

  if (.Platform$OS.type == "windows") {
    target_grid_txt_cmd <- system(sprintf("wsl wslpath -a '%s'", target_grid_txt), intern = TRUE)
    f_cmd               <- system(sprintf("wsl wslpath -a '%s'", f), intern = TRUE)
    linux_out_file_cmd  <- system(sprintf("wsl wslpath -a '%s'", linux_out_file), intern = TRUE)
  }

  cmd <- sprintf(
    "cdo -O %s,%s -setgridtype,curvilinear -selname,%s %s %s",
    method,
    target_grid_txt_cmd,
    internal_var,
    f_cmd,
    linux_out_file_cmd
  )

  if (.Platform$OS.type == "windows") {
    cmd <- paste("wsl", cmd)
  }

  message("Processing: ", var_file, " -> Variable: ", internal_var)
  system(cmd)
}
message("All files processed and regridded!")
# ---- End Regridding ---

# --- Create validation table ---
validation_results <- list()

for (f in fls.msk) {
    var_name <- tools::file_path_sans_ext(basename(f))
    linux_out_file <- file.path(regrid_subfolder, paste0(var_name, "_regrid.nc"))
    if (file.exists(linux_out_file)) {
      
      f_cmd <- f
      out_cmd <- linux_out_file
      
      if (.Platform$OS.type == "windows") {
        f_cmd <- system(sprintf("wsl wslpath -a '%s'", f), intern = TRUE)
        out_cmd <- system(sprintf("wsl wslpath -a '%s'", linux_out_file), intern = TRUE)
      }
      
      orig_raw <- system(paste0("cdo -s outputf,%g -fldmean ", f_cmd, " 2>/dev/null"), intern = TRUE)
      orig_mean_vec <- suppressWarnings(as.numeric(orig_raw))
      orig_mean_vec <- orig_mean_vec[!is.na(orig_mean_vec)]
      
      regrid_raw <- system(paste0("cdo -s outputf,%g -fldmean ", out_cmd, " 2>/dev/null"), intern = TRUE)
      regrid_mean_vec <- suppressWarnings(as.numeric(regrid_raw))
      regrid_mean_vec <- regrid_mean_vec[!is.na(regrid_mean_vec)]
      
      orig_mean   <- mean(orig_mean_vec, na.rm = TRUE)
      regrid_mean <- mean(regrid_mean_vec, na.rm = TRUE)
      
      pct_change <- if (!is.na(orig_mean) && orig_mean != 0) {
        ((regrid_mean - orig_mean) / orig_mean) * 100
      } else {
        NA_real_
      }
      
      status_check <- if (is.na(pct_change) || is.nan(pct_change) || is.infinite(pct_change)) {
        "NaN_ERROR"
      } else if (pct_change >= -10 && pct_change <= 10) {
        "PASS"
      } else {
        "CHECK"
      }
      
      validation_results[[var_name]] <- data.frame(
        Variable       = var_name,
        Original_Mean  = round(orig_mean, 4),
        Regridded_Mean = round(regrid_mean, 4),
        Pct_Change     = round(pct_change, 3),
        Status         = status_check,
        stringsAsFactors = FALSE
      )
    } else {
      warning("Regridded file missing for: ", var_name)
    }
}

if (length(validation_results) > 0) {
  final_table <- do.call(rbind, validation_results)
  print(final_table)
  write.csv(final_table, file.path(regrid_subfolder, "regridding_validation_report.csv"), row.names = FALSE)
} else {
  warning("No valid results were generated for the validation report.")
}
# --- End creation of validation table ---
