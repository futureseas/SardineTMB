# --- SETUP AND LIBRARIES ---
# The R Markdown chunk delimiters are REMOVED (```{r ...})

rm(list=ls())
gc()

# Define your base directory (MUST be the Linux path)
dir <- "/mnt/c/Users/steph/Documents/Sardine_SptTM"

# Note the package "here" does not play well with linux
load.lib <- c("ncdf4", "tidyverse", "data.table", "sf", "tools"#,
              )

install.lib <- load.lib[!load.lib %in% installed.packages()]

for(lib in install.lib) install.packages(lib, dependencies = TRUE)
sapply(load.lib,require,character = TRUE)

# Turn on s2 and use geos since you are computing distance and area
sf_use_s2(TRUE)
# --- END OF SETUP ---

# --- Create CRS object ---
# Since the data falls within two UTM zones and there is not an ideal projection to use, apply a custom equal area Albers projection
aea_custom <- "+proj=aea +lat_1=30 +lat_2=50 +lat_0=40 +lon_0=-124.5 +datum=WGS84 +units=m +no_defs"

## Create a CRS object
target_crs <- st_crs(aea_custom)

## Extract the two versions
proj_string <- target_crs$proj4string  # The PROJ version
wkt_string <- target_crs$wkt                 # The WKT version (long and complex)
# --- End create CRS object

# Create directory to save downloaded ROMS dataset
ROMSpath <- "Data/ROMS2.testMixedInput"
  
  if (!dir.exists(ROMSpath)) {
  dir.create(ROMSpath, recursive = TRUE)
    }

# Input ROMS files
fls <- list.files(file.path(dir, ROMSpath), pattern = "\\.nc$", full.names = TRUE)

###-----------------------------------------------------------------
# Note mld and bathymetry data has been already subset with an
# elevation mask applied. 
### If working with the published data, skip to line 378!
### ----------------------------------------------------------------

# --- Check structure ---
# Open file
res.check <- fls[1]

res.check <- nc_open(res.check)

# Load 2D lon/lat
lon.ck <- ncvar_get(res.check, "longitude")  # [lat, lon]
lat.ck <- ncvar_get(res.check, "latitude")   # [lat, lon]
  
# Transpose if needed (see earlier logic)
lon.ck  <- t(lon.ck)
lat.ck <- t(lat.ck)

# Calculate differences between adjacent grid cells
# (differences in degrees, approx only)
dx <- abs(lon.ck[, -1] - lon.ck[, -ncol(lon.ck)])  # difference in longitude
dy <- abs(lat.ck[-1, ] - lat.ck[-nrow(lat.ck), ])  # difference in latitude

# Average resolution
mean_dx <- mean(dx, na.rm = TRUE)
mean_dy <- mean(dy, na.rm = TRUE)

cat("Estimated average resolution:\n")
cat("Longitude:", round(mean_dx, 5), "°\n")
cat("Latitude :", round(mean_dy, 5), "°\n")

# Create a bounding box for projections
## 1. Define source bounds 
lon_range <- c(min(lon.ck), max(lon.ck))
lat_range <- c(min(lat.ck), max(lat.ck))

## 2. Creat a bounding box to polygon
bbox_wgs84 <- st_bbox(c(xmin = lon_range[1], xmax = lon_range[2],
			ymin = lat_range[1], ymax = lat_range[2]), 
			crs = st_crs(4326))

## 3. Transform bounding box to target CRS
bbox_proj <- st_as_sfc(bbox_wgs84) |>
                st_transform(crs = aea_custom) |>
                st_bbox()

## 4. Extract the new bounds in Metres
min_x <- bbox_proj["xmin"]
max_x <- bbox_proj["xmax"]
min_y <- bbox_proj["ymin"]
max_y <- bbox_proj["ymax"]

nc_close(res.check)
# --- End structure check ---

# --- Apply mask to each variable and time step ---
# Input ROMS files
fls <- list.files(file.path(dir, "Data", "ROMS2"), pattern = "\\.nc$", full.names = TRUE)

# Exclude bathyemetry file
# EXPLICIT FILTER: Remove any file path that contains the word "bathy"
fls <- fls[!grepl("bathy", basename(fls), ignore.case = TRUE)]

# Subfolder for masked output
masked_subfolder <- file.path(dir, "Data", "ROMS2", "Masked")
if (!dir.exists(masked_subfolder)) dir.create(masked_subfolder, recursive = TRUE)

for (fl in fls) {
  
  message("Processing file: ", basename(fl))

  # Note their is a misspelling in one of the variable names:
  ## Open, Rename, and CLOSE immediately to force the save to disk
  nc_temp <- nc_open(fl, write = TRUE)
  
  # Check if the old name exists before trying to rename (prevents errors on 2nd run)
  if ("u_vagv" %in% names(nc_temp$var)) {
      ncvar_rename(nc_temp, "u_vagv", "u_vavg")
      message("  Renamed 'u_vagv' to 'u_vavg'")
  }
  nc_close(nc_temp) # This "flushes" the change to the file header

  # Open again with renamed variable
  nc_in <- nc_open(fl)
 
  # Load mask and coordinates
  mask <- ncvar_get(nc_in, "mask")          # [xi_rho x eta_rho x level x time]
  lon <- ncvar_get(nc_in, "longitude")      # [xi_rho x eta_rho x level x time]
  lat <- ncvar_get(nc_in, "latitude")       # [xi_rho x eta_rho x level x time]
  time <- ncvar_get(nc_in, "time")          # [xi_rho x eta_rho x level x time]
  
  
  # Identify variables to process (exclude coordinates, level, mask & dz layer)
  var_names <- setdiff(names(nc_in$var), c("longitude", "latitude", "level", "time", "mask", "dz"))
  
  for (var_name in var_names) {
    message("  Masking variable: ", var_name)
    
    # Load variable
    ## Note ncdf4 inverts the dimension order
    var_data <- ncvar_get(nc_in, var_name)  # [xi_rho x eta_rho x time]
    
    # Remove level dimension if only 1 length (everything has already been integrated)
    var_data <- drop(var_data) # Turns [xi_rho, eta_rho, 1, time] into [xi_rho, eta_rho, time]	
    
    # Apply mask across all time slices
    for (t in seq_len(dim(var_data)[3])) {
      var_data[,,t][mask == 0] <- NA
    }
    
    # Define dimensions
    dim_time <- ncdim_def("time", "seconds since 1900-01-01 00:00:00", vals = time, unlim = TRUE)
    dim_lat  <- ncdim_def("eta_rho", "degrees_north", vals = seq_len(dim(lat)[2]))
    dim_lon  <- ncdim_def("xi_rho", "degrees_east", vals = seq_len(dim(lon)[1]))
    
    # Define masked variable
    var_def <- ncvar_def(
      var_name,
      units = nc_in$var[[var_name]]$units,
      dim = list(dim_lon, dim_lat, dim_time), # [xi_rho x eta_rho x time]
      missval = NA,
      prec = "float"
    )
    
    # Define curvilinear coordinates and boundaries
    ## 1. Define the 'nbnds' dimension (for the 4 corners of each cell)
    dim_nbnds <- ncdim_def("nbnds", "", 1:4)
    
    ## 2. Define Lon/Lat and their associated Bounds
    lon_def <- ncvar_def("lon", "degrees_east", list(dim_lon, dim_lat), prec = "float")
    lat_def <- ncvar_def("lat", "degrees_north", list(dim_lon, dim_lat), prec = "float")
    
    lon_bnds_def <- ncvar_def("lon_bnds", "degrees_east", list(dim_nbnds, dim_lon, dim_lat), prec = "double")
    lat_bnds_def <- ncvar_def("lat_bnds", "degrees_north", list(dim_nbnds, dim_lon, dim_lat), prec = "double")

    ### Construct output filename
    #Note the published data will have different names: e.g., mld_hindcast_2008.nc
    masked_fname <- paste0(tools::file_path_sans_ext(basename(fl)),
                           "_", var_name, "_masked.nc")
    masked_full <- file.path(masked_subfolder, masked_fname)
    
    ## 3. Create output NetCDF
    nc_out <- nc_create(masked_full, list(var_def, lon_def, lat_def, lon_bnds_def, lat_bnds_def))
    
    ## 4. Generate the bounds data using dynamic midpoints
    ### Dimensions: [nbnds=4, xi=176, eta=536]
    lon_bnds <- array(NA, dim = c(4, dim(lon)[1], dim(lon)[2]))
    lat_bnds <- array(NA, dim = c(4, dim(lat)[1], dim(lat)[2]))

    # Calculate the grid spacing for each cell based on its neighbors
    ## For cells on the edge, we pad the diffs so every cell gets a width/height
    calc_delta <- function(m) {
      d <- matrix(0, nrow=nrow(m), ncol=ncol(m))
      ### Internal differences
      d[1:(nrow(m)-1), ] <- diff(m) 
      ### Fill the last row with the same diff as the penultimate row
      d[nrow(m), ] <- d[nrow(m)-1, ]
      return(d * 0.5) # We want half-distance to the next center
    }

    d_lon <- calc_delta(lon)
    d_lat <- calc_delta(t(lat)) # Transpose to get lat diffs along the other axis
    d_lat <- t(d_lat)

    # Define Corners based on local grid spacing
    ## This creates a "tiled" mesh that follows the curvilinear rotation
    ### Corner 1: Bottom Left
    lon_bnds[1,,] <- lon - d_lon; lat_bnds[1,,] <- lat - d_lat
    ### Corner 2: Bottom Right
    lon_bnds[2,,] <- lon + d_lon; lat_bnds[2,,] <- lat - d_lat
    ### Corner 3: Top Right
    lon_bnds[3,,] <- lon + d_lon; lat_bnds[3,,] <- lat + d_lat
    ### Corner 4: Top Left
    lon_bnds[4,,] <- lon - d_lon; lat_bnds[4,,] <- lat + d_lat

    # 5. Write data
    ncvar_put(nc_out, var_name, var_data)
    ncvar_put(nc_out, "lon", lon)
    ncvar_put(nc_out, "lat", lat)
    ncvar_put(nc_out, "lon_bnds", lon_bnds)
    ncvar_put(nc_out, "lat_bnds", lat_bnds)

    # 6. Linking Attributes (Enhanced for CDO/CF-Compliance)

    # --- Define Coordinate Metadata ---
    # CDO needs to know these are definitely horizontal axes
    ncatt_put(nc_out, "lon", "standard_name", "longitude")
    ncatt_put(nc_out, "lon", "units", "degrees_east")
    ncatt_put(nc_out, "lon", "axis", "X")
    ncatt_put(nc_out, "lon", "bounds", "lon_bnds")

    ncatt_put(nc_out, "lat", "standard_name", "latitude")
    ncatt_put(nc_out, "lat", "units", "degrees_north")
    ncatt_put(nc_out, "lat", "axis", "Y")
    ncatt_put(nc_out, "lat", "bounds", "lat_bnds")

   # --- Link Data Variable to Coordinates ---
   # Ensure the variable knows it is mapped to these 2D coordinates
   ncatt_put(nc_out, var_name, "coordinates", "lon lat")

   # --- Global Attributes ---
   # This is often the "magic" fix for 'Source grid cell corner coordinates missing'
   ncatt_put(nc_out, 0, "Conventions", "CF-1.6")
   ncatt_put(nc_out, 0, "history", paste("Created", Sys.time()))


   nc_close(nc_out)

  }
  
  nc_close(nc_in)

 # Final message outside the loop
 message("All files processed and masked!")

}
# --- End mask function ---

# --- Structure the bathymetry file ---
# Input ROMS files
bathy <- file.path(dir, "Data", "ROMS2", "bathy_nep.nc")

# Open bathymetry file
bathy_in <- nc_open(bathy)

# Load coordinates
## Note that in this file xi_rho and eta_rho are named differently!!!
lon <- ncvar_get(bathy_in, "lon_rho")      # [lon_rho x lat_rho]
lat <- ncvar_get(bathy_in, "lat_rho")      # [lon_rho x lat_rho]

# Identify variables to process (exclude coordinates, level, mask & dz layer)
var_name <- setdiff(names(bathy_in$var), c("lon_rho", "lat_rho"))

# Load variable
## Note ncdf4 inverts the dimension order
var_data <- ncvar_get(bathy_in, var_name)  # [lon_rho x lat_rho]

# Define dimensions
dim_lat  <- ncdim_def("lat_rho", "degrees_north", vals = seq_len(dim(lat)[2]))
dim_lon  <- ncdim_def("lon_rho", "degrees_east", vals = seq_len(dim(lon)[1]))

# Define masked variable
var_def <- ncvar_def(
      var_name,
      units = bathy_in$var[[var_name]]$units,
      dim = list(dim_lon, dim_lat), # [lon_rho x lat_rho]
      missval = NA,
      prec = "float"
    )

# Define curvilinear coordinates and boundaries
## 1. Define the 'nbnds' dimension (for the 4 corners of each cell)
dim_nbnds <- ncdim_def("nbnds", "", 1:4)

## 2. Define Lon/Lat and their associated Bounds
lon_def <- ncvar_def("lon", "degrees_east", list(dim_lon, dim_lat), prec = "float")
lat_def <- ncvar_def("lat", "degrees_north", list(dim_lon, dim_lat), prec = "float")

lon_bnds_def <- ncvar_def("lon_bnds", "degrees_east", list(dim_nbnds, dim_lon, dim_lat), prec = "double")
lat_bnds_def <- ncvar_def("lat_bnds", "degrees_north", list(dim_nbnds, dim_lon, dim_lat), prec = "double")

### Construct output filename
bathy_fname <- paste0(tools::file_path_sans_ext(basename(bathy)), "_", var_name, "_masked.nc")
bathy_full <- file.path(masked_subfolder, bathy_fname)

## 3. Create output NetCDF
bathy_out <- nc_create(bathy_full, list(var_def, lon_def, lat_def, lon_bnds_def, lat_bnds_def))

## 4. Generate the bounds data using dynamic midpoints
### Dimensions: [nbnds=4, xi=176, eta=536]
lon_bnds <- array(NA, dim = c(4, dim(lon)[1], dim(lon)[2]))
lat_bnds <- array(NA, dim = c(4, dim(lat)[1], dim(lat)[2]))

### Calculate the grid spacing for each cell based on its neighbors
#### For cells on the edge, we pad the diffs so every cell gets a width/height
calc_delta <- function(m) {
      d <- matrix(0, nrow=nrow(m), ncol=ncol(m))
      # Internal differences
      d[1:(nrow(m)-1), ] <- diff(m)
      # Fill the last row with the same diff as the penultimate row
      d[nrow(m), ] <- d[nrow(m)-1, ]
      return(d * 0.5) # We want half-distance to the next center
    }

d_lon <- calc_delta(lon)
d_lat <- calc_delta(t(lat)) # Transpose to get lat diffs along the other axis
d_lat <- t(d_lat)

### Define Corners based on local grid spacing
#### This creates a "tiled" mesh that follows the curvilinear rotation
##### Corner 1: Bottom Left
lon_bnds[1,,] <- lon - d_lon; lat_bnds[1,,] <- lat - d_lat
##### Corner 2: Bottom Right
lon_bnds[2,,] <- lon + d_lon; lat_bnds[2,,] <- lat - d_lat
##### Corner 3: Top Right
lon_bnds[3,,] <- lon + d_lon; lat_bnds[3,,] <- lat + d_lat
##### Corner 4: Top Left
lon_bnds[4,,] <- lon - d_lon; lat_bnds[4,,] <- lat + d_lat

## 5. Write data
ncvar_put(bathy_out, var_name, var_data)
ncvar_put(bathy_out, "lon", lon)
ncvar_put(bathy_out, "lat", lat)
ncvar_put(bathy_out, "lon_bnds", lon_bnds)
ncvar_put(bathy_out, "lat_bnds", lat_bnds)

## 6. Linking Attributes (Enhanced for CDO/CF-Compliance)

# --- Define Coordinate Metadata ---
# CDO needs to know these are definitely horizontal axes
ncatt_put(bathy_out, "lon", "standard_name", "longitude")
ncatt_put(bathy_out, "lon", "units", "degrees_east")
ncatt_put(bathy_out, "lon", "axis", "X")
ncatt_put(bathy_out, "lon", "bounds", "lon_bnds")

ncatt_put(bathy_out, "lat", "standard_name", "latitude")
ncatt_put(bathy_out, "lat", "units", "degrees_north")
ncatt_put(bathy_out, "lat", "axis", "Y")
ncatt_put(bathy_out, "lat", "bounds", "lat_bnds")

# --- Link Data Variable to Coordinates ---
# Ensure the variable knows it is mapped to these 2D coordinates
ncatt_put(bathy_out, var_name, "coordinates", "lon lat")

# --- Global Attributes ---
# This is often the "magic" fix for 'Source grid cell corner coordinates missing'
ncatt_put(bathy_out, 0, "Conventions", "CF-1.6")
ncatt_put(bathy_out, 0, "history", paste("Created", Sys.time()))

nc_close(bathy_out)

nc_close(bathy_in)

 # Final message outside the loop
 message("Bathymetry processed and structured!")

# --- End structure the bathymetry file ---

###-----------------------------------------------------------------
### If files are already masked and updated, pick up here!
# Note this is where the codes begins if using the published data
### ----------------------------------------------------------------

# --- Create empty grid ---
# 0. Setup Path and Projection ---
## Define the file path FIRST so it's ready for Step 5
target_grid_nc <- "target_10km_grid.nc"

# 1. Define Target Grid Space (10km)
## We use meters because the target is Albers
x_vals <- seq(from = floor(min_x/10000)*10000, to = ceiling(max_x/10000)*10000, by = 10000)
y_vals <- seq(from = floor(min_y/10000)*10000, to = ceiling(max_y/10000)*10000, by = 10000)

# 2. Calculate 2D Lat/Lon Matrices (The CDO Workaround)
## This translates every X/Y meter point into a geographic degree point
grid_df     <- expand.grid(x = x_vals, y = y_vals)
grid_sf     <- st_as_sf(grid_df, coords = c("x", "y"), crs = target_crs)
grid_wgs84  <- st_transform(grid_sf, 4326)
coords      <- st_coordinates(grid_wgs84)

## Reshape to 2D matrices for NetCDF
lon_2d <- matrix(coords[,1], nrow = length(x_vals), ncol = length(y_vals))
lat_2d <- matrix(coords[,2], nrow = length(x_vals), ncol = length(y_vals))

# 3. Calculate 2D Corners (For Conservative Remapping)
## Each cell has 4 corners: [bottom-left, bottom-right, top-right, top-left]
dx <- c(-5000, 5000, 5000, -5000)
dy <- c(-5000, -5000, 5000, 5000)

## Create 3D arrays to hold 4 corners for every (x,y) cell
lon_corners <- array(NA, dim = c(4, length(x_vals), length(y_vals)))
lat_corners <- array(NA, dim = c(4, length(x_vals), length(y_vals)))

## This will have (4 * length(x) * length(y)) rows
all_corners_df <- do.call(rbind, lapply(1:4, function(i) {
  expand.grid(x = x_vals + dx[i], y = y_vals + dy[i])
}))

## Transform everything in ONE single batch
all_corners_sf <- st_as_sf(all_corners_df, coords = c("x", "y"), crs = target_crs)
all_corners_wgs84 <- st_transform(all_corners_sf, 4326)
all_coords <- st_coordinates(all_corners_wgs84)

## Map the big matrix back into the 3D arrays
### Each slice represents one corner (1 to 4)
n_points <- length(x_vals) * length(y_vals)

for (i in 1:4) {
  start_idx <- ((i - 1) * n_points) + 1
  end_idx   <- i * n_points
  
  lon_corners[i,,] <- matrix(all_coords[start_idx:end_idx, 1], 
                             nrow = length(x_vals), ncol = length(y_vals))
  lat_corners[i,,] <- matrix(all_coords[start_idx:end_idx, 2], 
                             nrow = length(x_vals), ncol = length(y_vals))
}

# 4. Define NetCDF Dimensions and Variables
dim_x  <- ncdim_def("x", "meters", x_vals)
dim_y  <- ncdim_def("y", "meters", y_vals)
dim_nv <- ncdim_def("nv", "", 1:4, create_dimvar = FALSE) # 4 corners

## Variables
var_lon      <- ncvar_def("lon", "degrees_east", list(dim_x, dim_y), prec="double")
var_lat      <- ncvar_def("lat", "degrees_north", list(dim_x, dim_y), prec="double")
var_lon_bnds <- ncvar_def("lon_bnds", "degrees_east", list(dim_nv, dim_x, dim_y), prec="double")
var_lat_bnds <- ncvar_def("lat_bnds", "degrees_north", list(dim_nv, dim_x, dim_y), prec="double")
var_crs      <- ncvar_def("crs", "", list(), prec="char")
var_dummy    <- ncvar_def("dummy", "none", list(dim_x, dim_y), prec="float")

# 5. Create and Write
nc <- nc_create(target_grid_nc, list(var_lon, var_lat, var_lon_bnds, var_lat_bnds, var_crs, var_dummy))

ncvar_put(nc, var_lon, lon_2d)
ncvar_put(nc, var_lat, lat_2d)
ncvar_put(nc, var_lon_bnds, lon_corners)
ncvar_put(nc, var_lat_bnds, lat_corners)
ncvar_put(nc, var_dummy, matrix(0, nrow=length(x_vals), ncol=length(y_vals)))

# 6. CRITICAL Metadata for Conservative Remapping
## This is the "Magic Fix" that tells CDO where the corners are
ncatt_put(nc, "lon", "bounds", "lon_bnds")
ncatt_put(nc, "lat", "bounds", "lat_bnds")

## Re-link dummy to coordinates
ncatt_put(nc, "dummy", "coordinates", "lon lat")
ncatt_put(nc, "dummy", "grid_mapping", "crs")

## Albers Projection Details (for GIS/Compatibility)
ncatt_put(nc, "crs", "grid_mapping_name", "albers_conical_equal_area")
ncatt_put(nc, "crs", "standard_parallel", c(30, 50))
ncatt_put(nc, "crs", "longitude_of_central_meridian", -124.5)
ncatt_put(nc, "crs", "latitude_of_projection_origin", 40)
ncatt_put(nc, "crs", "false_easting", 0.0)
ncatt_put(nc, "crs", "false_northing", 0.0)
ncatt_put(nc, "crs", "crs_wkt", wkt_string)
ncatt_put(nc, "crs", "proj4text", proj_string)

## X and Y Axis Attributes (CRITICAL for CDO)
ncatt_put(nc, "x", "standard_name", "projection_x_coordinate")
ncatt_put(nc, "x", "axis", "X")
ncatt_put(nc, "y", "standard_name", "projection_y_coordinate")
ncatt_put(nc, "y", "axis", "Y")

## Global Attributes
ncatt_put(nc, 0, "crs_wkt", wkt_string)
ncatt_put(nc, 0, "Conventions", "CF-1.6")

nc_close(nc)

## Create the text grid description CDO uses for the actual math
system(paste("cdo griddes", target_grid_nc, "> target_10km_grid.txt"))
# --- End grid creation ---

# --- Regrid data ---
# Output folder for re-gridded files
## IMPORTANT: All paths here must be valid Linux paths, which /mnt/c/... is.
regrid_subfolder <- file.path(dir, "Data", "ROMS2", "Regridded", "10x10km")
if (!dir.exists(regrid_subfolder)) dir.create(regrid_subfolder, recursive = TRUE)

# No conversion needed: target_grid_nc is already the WSL/Linux path
## We rename it to be clearer:
linux_target_grid <- target_grid_nc

# Load masked file list (if using published data change file path accordingly)
fls.msk <- list.files(file.path(dir, "Data", "ROMS2", "Masked"), pattern = "\\.nc$", full.names = TRUE)

# Loop over all masked files
for (f in fls.msk){
  
  # 1. Extract base filename without extension
  var_file <- tools::file_path_sans_ext(basename(f))
  
  # 2. Skip grid files and dz files immediately
  if (grepl("target_grid|target_10km_grid|dz_masked", var_file, ignore.case = TRUE)) {
    next
  }
  
  # 3. Initialize variable for safety
  internal_var <- NA_character_
  
  # --- IF-ELSE EXTRACTION LOGIC ---

  # CASE A: Bathymetry files (published or masked)
    if (grepl("bathymetry|^h$|h_masked", var_file, ignore.case = TRUE)) {
        internal_var <- "h" # ROMS NetCDF internal variable name for bathymetry
        
  # CASE B: Non-published data (e.g., "allvars_mld_hindcast_2008_mld_masked")
    } else if (grepl("_masked$", var_file) || grepl("^allvars_", var_file)) {
      # Extract everything between 4-digit year and "_masked"
      internal_var <- str_extract(var_file, "(?<=\\d{4}_).*(?=_masked)")
      
      # Fallback if regex missed (e.g., missing year format)
      if (is.na(internal_var)) {
        cleaned <- sub("^allvars_", "", sub("_masked$", "", var_file))
        internal_var <- sub(".*_\\d{4}_", "", cleaned)
        }
      
      # CASE C: Published data (e.g., "mld_hindcast_2008.nc")
    } else {
      # Extract the first section before the first underscore ("mld")
      internal_var <- sub("_.*", "", var_file)
    }
  
  # --- SAFETY CHECK ---
  # Prevents CDO from receiving an 'NA' variable name and breaking the loop
  if (is.na(internal_var) || internal_var == "") {
    warning(paste("Could not extract internal_var for file:", var_file, "- Skipping."))
    next
  }
  
  # 4. Define output file (already the correct Linux path)
  linux_out_file <- file.path(regrid_subfolder, paste0(var_file, "_regrid.nc"))
  
  # 5. Choose interpolation method
  method <- if (grepl("chla_vint|microzooplankton_vint|mesozooplankton_vint|Pzooplankton_vint|u_vavg|v_vavg|w_vavg",
                      internal_var, ignore.case = TRUE)) {
    "remapcon"  # Conservative remapping for flux/integrated quantities
  } else {
    "remapbil"  # Bilinear remapping for standard 2D state variables
  }
  
  # 6. Build a "Clean" CDO command
  ## Note: No commas between operators in CDO, only between operator and its parameter.
  ## We use '-setgridtype,curvilinear' to tell CDO how to read the X/Y axes
  ## We use '-selname' to be absolutely sure CDO targets the right variable
  cmd <- sprintf(
      #### CDO is called directly
      "cdo -O %s,target_10km_grid.txt -setgrid,target_10km_grid.txt -selname,%s %s %s",
      method,            # replaces the first  %s
	    internal_var,      # replaces the second %s
      ### Input file (Linux)
      f,                 # replaces the third  %s
      ### Output file (Linux)
	    linux_out_file     # replaces the fourth %s
    )
  
  # 7. Execute Command
  message("Processing: ", var_file, " -> Variable: ", internal_var)
  system(cmd)
}

# Final message outside the loop
message("All files processed and regridded!")
# ---- End Regridding

# --- Create validation table ---
# Initialize a results list to store validation data
validation_results <- list()

for (f in fls.msk) {
    var_name <- tools::file_path_sans_ext(basename(f))
    linux_out_file <- file.path(regrid_subfolder, paste0(var_name, "_regrid.nc"))
    
    if (file.exists(linux_out_file)) {
      if (file.exists(linux_out_file)) {
        # 1. Get field-weighted mean of original file (redirect stderr 2>/dev/null to prevent R warnings)
        orig_raw <- system(paste0("cdo -s outputf,%g -fldmean ", f, " 2>/dev/null"), intern = TRUE)
        orig_mean_vec <- suppressWarnings(as.numeric(orig_raw))
        orig_mean_vec <- orig_mean_vec[!is.na(orig_mean_vec)] # Clean non-numeric lines
        
        # 2. Get field-weighted mean of regridded file
        regrid_raw <- system(paste0("cdo -s outputf,%g -fldmean ", linux_out_file, " 2>/dev/null"), intern = TRUE)
        regrid_mean_vec <- suppressWarnings(as.numeric(regrid_raw))
        regrid_mean_vec <- regrid_mean_vec[!is.na(regrid_mean_vec)]
        
        # 3. Calculate overall average means across all timesteps safely
        orig_mean <- mean(orig_mean_vec, na.rm = TRUE)
        regrid_mean <- mean(regrid_mean_vec, na.rm = TRUE)
        
        # 4. Calculate metrics safely
        pct_change <- if (!is.na(orig_mean) && orig_mean != 0) {
          ((regrid_mean - orig_mean) / orig_mean) * 100
        } else {
          NA_real_
        }
        
        # 5. Evaluate pass/check status
        status_check <- if (is.na(pct_change) || is.nan(pct_change) || is.infinite(pct_change)) {
          "NaN_ERROR"
        } else if (pct_change >= -10 && pct_change <= 10) {
          "PASS"
        } else {
          "CHECK"
        }
        
        # Store in List
        validation_results[[var_name]] <- data.frame(
          Variable       = var_name,
          Original_Mean  = round(orig_mean, 4),
          Regridded_Mean = round(regrid_mean, 4),
          Pct_Change     = round(pct_change, 3),
          Status         = status_check,
          stringsAsFactors = FALSE # Correct spelling
        )
      } else {
        warning("Regridded file missing for: ", var_name)
      }
    }
    
    # Combine into a final table safely
    if (length(validation_results) > 0) {
      final_table <- do.call(rbind, validation_results)
      
      # Print the table to the console
      print(final_table)
      
      # Save to CSV for documentation
      write.csv(final_table, file.path(regrid_subfolder, "regridding_validation_report.csv"), row.names = FALSE)
    } else {
      message("No regridded files were found to validate.")
    }

# --- End creation of validation table ---

# Combine into a final table
final_table <- do.call(rbind, validation_results)

# Print the table to the console
print(final_table)

# Save to CSV for your documentation
write.csv(final_table, file.path(regrid_subfolder, "regridding_validation_report.csv"), row.names = FALSE)
# --- End creation of validation table ---
