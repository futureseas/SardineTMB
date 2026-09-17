#-------------------------------------------------------------------------------
# Make a transect spatial object
#-------------------------------------------------------------------------------

make_transect_sf <- function(df,
                             lat_start_col = "start_latitude", lon_start_col = "start_longitude",
                             lat_stop_col  = "stop_latitude",  lon_stop_col  = "stop_longitude",
                             id_col = "FID", crs = crs) {
  library(data.table)
  library(sf)
  
  df <- as.data.table(df)
  
  # Create list of LINESTRING coordinate matrices
  lines_list <- lapply(seq_len(nrow(df)), function(i) {
    start_lon <- df[[lon_start_col]][i]
    start_lat <- df[[lat_start_col]][i]
    stop_lon  <- df[[lon_stop_col]][i]
    stop_lat  <- df[[lat_stop_col]][i]
    
    # Check for NA values first. If any coordinate is NA, we return an empty geometry.
    # st_linestring() will handle the empty matrix gracefully.
    if (is.na(start_lon) || is.na(start_lat) || is.na(stop_lon) || is.na(stop_lat)) {
      return(st_linestring()) # Returns an empty matrix
    }
    
    # Check if start and end points are identical
    if (start_lon == stop_lon && start_lat == stop_lat) {
      # Return a matrix with a single coordinate pair
      return(matrix(c(start_lon, start_lat, start_lon + 1e-6, stop_lat + 1e-6 # adds small jitter to stop point
                      ), ncol = 2, byrow = TRUE))
    } else {
      # Return a matrix with the two coordinate pairs for a standard line
      return(matrix(c(
        start_lon, start_lat,
        stop_lon,  stop_lat
      ), ncol = 2, byrow = TRUE))
    }
  })
  
  # Convert to sf LINESTRING geometry
  # st_linestring() will take a list of matrices and create a list of `sfg` LINESTRING objects
  line_geoms <- st_sfc(lapply(lines_list, st_linestring), crs = crs)
  
  # Return sf object
  st_sf(data.frame(df), geometry = line_geoms)
}

#-------------------------------------------------------------------------------
# Define function to get sunrise/sunset
#-------------------------------------------------------------------------------

get_sun_times <- function(lat, lon, date) {
  times <- getSunlightTimes(date = date, 
                            lat = lat, 
                            lon = lon, 
                            keep = c("sunrise", "sunset"), 
                            # To update with desired tz
                            tz = "Etc/GMT+8")
  tibble(
    date = date,
    sunrise_pacific = times$sunrise,
    sunset_pacific = times$sunset,
    lat = lat,
    lon = lon
  )
}

#-------------------------------------------------------------------------------
# Spatial difference matrix -- Not used any more!
#-------------------------------------------------------------------------------


compute_spatial_dist_matrix <- function(tar_coords, knwn_coords) {
  # Ensure coordinates are matrices
  tar_coords <- as.matrix(tar_coords)
  knwn_coords <- as.matrix(knwn_coords)
  
  # Combine for pairwise distance computation
  combined_coords <- rbind(tar_coords, knwn_coords)
  
  # Compute full distance matrix
  full_dist <- as.matrix(dist(combined_coords))
  
  # Extract distances from tar_coords to knwn_coords
  n_tar <- nrow(tar_coords)
  n_knwn <- nrow(knwn_coords)
  
  full_dist[1:n_tar, (n_tar + 1):(n_tar + n_knwn)]
}

#-------------------------------------------------------------------------------
# Temporal difference matrix
#-------------------------------------------------------------------------------

## Note that time_vec(s) are POSIXct (or Date) timestamps
## units can be "secs", "mins", "hours", "days", or "weeks" — same as in difftime.
compute_time_diff_matrix <- function(time_vec1, time_vec2, units = "secs"){
  outer(time_vec1, time_vec2, 
        FUN = function(t1, t2) abs(as.numeric(difftime(t1, t2, units = units))))
}


