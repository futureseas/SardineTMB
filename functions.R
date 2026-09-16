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

#-------------------------------------------------------------------------------
# Model Comparison -- Not used any more!
#-------------------------------------------------------------------------------

# Function to extract criteria
model_summary <- function(fit, data) {
  preds <- predict(fit, type = "response")
  response <- data[[as.character(fit$call$formula[[2]])]]
  
  # Calculate Standard R²
  ss_res <- sum((response - preds$est)^2, na.rm = TRUE)
  ss_tot <- sum((response - mean(response, na.rm = TRUE))^2, na.rm = TRUE)
  r2Standard <- 1 - (ss_res / ss_tot)
  # Obtain R² of deired variable when comparing the fitted model's residual sum 
  # of squares to the total sum of squares from a null model that is forced through
  # the origin (0,0). Note this may skew the model results closer to 1, which really
  # assessing model fit.
  r2 <- summary(lm(Dnsty_kg_m2 ~ est, data = preds))$r.squared
  # Get mean absolute error of desired variable
  MAE = mean(abs(preds$Dnsty_kg_m2 - preds$est))
  
  c(
    ML = as.numeric(logLik(fit)),
    df = attr(logLik(fit), "df"),
    R2.SS_tot_mean = r2Standard,
    R2.SS_tot_0 = r2,
    MAE = MAE,
    AIC = AIC(fit)
  )
}

#--------------------------------------------------------
# Fit for each haul, a power function to derive length 
# weight relationships
#--------------------------------------------------------

# --- STEP 1: Global Model (Run once) ---
# This gives you the stable "Species Standard"
#global_model_data <- sardine.obs %>% 
#  filter(!is.na(weight), !is.na(standard_length))

#m_global <- lm(log(weight) ~ log(standard_length), data = global_model_data)
#a_global <- exp(coef(m_global)[1])
#b_global <- coef(m_global)[2]

# --- STEP 2: The Loop Function ---
#process_haul_stable <- function(haul_data, keys) {
#  
#  # Predict weight for every fish using the stable Global Model
#  # We use the 'hybrid' approach: use real weight if you have it, else predict
#  haul_data <- haul_data %>%
#    mutate(hybrid_weight = ifelse(!is.na(weight), 
#                                  weight, 
#                                  a_global * (standard_length ^ b_global)))
#  
#  # Calculate scaling factor to hit the reported CATCH_WEIGHT
#  actual_total <- unique(haul_data$CATCH_WEIGHT)
#  sum_hybrid <- sum(haul_data$hybrid_weight, na.rm = TRUE)
#  
#  if(!is.na(sum_hybrid) && sum_hybrid > 0) {
#    scaling_factor <- actual_total / sum_hybrid
#    
#    # Final Weight: The individual predicted weight adjusted for the total catch
#    haul_data <- haul_data %>%
#      mutate(final_weight = hybrid_weight * scaling_factor)
#  } else {
#    haul_data <- haul_data %>% mutate(final_weight = as.numeric(NA))
#  }
#  
#  return(haul_data)
#}
#
# --- STEP 3: Execute ---
#final_dataset <- partial.sardine.obs %>%
#  group_by(keys) %>%
#  group_modify(~ process_haul_stable(.x, .y)) %>% 
#  ungroup()
#
#final_dataset2 <- final_dataset %>%
#  rename(original.weight = weight) %>%
#  # rename predicted weights as original
#  rename(weight = hybrid_weight)
