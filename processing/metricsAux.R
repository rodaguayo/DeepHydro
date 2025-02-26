# Calculate FMS, FHV and FLV ----------------------------------------------------------------------

library(dplyr)

validate_inputs <- function(obs, sim) {
  if (length(obs) != length(sim)) {
    stop("Observed and simulated time series must have the same length")
  }
}

mask_valid <- function(obs, sim) {
  valid_indices <- !is.na(obs) & !is.na(sim)
  list(obs = obs[valid_indices], sim = sim[valid_indices])
}

get_fdc <- function(series) {
  sort(series, decreasing = TRUE)
}


# FMS -----
fdc_fms <- function(obs, sim, lower = 0.2, upper = 0.7) {

  validate_inputs(obs, sim)
  
  masked_data <- mask_valid(obs, sim)
  obs <- masked_data$obs
  sim <- masked_data$sim
  
  if (length(obs) < 1) {
    return(NA)
  }
  
  if (any(lower <= 0, upper <= 0, lower >= 1, upper >= 1)) {
    stop("upper and lower have to be in range ]0,1[")
  }
  
  if (lower >= upper) {
    stop("The lower threshold has to be smaller than the upper.")
  }
  
  obs <- get_fdc(obs)
  sim <- get_fdc(sim)
  
  sim[sim <= 0] <- 1e-6
  obs[obs == 0] <- 1e-6
  
  qsm_lower <- log(sim[round(lower * length(sim))])
  qsm_upper <- log(sim[round(upper * length(sim))])
  qom_lower <- log(obs[round(lower * length(obs))])
  qom_upper <- log(obs[round(upper * length(obs))])
  
  fms <- ((qsm_lower - qsm_upper) - (qom_lower - qom_upper)) / (qom_lower - qom_upper + 1e-6)
  
  return(fms * 100)
}

# FHV -----
fdc_fhv <- function(obs, sim, h = 0.02) {
  validate_inputs(obs, sim)
  
  masked_data <- mask_valid(obs, sim)
  obs <- masked_data$obs
  sim <- masked_data$sim
  
  if (length(obs) < 1) {
    return(NA)
  }
  
  if (h <= 0 || h >= 1) {
    stop("h has to be in range ]0,1[. Consider small values, e.g. 0.02 for 2% peak flows")
  }
  
  obs <- get_fdc(obs)
  sim <- get_fdc(sim)
  
  obs <- obs[1:round(h * length(obs))]
  sim <- sim[1:round(h * length(sim))]
  
  fhv <- sum(sim - obs) / sum(obs)
  
  return(fhv * 100)
}

# FLV -----
fdc_flv <- function(obs, sim, l = 0.3) {
  validate_inputs(obs, sim)
  
  masked_data <- mask_valid(obs, sim)
  obs <- masked_data$obs
  sim <- masked_data$sim
  
  if (length(obs) < 1) {
    return(NA)
  }
  
  if (l <= 0 || l >= 1) {
    stop("l has to be in range ]0,1[. Consider values like 0.3 for 30% low flows")
  }
  
  obs <- get_fdc(obs)
  sim <- get_fdc(sim)
  
  sim[sim <= 0] <- 1e-6
  obs[obs == 0] <- 1e-6
  
  obs <- obs[(length(obs) - round(l * length(obs)) + 1):length(obs)]
  sim <- sim[(length(sim) - round(l * length(sim)) + 1):length(sim)]
  
  obs <- log(obs)
  sim <- log(sim)
  
  qsl <- sum(sim - min(sim))
  qol <- sum(obs - min(obs))
  
  flv <- -1 * (qsl - qol) / (qol + 1e-6)
  
  return(flv * 100)
  }