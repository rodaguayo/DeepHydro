# Historical and future climate for all catchments  -----------------------------------------------------------

rm(list=ls())
cat("\014")  

library("foreach")
library("doParallel")
library("exactextractr")

library("terra")
library("sf")
library("arrow")

# config terra cache for better performance
terra::gdalCache(30000) 

setwd("/home/rooda/OneDrive/Projects/DeepHydro/")
source("processing/elevationZones.R")

path_pmet <- "/home/rooda/OneDrive/Projects/PatagoniaMet/"
path_disk <- "/home/rooda/Pipeline/DeepHydro/CLIMATE/catchments/"

# function to extract time series for each catchment
ts_extract <- function(stack, shape) {
  ts <- round(t(exact_extract(stack, st_as_sf(shape), "mean", progress = T, max_cells_in_memory = 3e+08)), 2)
  
  # Check if there are any NA values
  if (any(is.na(ts))) {
    stop("Error: Extracted time series contains NA values.")
  }
  
  colnames(ts) <- shape$gauge_id
  rownames(ts) <- as.character(time(stack))
  ts <- as.data.frame(ts)
  ts <- cbind(date = row.names(ts), ts)
  ts$date <- as.Date(ts$date)
  return(ts)
}

# glaciers that are modeled
rgi7_hydro      <- vect("data/GIS/RGI7_Hydro.gpkg")

# dem 
dem <- rast(paste0(path_pmet, "GIS/dem_patagonia1.tif"))
dem <- subst(dem, NA, 0)
dem <- aggregate(dem, fact=3, fun="mean", cores = 20)
dem <- dem + init(dem, runif) # add random noise

# catchments: PMET dataset 
catchments_pmet <- vect(paste0(path_pmet, "GIS/Basins_PMETobs_dev.shp"))
basin_selection <- vect("data/GIS/Basins_PMETobs_points_subset.gpkg")
catchments_pmet <- catchments_pmet[catchments_pmet$gauge_id %in% basin_selection$gauge_id, ] 
catchments_pmet_ng <- erase(catchments_pmet, rgi7_hydro)

catchments_pmet_eb <- vect()
for(basin in 1:length(catchments_pmet_ng)) {
  
  # get the zones
  basin_ds <- elevationZones(x=catchments_pmet_ng[basin,], dem=dem, max.zones = 10, min.elevZ = 300)
  basin_i  <- as.polygons(basin_ds$zonesRaster, na.rm = TRUE, dissolve = TRUE)
  
  # assign id and zone
  basin_i$gauge_id <- paste0(catchments_pmet_ng[basin,]$gauge_id, "_", basin_ds$table$No.zone)
  catchments_pmet_eb <- rbind(catchments_pmet_eb, basin_i)
  print(basin)
}
  
# catchment: all 
catchments_all  <- vect("data/GIS/Basins_Patagonia_all.gpkg")
catchments_all_ng  <- erase(catchments_all, rgi7_hydro) 

## TODO this does not have an impact as the glacier runoff is weighted by the glacier area (??)
to_recover <- catchments_all$gauge_id[!(catchments_all$gauge_id %in% catchments_all_ng$gauge_id)]
to_recover <- catchments_all[catchments_all$gauge_id %in% to_recover]
catchments_all_ng <- rbind(catchments_all_ng, to_recover) # only 2

catchments_all_eb <- vect()
for(basin in 1:length(catchments_all_ng)) {
  
  # get the zones (if this fails, restart the session)
  basin_ds <- elevationZones(x=catchments_all_ng[basin,], dem=dem, max.zones = 10, min.elevZ = 300)
  basin_i  <- as.polygons(basin_ds$zonesRaster, na.rm = TRUE, dissolve = TRUE)
  
  # assign id and zone
  basin_i$gauge_id <- paste0(catchments_all_ng[basin,]$gauge_id, "_", basin_ds$table$No.zone)
  catchments_all_eb <- rbind(catchments_all_eb, basin_i)
  print(basin)
}

# Baseline climate: PMET -------------------------------------------------------------------------

pp   <- rast(paste0(path_pmet, "data/Zenodo/v11/PP_PMETsim_1980_2020_v11d.nc"))
pp   <- pp[[as.Date("1996-01-01") <= time(pp) &  time(pp) <= as.Date("2019-12-31")]]
pp   <- focal(pp, w=3, fun=mean, na.policy="only", na.rm=T) # cleaner fix for 1 basin
write_parquet(ts_extract(pp, catchments_pmet),    paste0(path_disk, "PP_ref_PMET_basins_full.parquet"))
write_parquet(ts_extract(pp, catchments_pmet_eb), paste0(path_disk, "PP_ref_PMET_basins_eb.parquet"))
write_parquet(ts_extract(pp, catchments_all),     paste0(path_disk, "PP_ref_all_basins_full.parquet"))
write_parquet(ts_extract(pp, catchments_all_eb),  paste0(path_disk, "PP_ref_all_basins_eb.parquet"))

t2m  <- rast(paste0(path_pmet, "data/Zenodo/v11/Tavg_PMETsim_1980_2020_v11d.nc"))
t2m   <- t2m[[as.Date("1996-01-01") <= time(t2m) &  time(t2m) <= as.Date("2019-12-31")]]
write_parquet(ts_extract(t2m, catchments_pmet),    paste0(path_disk, "T2M_ref_PMET_basins_full.parquet"))
write_parquet(ts_extract(t2m, catchments_pmet_eb), paste0(path_disk, "T2M_ref_PMET_basins_eb.parquet"))
write_parquet(ts_extract(t2m, catchments_all),     paste0(path_disk, "T2M_ref_all_basins_full.parquet"))
write_parquet(ts_extract(t2m, catchments_all_eb),  paste0(path_disk, "T2M_ref_all_basins_eb.parquet"))

tmax  <- rast(paste0(path_pmet, "data/Zenodo/v11/Tmax_PMETsim_1980_2020_v11d.nc"))
tmax  <- tmax[[as.Date("1996-01-01") <= time(tmax) &  time(tmax) <= as.Date("2019-12-31")]]
write_parquet(ts_extract(tmax, catchments_pmet),    paste0(path_disk, "TMAX_ref_PMET_basins_full.parquet"))
write_parquet(ts_extract(tmax, catchments_pmet_eb), paste0(path_disk, "TMAX_ref_PMET_basins_eb.parquet"))
write_parquet(ts_extract(tmax, catchments_all),     paste0(path_disk, "TMAX_ref_all_basins_full.parquet"))
write_parquet(ts_extract(tmax, catchments_all_eb),  paste0(path_disk, "TMAX_ref_all_basins_eb.parquet"))

tmin  <- rast(paste0(path_pmet, "data/Zenodo/v11/Tmin_PMETsim_1980_2020_v11d.nc"))
tmin  <- tmin[[as.Date("1996-01-01") <= time(tmin) &  time(tmin) <= as.Date("2019-12-31")]]
write_parquet(ts_extract(tmin, catchments_pmet),    paste0(path_disk, "TMIN_ref_PMET_basins_full.parquet"))
write_parquet(ts_extract(tmin, catchments_pmet_eb), paste0(path_disk, "TMIN_ref_PMET_basins_eb.parquet"))
write_parquet(ts_extract(tmin, catchments_all),     paste0(path_disk, "TMIN_ref_all_basins_full.parquet"))
write_parquet(ts_extract(tmin, catchments_all_eb),  paste0(path_disk, "TMIN_ref_all_basins_eb.parquet"))

pet  <- rast(paste0(path_pmet, "data/Evaporation/Ep_PMETsim_1980_2020d_dev.nc"))
pet  <- pet[[as.Date("1996-01-01") <= time(pet) &  time(pet) <= as.Date("2019-12-31")]]
write_parquet(ts_extract(pet, catchments_pmet),    paste0(path_disk, "PET_ref_PMET_basins_full.parquet"))
write_parquet(ts_extract(pet, catchments_pmet_eb), paste0(path_disk, "PET_ref_PMET_basins_eb.parquet"))
write_parquet(ts_extract(pet, catchments_all),     paste0(path_disk, "PET_ref_all_basins_full.parquet"))
write_parquet(ts_extract(pet, catchments_all_eb),  paste0(path_disk, "PET_ref_all_basins_eb.parquet"))

# Climate projections ----------------------------------------------------------------------------------

gcms  <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", "MRI-ESM2-0")
gcms  <- c(paste(gcms, "ssp126", sep = "_"), paste(gcms, "ssp585", sep = "_"))

variables <- c("PP", "T2M", "TASMAX", "TASMIN", "PET")

for (variable in variables) {
  for (gcm in gcms) {
    
    pattern <- paste(variable, gcm, sep = "_")
    stack   <- rast(list.files("~/Pipeline/DeepHydro/CLIMATE/future_bias_corrected", pattern, full.names = T))
    names(stack) <- time(stack)
    
    if (variable == "PP") {
      stack   <- focal(stack, w=3, fun=mean, na.policy="only", na.rm=T) 
      }
    
    write_parquet(ts_extract(stack, catchments_all),    paste(paste0(path_disk, variable), gcm,  "all_basins_full.parquet", sep = "_"))
    write_parquet(ts_extract(stack, catchments_all_eb), paste(paste0(path_disk, variable), gcm,  "all_basins_eb.parquet",   sep = "_"))
    cat(gcm, variable)
  }
}
