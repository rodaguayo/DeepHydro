# Regionalization for process-based models  -------------------------------------------------------

rm(list=ls())
cat("\014")  

library("hydroTSM")
library("reshape2")
library("terra")
library("stars")
library("dplyr")
library("arrow")
library("xts")

setwd("/home/rooda/OneDrive/Projects/DeepHydro/")
source("modelling/Process-based_models.R")

path_pmet <- "/home/rooda/OneDrive/Projects/PatagoniaMet/"
path_disk <- "/home/rooda/Pipeline/DeepHydro/"

# 1. Data ---------------------------------------------------------------------

## Basin data: PMET-obs v1.1
basin_selection <- vect("data/GIS/Basins_PMETobs_points_subset.gpkg")
basin_centroids <- vect(paste0(path_pmet, "GIS/Basins_PMETobs_dev.shp"))
basin_centroids <- geom(centroids(basin_centroids), df = TRUE)

q_metadata_pmet      <- read.csv(paste0("data/Attributes_all_basins_pmet.csv"), row.names = "gauge_id")
q_metadata_pmet$cenlat <- basin_centroids$y
q_metadata_pmet$cenlon <- basin_centroids$x
q_metadata_pmet <- q_metadata_pmet[rownames(q_metadata_pmet) %in% basin_selection$gauge_id, ] 

q_metadata_all <- read.csv("data/Attributes_all_basins.csv", row.names =  "gauge_id")
q_metadata_all$glacier_dhdt[is.na(q_metadata_all$glacier_dhdt)] <- 0

## Parameters (all training)
params <- list(TUWmodel = read.csv("results/performance/Historical_all_TUWmodel.csv", row.names = "gauge_id"),
               GR4J     = read.csv("results/performance/Historical_all_GR4J.csv", row.names = "gauge_id"))

# 2. Regionalization ----------------------------------------------------------
variables = c('total_area', 'elev_median', 'slope_mean',
              'forest_cover', 'lake_cover', "herbaceous_veg_cover", 'glacier_cover',
              'p_mean_PMET', 'pet_mean_PMET', 'high_prec_freq_PMET')

# normalize between 0 and 1, check correlation using 
q_metadata <- rbind(q_metadata_all[,variables], q_metadata_pmet[,variables])
q_metadata_norm <- apply(q_metadata, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

# 3. Historical runoff (2000-2019) --------------------------------------------------------------

# 2000-2019
period          <- seq(as.POSIXct("1997-01-01", tz = "UTC"), as.POSIXct("2019-12-31", tz = "UTC"), by = "day")
period_warm_up  <- seq(as.POSIXct("1997-01-01", tz = "UTC"), as.POSIXct("1999-12-31", tz = "UTC"), by = "day")
period_calib    <- seq(as.POSIXct("2000-01-01", tz = "UTC"), as.POSIXct("2019-12-31", tz = "UTC"), by = "day")

## Climate: PMETsim v1.1

pp   <- read_parquet(paste0(path_disk, "CLIMATE/catchments/PP_ref_all_basins_eb.parquet"),  row.names = 1)  
t2m  <- read_parquet(paste0(path_disk, "CLIMATE/catchments/T2M_ref_all_basins_eb.parquet"), row.names = 1)
pet  <- read_parquet(paste0(path_disk, "CLIMATE/catchments/PET_ref_all_basins_eb.parquet"), row.names = 1)
pp   <- pp[as.POSIXct(pp$date)   %in% period, ] 
t2m  <- t2m[as.POSIXct(t2m$date) %in% period, ] 
pet  <- pet[as.POSIXct(pet$date) %in% period, ] 

## Glacier runoff: OGGM v1.5.4
q_glacier <- read.csv("results/runoff/glacier_runoff_historical_all.csv", row.names = 1)
q_glacier[ ,  setdiff( row.names(q_metadata_all), colnames(q_glacier))] <- 0
q_glacier <- q_glacier[rownames(q_glacier) %in% as.character(period_calib), ] 
q_glacier <- xts(q_glacier, order.by = as.Date(rownames(q_glacier)))

q_sim_TUWmodel <- data.frame(matrix(ncol = nrow(q_metadata_all), nrow = length(period_calib))) 
colnames(q_sim_TUWmodel) <- row.names(q_metadata_all)
q_sim_GR4J     <- data.frame(matrix(ncol = nrow(q_metadata_all), nrow = length(period_calib))) 
colnames(q_sim_GR4J) <- row.names(q_metadata_all)

for (basin in rownames(q_metadata_all)) { 
  
  basin_i  <- q_metadata_norm[basin,] 
  basins_i <- sweep(q_metadata_norm[rownames(q_metadata_pmet),], 2, basin_i)
  basins_i <- rowSums(abs(basins_i))
  best_basin <- names(which.min(basins_i))
  
  meteo_i <- hydroclimate_basin(basin, obs = FALSE)
  
  q_sim_i <- hydro_TUWmodel(param.values = params$TUWmodel[best_basin,][param_names_TUWmodel], 
                          meteo_i$pp_i, meteo_i$t2m_i,  meteo_i$pet_i, meteo_i$q_glacier_i,
                          calibration = FALSE, obs = NULL, zones_i = meteo_i$zones_i)
  q_sim_TUWmodel[,basin] <- q_sim_i
  
  q_sim_i <- hydro_GR4J(param.values = params$GR4J[best_basin,][param_names_GR4J], 
                      meteo_i$pp_i,  meteo_i$t2m_i,  meteo_i$pet_i,  
                      meteo_i$q_glacier_i,  meteo_i$zones_i,
                      calibration = FALSE, obs = NULL)
  q_sim_GR4J[,basin] <- q_sim_i
  print(basin)
}

q_sim_GR4J     <- cbind(date = period_calib, q_sim_GR4J)
q_sim_TUWmodel <- cbind(date = period_calib, q_sim_TUWmodel)

write.csv(q_sim_GR4J, "results/runoff/total_runoff_historical_GR4J_all.csv", row.names = FALSE)
write.csv(q_sim_TUWmodel, "results/runoff/total_runoff_historical_TUWmodel_all.csv", row.names = FALSE)

# 4. Future runoff (2021-2099) -----------------------------------------------------------------

period         <- seq(as.POSIXct("2021-01-01", tz = "UTC"), as.POSIXct("2098-12-31", tz = "UTC"), by = "day")
period_warm_up <- seq(as.POSIXct("2021-01-01", tz = "UTC"), as.POSIXct("2021-12-31", tz = "UTC"), by = "day")
period_calib   <- seq(as.POSIXct("2022-01-01", tz = "UTC"), as.POSIXct("2098-12-31", tz = "UTC"), by = "day")

gcms  <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", "MRI-ESM2-0")
gcms  <- c(paste(gcms, "ssp126", sep = "_"), paste(gcms, "ssp585", sep = "_"))

q_glacier_future <- read_mdim("results/runoff/glacier_runoff_future_all.nc", raster = c("gcm", "ssp"))

for (scenario in gcms) {
  
  # climate projections
  pp   <- read_parquet(paste0(path_disk, "CLIMATE/catchments/PP_",  scenario, "_all_basins_eb.parquet"), row.names = 1)
  pp   <- pp[as.POSIXct(pp$date)   %in% period, ]
  pp[pp < 0] <- 0 # 0.001 % of all values
  
  t2m  <- read_parquet(paste0(path_disk, "CLIMATE/catchments/T2M_", scenario, "_all_basins_eb.parquet"), row.names = 1)
  t2m  <- t2m[as.POSIXct(t2m$date)   %in% period, ] 
  
  pet  <- read_parquet(paste0(path_disk, "CLIMATE/catchments/PET_", scenario, "_all_basins_eb.parquet"), row.names = 1)
  pet  <- pet[as.POSIXct(pet$date)   %in% period, ] 

  q_glacier <- filter(q_glacier_future, ssp == strsplit(scenario, "_")[[1]][2], gcm == strsplit(scenario, "_")[[1]][1])
  q_glacier <- as.data.frame(q_glacier)[c("time", "rgi_id", "glacier_runoff")]
  
  q_glacier <- reshape2::dcast(q_glacier, time ~ rgi_id, value.var = "glacier_runoff")
  q_glacier[ ,  setdiff( row.names(q_metadata_all), colnames(q_glacier))] <- 0
  q_glacier <- xts(q_glacier, order.by = as.Date(q_glacier$time))
  q_glacier <- q_glacier[period_calib]
  
  # to fill
  q_sim_TUWmodel <- data.frame(matrix(ncol = nrow(q_metadata_all), nrow = length(period_calib))) 
  colnames(q_sim_TUWmodel) <- row.names(q_metadata_all)
  q_sim_GR4J     <- data.frame(matrix(ncol = nrow(q_metadata_all), nrow = length(period_calib))) 
  colnames(q_sim_GR4J) <- row.names(q_metadata_all)
  
  for (basin in rownames(q_metadata_all)) { 
    
    basin_i  <- q_metadata_norm[basin,] 
    basins_i <- sweep(q_metadata_norm[rownames(q_metadata_pmet),], 2, basin_i)
    basins_i <- rowSums(abs(basins_i))
    best_basin <- names(which.min(basins_i))
    
    meteo_i <- hydroclimate_basin(basin, obs = FALSE)
    
    q_sim_i <- hydro_TUWmodel(param.values = params$TUWmodel[best_basin,][param_names_TUWmodel], 
                              meteo_i$pp_i, meteo_i$t2m_i,  meteo_i$pet_i, meteo_i$q_glacier_i,
                              calibration = FALSE, obs = NULL, zones_i = meteo_i$zones_i)
    q_sim_TUWmodel[,basin] <- q_sim_i
    
    q_sim_i <- hydro_GR4J(param.values = params$GR4J[best_basin,][param_names_GR4J], 
                          meteo_i$pp_i,  meteo_i$t2m_i,  meteo_i$pet_i,  
                          meteo_i$q_glacier_i,  meteo_i$zones_i,
                          calibration = FALSE, obs = NULL)
    q_sim_GR4J[,basin] <- q_sim_i
    print(basin)
  }
  
  q_sim_TUWmodel <- cbind(date = period_calib, round(q_sim_TUWmodel, 3))
  q_sim_GR4J     <- cbind(date = period_calib, round(q_sim_GR4J, 3))
  
  write_parquet(q_sim_TUWmodel, paste("results/runoff/total_runoff_future", scenario, "TUWmodel_all.parquet", sep = "_"))
  write_parquet(q_sim_GR4J,     paste("results/runoff/total_runoff_future", scenario, "GR4J_all.parquet", sep = "_"))
  print(scenario)
  
}


