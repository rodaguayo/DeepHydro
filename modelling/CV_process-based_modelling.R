# Leave-one-out Cross-Validation  ----------------------------------------------------------------

rm(list=ls())
cat("\014")  

library("hydroTSM")
library("terra")
library("xts")
library("arrow")

setwd("/home/rooda/OneDrive/Projects/DeepHydro/")
source("modelling/Process-based_models.R")
source("processing/metricsAux.R")

path_pmet <- "/home/rooda/OneDrive/Projects/PatagoniaMet/"
path_disk <- "/home/rooda/Pipeline/DeepHydro/"

period  <- seq(as.POSIXct("1997-01-01", tz = "UTC"), as.POSIXct("2019-12-31", tz = "UTC"), by = "day")
period_warm_up  <- seq(as.POSIXct("1997-01-01", tz = "UTC"), as.POSIXct("1999-12-31", tz = "UTC"), by = "day")
period_calib    <- seq(as.POSIXct("2000-01-01", tz = "UTC"), as.POSIXct("2019-12-31", tz = "UTC"), by = "day")

# 1. Data ----------------------------------------------------------------------------------------

## Basin data: PMET-obs v1.1
basin_selection <- vect("data/GIS/Basins_PMETobs_points_subset.gpkg")
q_shape         <- vect(paste0(path_pmet, "GIS/Basins_PMETobs_dev.shp"))
q_shape         <- q_shape[q_shape$gauge_id %in% basin_selection$gauge_id, ] 
basin_centroids <- geom(centroids(q_shape), df = TRUE)

q_metadata       <- read.csv(paste0("data/Attributes_all_basins_pmet.csv"), row.names = "gauge_id")
q_metadata       <- q_metadata[rownames(q_metadata) %in% basin_selection$gauge_id, ] 

q_metadata$kfold_pur  <- basin_selection$kfold_pur
q_metadata$kfold_pub  <- basin_selection$kfold_pub

q_obs      <- read.csv(paste0(path_pmet, "data/Zenodo/v11/Q_PMETobs_1950_2020_v11d.csv"), row.names = "Date")
q_obs      <- q_obs[rownames(q_obs) %in% as.character(period_calib), ] 
q_obs      <- q_obs[,basin_selection$gauge_id]

## Climate: PMETsim v1.1
pp   <- read_parquet(paste0(path_disk, "CLIMATE/catchments/PP_ref_PMET_basins_eb.parquet"))  
t2m  <- read_parquet(paste0(path_disk, "CLIMATE/catchments/T2M_ref_PMET_basins_eb.parquet"))
pet  <- read_parquet(paste0(path_disk, "CLIMATE/catchments/PET_ref_PMET_basins_eb.parquet"))
pp   <- pp[as.POSIXct(pp$date)   %in% period, ] 
t2m  <- t2m[as.POSIXct(t2m$date) %in% period, ] 
pet  <- pet[as.POSIXct(pet$date) %in% period, ] 

## Glacier runoff: OGGM v1.5.4
q_glacier <- read.csv("results/runoff/glacier_runoff_historical_pmet.csv", row.names = 1)
q_glacier <- q_glacier[rownames(q_glacier) %in% as.character(period_calib), ] 
q_glacier <- xts(q_glacier, order.by = as.Date(rownames(q_glacier)))

## Parameters (all training)
params <- list(TUWmodel = read.csv("results/performance/Historical_all_TUWmodel.csv", row.names = "gauge_id"),
               GR4J     = read.csv("results/performance/Historical_all_GR4J.csv", row.names = "gauge_id"))

# 2. Cross-validation  --------------------------------------------------------------------------
variables = c('total_area', 'elev_median', 'slope_mean',
              'forest_cover', 'lake_cover', "herbaceous_veg_cover", 'glacier_cover',
              'p_mean_PMET', 'pet_mean_PMET', 'high_prec_freq_PMET')

# normalize 0-1
q_metadata_norm <- q_metadata[,variables]
q_metadata_norm <- apply(q_metadata_norm, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
#lares::corr_cross(as.data.frame(q_metadata_norm)) #  check correlation using 

# 2.1 PUR  --------------------------------------------------------------------------------------
metric_dataset <- data.frame(matrix(ncol = 7, nrow = 0)) 
df_sim_GR4J <- data.frame(matrix(ncol = 0, nrow = length(period_calib))) 
df_sim_TUWmodel <- data.frame(matrix(ncol = 0, nrow = length(period_calib))) 

for (basin in rownames(q_metadata)) { 

  basin_i  <- q_metadata_norm[basin,]
  
  # subset k fold and search for the minimum
  q_metadata_norm_sub <- q_metadata_norm[!(q_metadata$kfold_pur %in% q_metadata[basin,]$kfold_pur), ]
  basins_i <- sweep(q_metadata_norm_sub, 2, basin_i)
  basins_i <- rowSums(abs(basins_i))
  best_basin <- names(which.min(basins_i))

  meteo_i <- hydroclimate_basin(basin, obs = TRUE)
  
  TUWmodel_i <- hydro_TUWmodel(param.values = params$TUWmodel[best_basin,][param_names_TUWmodel], 
                               calibration = FALSE,
                               meteo_i$pp_i, meteo_i$t2m_i, meteo_i$pet_i,  
                               meteo_i$q_glacier_i, meteo_i$q_obs_i, meteo_i$zones_i)
  
  metric_dataset <- rbind(metric_dataset, c("TUWmodel", basin, TUWmodel_i$metrics))
  df_sim_TUWmodel <- cbind(df_sim_TUWmodel, TUWmodel_i$q_sim)
  
  GR4J_i <- hydro_GR4J(param.values = params$GR4J[best_basin,][param_names_GR4J], 
                      meteo_i$pp_i,  meteo_i$t2m_i,  meteo_i$pet_i,  
                      meteo_i$q_glacier_i,  meteo_i$q_obs_i, meteo_i$zones_i, 
                      calibration = FALSE)
  metric_dataset <- rbind(metric_dataset, c("GR4J", basin, GR4J_i$metrics))
  df_sim_GR4J <- cbind(df_sim_GR4J, GR4J_i$q_sim)

  print(basin)
}

colnames(df_sim_GR4J)     <- rownames(q_metadata)
colnames(df_sim_TUWmodel) <- rownames(q_metadata)
rownames(df_sim_GR4J)     <- period_calib
rownames(df_sim_TUWmodel) <- period_calib
write.csv(df_sim_GR4J,     "results/runoff/total_runoff_historical_GR4J_pmet_CV_PUR.csv", row.names = TRUE)
write.csv(df_sim_TUWmodel, "results/runoff/total_runoff_historical_TUWmodel_pmet_CV_PUR.csv", row.names = TRUE)

colnames(metric_dataset) <- c("Model", "Basin", metric_names)
write.csv(metric_dataset, "results/performance/Historical_CV_PUR_process_based.csv", row.names = FALSE)

# 2.2 PUB  --------------------------------------------------------------------------------------
metric_dataset <- data.frame(matrix(ncol = 7, nrow = 0)) 
df_sim_GR4J <- data.frame(matrix(ncol = 0, nrow = length(period_calib))) 
df_sim_TUWmodel <- data.frame(matrix(ncol = 0, nrow = length(period_calib))) 

for (basin in rownames(q_metadata)) { 
  
  basin_i  <- q_metadata_norm[basin,]
  
  # subset k fold and search for the minimum
  q_metadata_norm_sub <- q_metadata_norm[!(q_metadata$kfold_pub %in% q_metadata[basin,]$kfold_pub), ]
  basins_i <- sweep(q_metadata_norm_sub, 2, basin_i)
  basins_i <- rowSums(abs(basins_i))
  best_basin <- names(which.min(basins_i))
  
  meteo_i <- hydroclimate_basin(basin, obs = TRUE)
  
  TUWmodel_i <- hydro_TUWmodel(param.values = params$TUWmodel[best_basin,][param_names_TUWmodel], 
                               calibration = FALSE,
                               meteo_i$pp_i, meteo_i$t2m_i, meteo_i$pet_i,  
                               meteo_i$q_glacier_i, meteo_i$q_obs_i, meteo_i$zones_i)
  
  metric_dataset <- rbind(metric_dataset, c("TUWmodel", basin, TUWmodel_i$metrics))
  df_sim_TUWmodel <- cbind(df_sim_TUWmodel, TUWmodel_i$q_sim)
  
  GR4J_i <- hydro_GR4J(param.values = params$GR4J[best_basin,][param_names_GR4J], 
                       meteo_i$pp_i,  meteo_i$t2m_i,  meteo_i$pet_i,  
                       meteo_i$q_glacier_i,  meteo_i$q_obs_i, meteo_i$zones_i, 
                       calibration = FALSE)
  metric_dataset <- rbind(metric_dataset, c("GR4J", basin, GR4J_i$metrics))
  df_sim_GR4J <- cbind(df_sim_GR4J, GR4J_i$q_sim)
  
  print(basin)
}

colnames(df_sim_GR4J)     <- rownames(q_metadata)
colnames(df_sim_TUWmodel) <- rownames(q_metadata)
rownames(df_sim_GR4J)     <- period_calib
rownames(df_sim_TUWmodel) <- period_calib
write.csv(df_sim_GR4J,     "results/runoff/total_runoff_historical_GR4J_pmet_CV_PUB.csv", row.names = TRUE)
write.csv(df_sim_TUWmodel, "results/runoff/total_runoff_historical_TUWmodel_pmet_CV_PUB.csv", row.names = TRUE)

colnames(metric_dataset) <- c("Model", "Basin", metric_names)
write.csv(metric_dataset, "results/performance/Historical_CV_PUB_process_based.csv", row.names = FALSE)

