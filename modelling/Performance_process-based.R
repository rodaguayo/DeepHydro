# Hydrological modelling using TUWmodel  -----------------------------------------------------------

rm(list=ls())
cat("\014")  

library("hydroTSM")
library("hydroPSO")
library("terra")
library("arrow")
library("xts")

setwd("/home/rooda/OneDrive/Projects/DeepHydro/")
source("modelling/Process-based_models.R")
source("processing/metricsAux.R")

options(xts_check_TZ = FALSE)
path_pmet <- "/home/rooda/OneDrive/Projects/PatagoniaMet/"
path_disk <- "/home/rooda/Pipeline/DeepHydro/CLIMATE/catchments/"

# time period
period  <- seq(as.POSIXct("1997-01-01", tz = "UTC"), as.POSIXct("2019-12-31", tz = "UTC"), by = "day")
period_warm_up  <- seq(as.POSIXct("1997-01-01", tz = "UTC"), as.POSIXct("1999-12-31", tz = "UTC"), by = "day")
period_calib    <- seq(as.POSIXct("2000-01-01", tz = "UTC"), as.POSIXct("2019-12-31", tz = "UTC"), by = "day")

# Basin data: PMET-obs v1.1
basin_selection <- vect("data/GIS/Basins_PMETobs_points_subset.gpkg")
q_metadata      <- read.csv(paste0(path_pmet, "data/Zenodo/v11/Q_PMETobs_v11_metadata.csv"), row.names = "gauge_id")
q_metadata$gauge_code <- as.character(q_metadata$gauge_code)
q_metadata      <- q_metadata[rownames(q_metadata) %in% basin_selection$gauge_id, ] 

q_obs      <- read.csv(paste0(path_pmet, "data/Zenodo/v11/Q_PMETobs_1950_2020_v11d.csv"), row.names = "Date")
q_obs      <- q_obs[rownames(q_obs) %in% as.character(period_calib), ] 
q_obs      <- q_obs[,basin_selection$gauge_id]

# Climate and glacier runoff data

# PMETsim v1.1
pp   <- read_parquet(paste0(path_disk, "PP_ref_PMET_basins_eb.parquet"))
t2m  <- read_parquet(paste0(path_disk, "T2M_ref_PMET_basins_eb.parquet"))
pet  <- read_parquet(paste0(path_disk, "PET_ref_PMET_basins_eb.parquet"))
pp   <- pp[as.POSIXct(pp$date)   %in% period, ] 
t2m  <- t2m[as.POSIXct(t2m$date) %in% period, ] 
pet  <- pet[as.POSIXct(pet$date) %in% period, ] 

# simulations from OGGM
q_glacier <- read.csv("results/runoff/glacier_runoff_historical_pmet.csv", row.names = 1)
q_glacier <- q_glacier[rownames(q_glacier) %in% as.character(period_calib), ] 
q_glacier <- xts(q_glacier, order.by = as.Date(rownames(q_glacier)))

# Parameters: Lower and upper bound
lower_param_TUWmodel        <-  c(0.9,    0.0,   1.0,  -3.0,   -2.0,     0.0,     0,    0 ,  0.1,    2,   30,      1,       0,      0,        0)
upper_param_TUWmodel        <-  c(1.5,    5.0,   4.0,   1.01,   2.0,     1.0,   600,    20,    2,   30,  250,    100,       8,     30,       50)
names(upper_param_TUWmodel) <- param_names_TUWmodel
names(lower_param_TUWmodel) <- param_names_TUWmodel

lower_param_GR4J <- c(100, -10,  100,     1,  0,   0)
upper_param_GR4J <- c(3000,  3, 3000,    20,  1,   10)
names(upper_param_GR4J) <- param_names_GR4J
names(lower_param_GR4J) <- param_names_GR4J

# dataset to save performance and parameters
params_TUWmodel <- data.frame(matrix(ncol = 8 + length(lower_param_TUWmodel), nrow = 0)) 
df_sim_TUWmodel <- data.frame(matrix(ncol = 0, nrow = length(period_calib))) 

params_GR4J <- data.frame(matrix(ncol = 8 + length(lower_param_GR4J), nrow = 0)) 
df_sim_GR4J <- data.frame(matrix(ncol = 0, nrow = length(period_calib))) 

for (basin in rownames(q_metadata)) {
  
  meteo_i <- hydroclimate_basin(basin, obs = TRUE)
  
  # TUWmodel
  set.seed(123)
  out <- hydroPSO(fn="hydromodInR", lower=lower_param_TUWmodel, upper=upper_param_TUWmodel, 
                  method="spso2011", model.FUN="hydro_TUWmodel", 
                  control = list(write2disk=TRUE, MinMax="max", npart=32, maxit=100, REPORT=10, 
                                 reltol=1E-6, parallel = "parallel", par.nnodes = 32), 
                  model.FUN.args = list(obs = meteo_i$q_obs_i, 
                                        pp_i    = meteo_i$pp_i, 
                                        t2m_i   = meteo_i$t2m_i, 
                                        pet_i   = meteo_i$pet_i, 
                                        zones_i = meteo_i$zones_i,
                                        q_glacier_i = meteo_i$q_glacier_i,
                                        calibration = TRUE))
  
  TUWmodel_i <- hydro_TUWmodel(param.values = out$par, calibration = FALSE,
                               meteo_i$pp_i, meteo_i$t2m_i, meteo_i$pet_i,  
                               meteo_i$q_glacier_i, meteo_i$q_obs_i, meteo_i$zones_i)
  params_TUWmodel <- rbind(params_TUWmodel, c(basin, TUWmodel_i$metrics, out$par))
  df_sim_TUWmodel <- cbind(df_sim_TUWmodel, TUWmodel_i$q_sim)
  
  # GR4J
  set.seed(123)
  out <- hydroPSO(fn="hydromodInR", lower=lower_param_GR4J, upper=upper_param_GR4J, 
                  method="spso2011", model.FUN="hydro_GR4J", 
                  control = list(write2disk=TRUE, MinMax="max", npart=32, maxit=100, REPORT=10, 
                                 reltol=1E-6, parallel = "parallel", par.nnodes = 32), 
                  model.FUN.args = list(obs = meteo_i$q_obs_i, 
                                        pp_i    = meteo_i$pp_i, 
                                        t2m_i   = meteo_i$t2m_i, 
                                        pet_i   = meteo_i$pet_i, 
                                        zones_i = meteo_i$zones_i,
                                        q_glacier_i = meteo_i$q_glacier_i,
                                        calibration = TRUE))
  
  GR4J_i <- hydro_GR4J(param.values = out$par, calibration = FALSE,
                       meteo_i$pp_i,  meteo_i$t2m_i,  meteo_i$pet_i,  
                       meteo_i$q_glacier_i,  meteo_i$q_obs_i, meteo_i$zones_i)
  
  params_GR4J <- rbind(params_GR4J, c(basin, GR4J_i$metrics, out$par))
  df_sim_GR4J <- cbind(df_sim_GR4J, GR4J_i$q_sim)
  
  
  ## final plot
  plot(meteo_i$q_obs_i, type = "l")
  lines(TUWmodel_i$q_sim, type = "l", col = 2)
  lines(GR4J_i$q_sim, type = "l", col = 3)
  
  cat(basin, q_metadata[basin,]$gauge_name)
}

colnames(df_sim_TUWmodel) <- rownames(q_metadata)
rownames(df_sim_TUWmodel) <- period_calib
write.csv(df_sim_TUWmodel, "results/runoff/total_runoff_historical_TUWmodel_pmet.csv", row.names = TRUE)
colnames(df_sim_GR4J) <- rownames(q_metadata)
rownames(df_sim_GR4J) <- period_calib
write.csv(df_sim_GR4J, "results/runoff/total_runoff_historical_GR4J_pmet.csv", row.names = TRUE)

colnames(params_TUWmodel) <- c("gauge_id", metric_names, names(upper_param_TUWmodel))
write.csv(params_TUWmodel, "results/performance/Historical_all_TUWmodel.csv", row.names = FALSE)
colnames(params_GR4J) <- c("gauge_id", metric_names, names(upper_param_GR4J))
write.csv(params_GR4J, "results/performance/Historical_all_GR4J.csv", row.names = FALSE)
