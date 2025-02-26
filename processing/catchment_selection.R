# Catchment selection  --------------------------------------------------------------

rm(list=ls())
cat("\014")  

library("caret")
library("terra")
library("sf")

options(xts_check_TZ = FALSE)
setwd("/home/rooda/OneDrive/Projects/DeepHydro/")
path_pmet <- "/home/rooda/OneDrive/Projects/PatagoniaMet/"

# threshold
threshold <- 0.75

# time period
period_calib    <- seq(as.POSIXct("2000-01-01", tz = "UTC"), 
                       as.POSIXct("2019-12-31", tz = "UTC"), by = "day")

# Basin data
q_metadata <- read.csv(paste0(path_pmet, "data/Zenodo/v11/Q_PMETobs_v11_metadata.csv"), row.names = "gauge_id")
q_metadata$gauge_code <- as.character(q_metadata$gauge_code)
q_obs      <- read.csv(paste0(path_pmet, "data/Zenodo/v11/Q_PMETobs_1950_2020_v11d.csv"), row.names = "Date")
q_obs      <- q_obs[rownames(q_obs) %in% as.character(period_calib), ] 
q_shape    <- vect(paste0(path_pmet, "GIS/Basins_PMETobs_points_dev.shp"))

# subset based on data available
basin_select <- colSums(!is.na(q_obs))/nrow(q_obs)
q_metadata   <- q_metadata[(basin_select > threshold),] 
q_shape   <- q_shape[(basin_select > threshold),] 

# subset based on dam presence
basin_select <- q_metadata$dam < 1
q_metadata <- q_metadata[basin_select,]
q_shape   <- q_shape[basin_select,] 


# split into 10 groups
q_shape$kfold_pur <- as.integer(cut(1:length(q_shape$gauge_id), breaks = 10, labels = FALSE))

# split into 10 groups (randomly)
set.seed(123)
q_shape <- q_shape[sample(nrow(q_shape)), ]
q_shape$kfold_pub <- as.integer(cut(1:length(q_shape$gauge_id), breaks = 10, labels = FALSE))
q_shape <- q_shape[order(q_shape$gauge_lat, decreasing = T), ]

# save
writeVector(q_shape, "data/GIS/Basins_PMETobs_points_subset.gpkg", overwrite = T, options=NULL)
