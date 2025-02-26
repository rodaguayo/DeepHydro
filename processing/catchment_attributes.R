# Calculate catchment-wide attribute -------------------------------------------------------------

rm(list=ls())
cat("\014")  

library("exactextractr")
library("zoo")
library("sf")
library("terra")

terra::gdalCache(40000)

setwd("/home/rooda/OneDrive/Projects/DeepHydro/")
path_pmet <- "/home/rooda/OneDrive/Projects/PatagoniaMet/"

# Basins to simulate
basin_shp  <- st_read("data/GIS/Basins_Patagonia_all.gpkg")
basin_shp  <- st_read(paste0(path_pmet, "GIS/Basins_PMETobs_dev.shp"))

# Reference period
period <- c(as.POSIXct("1990-01-01"), as.POSIXct("2019/12/31"))

# a) Topographic attributes ---------------------------------------------------------------------

basin_shp$total_area  <- expanse(vect(basin_shp), unit="km")
basin_shp$total_perim <- perim(vect(basin_shp))
basin_shp$gc_index <- basin_shp$total_perim/(2*(sqrt(pi*basin_shp$total_area)))

dem   <- rast(paste0(path_pmet, "GIS/dem_patagonia3f.tif"))
dem   <- subst(dem, NA, 0)  # NAs to sea level (= 0)
basin_shp$elev_mean   <- round(exact_extract(dem,   basin_shp, "mean"),   1)
basin_shp$elev_median <- round(exact_extract(dem,   basin_shp, "median"), 1)
basin_shp$elev_max    <- round(exact_extract(dem,   basin_shp, "max"),    1)
basin_shp$elev_sd     <- round(exact_extract(dem,   basin_shp, "stdev"),  1)

slope <- terrain(dem, v='aspect', unit='degrees')
basin_shp$slope_mean  <- round(exact_extract(slope, basin_shp, "median"), 1)

# b) Land cover ----------------------------------------------------------------------------------

# leaf area index (LAI)
lai_data   <- rast(paste0(path_pmet, "GIS/LAI_MOD15A2H_climatology.nc"))
lai_data   <- round(exact_extract(lai_data, basin_shp, "mean"), 3)
basin_shp$lai_max  <- apply(lai_data, 1, max, na.rm=TRUE)
basin_shp$lai_diff <- apply(lai_data, 1, max, na.rm=TRUE) - apply(lai_data, 1, min, na.rm=TRUE)

land_cover_data   <- rast(paste0(path_pmet, "GIS/Land_cover/Land_cover_copernicus_100m_2019.tif"))
land_cover_data   <- classify(land_cover_data, matrix(c(100, 199, 10), ncol = 3, byrow = TRUE))
land_cover_data   <- extract(land_cover_data, basin_shp, fun="table", na.rm=TRUE, exact=FALSE)
land_cover_data$ID <- NULL
land_cover_data <- land_cover_data*100/rowSums(land_cover_data)

basin_shp$forest_cover <- round(land_cover_data[["10"]], 3) # All forest classes
basin_shp$shrubs_veg_cover <- round(land_cover_data[["20"]], 3)
basin_shp$herbaceous_veg_cover <- round(land_cover_data[["30"]], 3) 
basin_shp$sparse_veg_cover <- round(land_cover_data[["60"]], 3)

lake_cover <- rast(paste0(path_pmet, "GIS/lc_water_500m.tif"))
lake_cover <- subst(lake_cover, NA, 0)  # NAs to 0 %
basin_shp$lake_cover   <- round(exact_extract(lake_cover, basin_shp, "mean"), 1)

# glacier features (RGI7)
glaciers  <- vect(paste0(path_pmet, "GIS/Glaciers/RGI7.shp"))
glaciers  <- subset(glaciers, glaciers$cenlat <= -40.5)
glaciers  <- rasterize(glaciers, dem, background = 0) * 100
basin_shp$glacier_cover <- round(exact_extract(glaciers,   basin_shp, "mean"), 1)

# Glacier change (elevation in m y-1)
dh_dt <- rast(paste0(path_pmet, "GIS/Glaciers/dhdt_2021.tif"))
dh_dt <- project(dh_dt, crs(basin_shp))
dh_dt <- mask(dh_dt, vect(basin_shp), overwrite=TRUE)
dh_dt <- crop(dh_dt, vect(basin_shp))*1000 # from m to mm
dh_dt <- dh_dt / 1.091 # from ice to water
dh_dt[is.na(dh_dt)] <- 0

basin_shp$glacier_dhdt <- round(exact_extract(dh_dt, basin_shp, "mean"), 1)
basin_shp$glacier_dhdt[is.na(basin_shp$glacier_dhdt)] <- 0

## c) Climatological -------------------------------------------------------------------

# p_mean from PMET v1.1
pp_stack <- rast(paste0(path_pmet, "Zenodo/v11/PP_PMETsim_1980_2020_v11d.nc"))
terra::time(pp_stack) <- as.POSIXct(time(pp_stack), tz= "UTC") 
pp_stack <- subset(pp_stack, which(time(pp_stack) > period[1] & time(pp_stack) <= period[2]))

stack <- mean(tapp(pp_stack, strftime(time(pp_stack), format="%Y"), fun = sum, na.rm = TRUE))
basin_shp$p_mean_PMET <- round(exact_extract(stack, basin_shp, "mean"), 1)

stack <- tapp(pp_stack, strftime(time(pp_stack), format="%Y-%m"), fun = sum, na.rm = TRUE)
stack <- tapp(stack, rep(seq(1,12),30), fun = mean, na.rm = TRUE)
stack <- app(stack^2, sum) * 100 / app(stack, sum)^2
basin_shp$p_pci_PMET <- round(exact_extract(stack, basin_shp, "mean"), 1)

# pet_mean from PMET v1.1 (Hargreaves eq based on temp)
pet_stack <- rast(paste0(path_pmet, "Data/Evaporation/Ep_PMETsim_1980_2020d_dev.nc"))
terra::time(pet_stack) <- as.POSIXct(time(pet_stack), tz= "UTC") 
pet_stack <- subset(pet_stack, which(time(pet_stack) > period[1] & time(pet_stack) <= period[2]))
stack <- mean(tapp(pet_stack, strftime(time(pet_stack),format="%Y"), fun = sum, na.rm = TRUE))
basin_shp$pet_mean_PMET <- round(exact_extract(stack, basin_shp, "mean"), 1)
basin_shp$aridity_PMET  <- round(basin_shp$p_mean_PMET/basin_shp$pet_mean_PMET, 3)

# solid pp (< 0◦C)
t2m_stack <- rast(paste0(path_pmet, "Zenodo/v11/Tavg_PMETsim_1980_2020_v11d.nc"))
terra::time(t2m_stack) <- as.POSIXct(time(t2m_stack), tz= "UTC") 
t2m_stack <- subset(t2m_stack, which(time(t2m_stack) > period[1] & time(t2m_stack) <= period[2]))
t2m_stack <- resample(t2m_stack, pp_stack, method = "near", threads = 24)
stack <- ifel(t2m_stack > 0, 0, pp_stack)
stack <- mean(tapp(stack, strftime(time(stack), format="%Y"), fun = sum, na.rm = TRUE))
basin_shp$snow_PMET <- round(exact_extract(stack, basin_shp, "mean"), 1)
basin_shp$frac_snow_PMET     <- round(basin_shp$snow_PMET/basin_shp$p_mean_PMET, 3)

# pp freq
stack <- sum(pp_stack < 1, na.rm = TRUE)
stack <- nlyr(pp_stack) / stack
basin_shp$low_prec_freq_PMET <- round(exact_extract(stack, basin_shp, "mean"), 3)

pp_stack_mean <- mean(pp_stack)
stack <- sum(pp_stack > pp_stack_mean*5, na.rm = TRUE)
stack <- nlyr(pp_stack) / stack
basin_shp$high_prec_freq_PMET <- round(exact_extract(stack, basin_shp, "mean"), 3)

## d) Hydrological ------------------------------------------------------------
stack  <- rast(paste0(path_pmet, "GIS/WB_DGA_2020_Qrunoff.tif"))
basin_shp$Qrunoff_DGA <- exact_extract(stack,   basin_shp, "mean")

# specific problem in Y00003508
basin_shp[is.na(basin_shp)] <- 0

writeVector(vect(basin_shp), "GIS/Basins_Patagonia_all_data.gpkg", overwrite=TRUE)
write.csv(st_drop_geometry(basin_shp), "data/Attributes_all_basins.csv", row.names = FALSE)
write.csv(st_drop_geometry(basin_shp), "data/Attributes_all_basins_pmet.csv", row.names = FALSE)
