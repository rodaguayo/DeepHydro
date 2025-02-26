# Create glacier dataset based on RGI7 ------------------------------------------------------------

rm(list=ls())
cat("\014")  

library("exactextractr")
library("stringr")
library("terra")
library("dplyr")
library("sf")

setwd("/home/rooda/OneDrive/Projects/DeepHydro/")
path_pmet <- "/home/rooda/OneDrive/Projects/PatagoniaMet/"
path_data <- "/home/rooda/Datasets/GLACIERS/"

# read file from catchment_delimitation.R
all_basins <- st_read("data/GIS/Basins_Patagonia_all.gpkg")
all_basins <- st_transform(all_basins, 32718) # UTM 18S for the glaciers

# read file from catchment_delimitation.R
pmet_basins <- st_read(paste0(path_pmet, "GIS/Basins_PMETobs_int_dev.shp"))
pmet_basins$gauge_code <- read.csv(paste0(path_pmet, "Zenodo/v11/Q_PMETobs_v11_metadata.csv"))$gauge_code
pmet_basins <- st_transform(pmet_basins, 32718) # UTM 18S for the glaciers

# dem file
dem <- rast(paste0(path_pmet, "GIS/dem_patagonia3f.tif"))
dem <- project(dem, "epsg:32718", method = "bilinear")
dem <- subst(dem, NA, 0)  # NAs to sea level (= 0)

# 1. RGI glaciers: Selection and basin_id assignment ----------------------------------------------
RGI  <- st_read(paste0(path_data, "RGI/RGI7_17.shp"))
RGI  <- subset(RGI, RGI$cenlat < -40.7)
RGI  <- subset(RGI, RGI$cenlon < -68.35) # PMET extent
RGI  <- RGI[,c("rgi_id", "src_date", "cenlon", "cenlat")]
RGI  <- st_transform(RGI, 32718) # UTM 18S
RGI$area_km2 <- expanse(vect(RGI), unit="km")

## assignment based on terminus location 
RGI_c <- extract(dem, as.lines(vect(RGI)), xy = TRUE)
RGI_c <- RGI_c %>% group_by(ID) %>% slice_min(order_by = dem_patagonia3f)
RGI_c <- aggregate(RGI_c, by = list(RGI_c$ID), FUN = mean) # mean if there area several options
RGI_c <- st_as_sf(RGI_c, coords = c("x","y"), crs = 32718)
RGI_c <- st_join(RGI_c["ID"], all_basins[,c("gauge_id", "ID_Zone", "geometry")], join = st_within) 
RGI$ID_basin <- RGI_c$gauge_id 
RGI$ID_Zone <- RGI_c$ID_Zone 
sum(!is.na(RGI$ID_basin))

### the same but for PMET
RGI_c <- extract(dem, as.lines(vect(RGI)), xy = TRUE)
RGI_c <- RGI_c %>% group_by(ID) %>% slice_min(order_by = dem_patagonia3f)
RGI_c <- aggregate(RGI_c, by = list(RGI_c$ID), FUN = mean) # mean if there area several options
RGI_c <- st_as_sf(RGI_c, coords = c("x","y"), crs = 32718)
RGI_c <- st_join(RGI_c["ID"], pmet_basins[,c("gauge_code", "geometry")], join = st_within) 
RGI$ID_PMET <- RGI_c$gauge_code 
RGI$ID_PMET <- as.character(RGI$ID_PMET)

## assignment based on centroid location
RGI_c <- st_centroid(RGI["geometry"]) 
RGI_c <- st_join(RGI_c, all_basins[,c("gauge_id", "ID_Zone", "geometry")], join = st_within)
RGI$ID_basin[is.na(RGI$ID_basin)] <- RGI_c$gauge_id[is.na(RGI$ID_basin)]  # Assign catchment
RGI$ID_Zone[is.na(RGI$ID_Zone)] <- RGI_c$ID_Zone[is.na(RGI$ID_Zone)]

# select glaciers with a id_basin
RGI <- subset(RGI, !is.na(RGI$ID_basin)) 

# 2. Glacier attributes ---------------------------------------------------------------------------

## 2.1 Volume data from and Millan et al. 2022 ----------------------------------------------------
vol_M22  <- rast(paste0(path_pmet, "GIS/Glaciers/Thickness_2022.tif"))

RGI$vol_M22  <- exact_extract(vol_M22,  RGI, "sum") * res(vol_M22)[1] * res(vol_M22)[2] * 1e-9
RGI$vol_M22c <- round(exact_extract(not.na(vol_M22),  RGI, "mean"), 3)
RGI$vol_M22[RGI$vol_M22c < 0.5] <- NA
RGI$vol_M22[RGI$vol_M22 == 0]   <- NA

for (i in sort(unique(RGI$ID_Zone))) { # VAS based on M22 data 
  model <- lm(log(RGI$vol_M22)[RGI$ID_Zone == i & RGI$vol_M22c > 0.9] ~ log(RGI$area_km2)[RGI$ID_Zone == i & RGI$vol_M22c > 0.9] )
  RGI$vol_M22[RGI$ID_Zone == i & is.na(RGI$vol_M22)] <- exp(coef(model)[1]) * RGI$area_km2[RGI$ID_Zone == i & is.na(RGI$vol_M22)] ** coef(model)[2]
}

## 2.2 dmdtda: specific-mass change rate in meters water-equivalent per year ----------------------
dhdt_21   <- rast(paste0(path_pmet, "GIS/Glaciers/dhdt_2021.tif"))

RGI$dmdtda_21  <- exact_extract(dhdt_21,  RGI, "mean") * 0.850 
RGI$dmdtda_21c <- round(exact_extract(not.na(dhdt_21),  RGI, "mean"), 3)
RGI$dmdtda_21[RGI$dmdtda_21c < 0.9] <- NA

dmdtda_RGI_mean <- sapply(split(RGI, RGI$ID_Zone), function(d) weighted.mean(d$dmdtda_21, w = d$area_km2, na.rm = T))

for (i in sort(unique(RGI$ID_Zone))) {
  RGI$dmdtda_21[RGI$ID_Zone == i & is.na(RGI$dmdtda_21)] <- dmdtda_RGI_mean[i]
}

# 3. Cook RGI df ---------------------------------------------------------------------------------

RGI <- st_transform(RGI, 4326) # go back to WGS84 (RGI format)

# source year in RGI format
RGI$BgnDate <- paste0(substr(RGI$src_date, 0,4), substr(RGI$src_date, 6,7), substr(RGI$src_date, 9,10)) 

# for OGGM
RGI["GLIMSId"]  <- '' 
RGI['O1Region'] <- '17'
RGI['O2Region'] <- '02'
RGI['RGIId'] <- paste0("RGI70-17.", substr(RGI$rgi_id, 19,24)) 
RGI['vol_M22c'] <- NULL
RGI['dmdtda_21c'] <- NULL

# save using terra (problems in sf)
st_write(st_zm(RGI), "data/GIS/RGI7_Hydro.gpkg", append = FALSE)
