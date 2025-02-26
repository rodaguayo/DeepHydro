# Calculate equal-area elevation bands -----------------------------------------------------------


by.areas <- function(x, n.zones){ # .byAreas
  
  sorted <- sort(x)
  range <- length(sorted) / n.zones
  
  cuts <- sorted[c(1, 1:n.zones * range)]
  cuts <- round(cuts, 1)
  
  median.elevation <- zoo::rollapply(cuts, width = 2, by = 1, median)
  
  table <- data.frame(No.zone = 1:n.zones, Min = cuts[1:n.zones],
                      Max = cuts[2:(n.zones+1)], Median = median.elevation)
  
  cut.values <- cut(x, breaks = cuts)
  
  vals <- as.numeric(ftable(cut.values))
  
  number.values <- length(which(!is.na(cut.values)))
  area <- round(vals / number.values, 5)
  
  return(list("cuts" = cuts,"median.elevation"= median.elevation,
              "area" = area, "table" = table))
  
} 

elevationZones <- function(x, dem, max.zones = 10,
                           min.elevZ = 200, elev.thres = NULL) {
  
  # Masking the DEM to the catchment boundary
  if (terra::ext(x) != terra::ext(dem)) dem <-
      crop(dem, ext(x), snap = "OUT")
  dem <- terra::mask(dem, x)
  
  # Getting the minimum and maximum elevations
  elev.min <- min(terra::values(dem), na.rm = TRUE)
  elev.max <- max(terra::values(dem), na.rm = TRUE)
  
  ### Here should enter the separation of  cases
  #  1) elev.min > elev.thres and (elev.diff / min.elevZ) < max.zones -> @200 m
  #  2) elev.min > elev.thres and (elev.diff / min.elevZ) > max.zones ->
  #     all elevation zones by area
  #  3) elev.min < elev.thres < elev.max ->
  #     elev.zones + one from threshold to min
  
  # Checking the number of elevation zones for the catchment
  elev.diff <- elev.max - elev.min
  
  if(is.null(elev.thres)) 
    elev.thres <- -1
  
  if( elev.thres < elev.min |
      elev.thres > elev.max |
      (elev.thres - elev.min) < min.elevZ ){
    
    if ( (elev.diff / min.elevZ) < max.zones ){
      # CASE 1, elev.thres = NULL
      max.zones   <- ceiling(elev.diff / min.elevZ)  
      elev.ranges <- by.areas(values(dem), max.zones)
      
    } else {
      # CASE 2, elev.thres = NULL
      elev.ranges <- by.areas(values(dem), max.zones)
      
    }
    
  } else {
    min.zonediff <- elev.thres - elev.min
    # Excluding first zone
    dem.rest                         <- dem
    dem.rest[dem.rest <= elev.thres] <- NA
    
    # Checking the elevation difference fron threshold up
    elev.diff <- elev.max - elev.thres
    
    # Getting the elevation zones from the threshold elevation up
    if ( (elev.diff / min.elevZ) < max.zones - 1 ){
      # CASE 1, elev.thres = NULL
      max.zones   <- ceiling(elev.diff / min.elevZ)  
      elev.ranges <- by.areas(values(dem.rest), max.zones)
      
    } else {
      # CASE 2, elev.thres = NULL
      elev.ranges <- by.areas(values(dem.rest), max.zones - 1)
      
    }

    # Adding first zone to cuts from results
    elev.ranges$cuts <-c(elev.min, elev.ranges$cuts)
    # Calculating median elevation of first zone
    elev.first                          <- dem
    elev.first[elev.first > elev.thres] <- NA
    
    elev.med <- median(values(elev.first), na.rm = TRUE)
    elev.ranges$median.elevation <-c(elev.med, elev.ranges$median.elevation)
    
    # Calculating area of first zone
    area.others <- length(which(!is.na(values(dem.rest))))
    area.new    <- length(which(!is.na(values(elev.first))))
    area.ratio  <- area.new / area.others
    area.ratio  <- c(area.ratio, elev.ranges$area)
    area.ratio  <- area.ratio / sum(area.ratio)
    
    elev.ranges$area <- area.ratio
    
    # Calculating first zone in the table
    elev.ranges$table <-
      data.frame(No.zone = 1:length(elev.ranges$median.elevation),
                 Min = c(elev.min, elev.ranges$table$Min),
                 Max = c(elev.ranges$cuts[2], elev.ranges$table$Max),
                 Median = elev.ranges$median.elevation)
    
  }
  
  # Rasterizing the elevation zones
  elev.zones <- terra::classify(dem, elev.ranges$cuts, include.lowest=TRUE, brackets=TRUE)
  
  # Function to create a mask per elevation zone
  rasterize.zones <- function(el.zone, dem, elev.zones){
    
    rast         <- dem
    values(rast) <- NA
    
    rast[which(values(elev.zones) == el.zone)] <- 1
    
    return(rast)
  }
  
  n.zones <- max(elev.ranges$table$No.zone)
  zones.stack <- elev.zones

  elev.ranges$table$Difference <- elev.ranges$table$Max - elev.ranges$table$Min
  
  # Adding the area in the table
  elev.ranges$table$Area.perc <- round(elev.ranges$area * 100, 2)
  
  # Creating the final output
  elev.ranges$zonesRaster <- zones.stack
  
  return(elev.ranges)
  
}
