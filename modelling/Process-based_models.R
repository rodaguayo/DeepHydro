# GR4J-CemaNeige and TUWmodel + aux functions

library("hydroGOF")
library("TUWmodel")
library("hydromad")
library("airGR")

metric_names <- c("KGE", "r", "Beta", "Gamma", "NSE",  "FMS", "FLV", "FHV")
param_names_TUWmodel <-  c("SCF", "DDF", "Tr",  "Ts",  "Tm", "LPrat",  "FC", 
                           "BetaP", "k0", "k1", "k2", "lsuz", "cperc", "bmax", "croute")

param_names_GR4J <- c("GR4J_X1", "GR4J_X2", "GR4J_X3",  "GR4J_X4", "CemaNeige_X1", "CemaNeige_X2")

## hydroclimate data for each basin 
hydroclimate_basin <- function(basin, obs) {
  
  pp_i    <- as.matrix(pp[grepl(basin, colnames(pp))])
  t2m_i   <- as.matrix(t2m[grepl(basin, colnames(t2m))])
  pet_i   <- as.matrix(pet[grepl(basin, colnames(pet))])
  zones_i <- rep(1/ncol(pp_i), ncol(pp_i))
  area_i  <- q_metadata[basin,]$total_area
  
  q_glacier_i <- as.numeric(q_glacier[,basin])
  q_glacier_i <- (q_glacier_i*1000*86400)/(area_i*10^6)
  
  if (obs) {
  q_obs_i     <- as.numeric(q_obs[,basin]) 
  q_obs_i     <- (q_obs_i*1000*86400)/(area_i*10^6)
  
  return(list(q_obs_i = q_obs_i, q_glacier_i = q_glacier_i, 
              pp_i = pp_i, t2m_i = t2m_i, pet_i = pet_i, zones_i = zones_i))
  } else {
    return(list(q_glacier_i = q_glacier_i, pp_i = pp_i, t2m_i = t2m_i, 
                pet_i = pet_i, zones_i = zones_i))
  }
}

## TUWmodel
hydro_TUWmodel <- function(param.values, pp_i, t2m_i, pet_i, q_glacier_i, obs = NULL, zones_i, calibration) {
  
  q_sim_i <- TUWmodel(prec = pp_i, airt = t2m_i, ep = pet_i, area = zones_i, param=as.numeric(param.values))
  q_sim_i <- as.numeric(q_sim_i$q)[period %in% period_calib]
  
  q_glacier_i <- TUWmodel(prec = as.numeric(q_glacier_i), param=as.numeric(param.values),
                          airt =  rep(30, length(period_calib)), ep =  rep(0, length(period_calib)))
  q_sim_i <- q_sim_i + as.numeric(q_glacier_i$q)
  
  if (is.null(obs)) {return(list(q_sim = q_sim_i))}
  
  NSE_i   <- NSE(sim=q_sim_i, obs=obs, epsilon="none", out.type="full", na.rm=TRUE)
  if (calibration == TRUE) { return(list("GoF" = NSE_i, "sim" = q_sim_i)) }
  else {
    KGE_i   <- KGE(sim=q_sim_i, obs=obs, method="2012", out.type="full", na.rm=TRUE)
    fdc_fms_i <- fdc_fms(obs = obs, sim = q_sim_i, lower = 0.2, upper = 0.7)
    fdc_flv_i <- fdc_flv(obs = obs, sim = q_sim_i, l = 0.9)
    fdc_fhv_i <- fdc_fhv(obs = obs, sim = q_sim_i, h = 0.02)
    metrics_i   <- as.numeric(c(KGE_i$KGE.value, KGE_i$KGE.elements, NSE_i, fdc_fms_i, fdc_flv_i, fdc_fhv_i))
    return(TUWmodel_i <- list(metrics = metrics_i, q_sim = q_sim_i))
    }
}

## GR4J
precipitation_fraction_df <- function(temperature_df) {
  snow_fraction_df <- temperature_df
  snow_fraction_df[temperature_df <= 0] <- 1
  snow_fraction_df[temperature_df >= 2] <- 0
  between_indices <- temperature_df > 0 & temperature_df < 2
  snow_fraction_df[between_indices] <- (2 - temperature_df[between_indices]) / 2
  
  return(snow_fraction_df)
}

hydro_GR4J <- function(param.values, pp_i, t2m_i, pet_i, q_glacier_i, obs = NULL, zones_i, calibration) {
  
  pp_i    <- as.list(as.data.frame(pp_i))
  names(pp_i) <- paste0("L", seq(1, length(pp_i), 1))
  t2m_i   <- as.list(as.data.frame(t2m_i))
  names(t2m_i) <- paste0("L", seq(1, length(t2m_i), 1))
  
  InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, 
                                   DatesR = period, verbose = F,
                                   NLayers   = length(zones_i), # number of layers
                                   HypsoData = rep(0, 101), ZInputs = 0, 
                                   Precip    = rowMeans(as.data.frame(pp_i)), 
                                   TempMean  = rowMeans(as.data.frame(t2m_i)),
                                   PotEvap   = rowMeans(pet_i)) 
  
  InputsModel$LayerPrecip   <- pp_i
  InputsModel$LayerTempMean <- t2m_i
  t2m_i <- as.data.frame(t2m_i)
  InputsModel$LayerFracSolidPrecip <- as.list(precipitation_fraction_df(t2m_i))
  
  RunOptions  <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J, InputsModel = InputsModel, warnings = F,
                                  IsHyst = FALSE,
                                  Outputs_Cal = c("Qsim", "PliqAndMelt"), 
                                  IndPeriod_Run    = seq(1, length(period))[period %in% period_calib], 
                                  IndPeriod_WarmUp = seq(1, length(period))[period %in% period_warm_up])  
  
  q_sim_i   <- RunModel_CemaNeigeGR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = as.numeric(param.values))$Qsim
  q_glacier_i <- gr4jrouting.sim(q_glacier_i, x2 = param.values[["GR4J_X2"]], x3 = param.values[["GR4J_X3"]], x4 = param.values[["GR4J_X4"]], split = 0.9)
  q_sim_i   <- q_sim_i + q_glacier_i
  
  if (is.null(obs)) {return(list(q_sim = q_sim_i))}
  
  NSE_i   <- NSE(sim=q_sim_i, obs=obs, epsilon="none", out.type="full", na.rm=TRUE)
  if (calibration == TRUE) { return(list("GoF" = NSE_i, "sim" = q_sim_i)) }
  else {
  KGE_i   <- KGE(sim=q_sim_i, obs=obs, method="2012", out.type="full", na.rm=TRUE)
  fdc_fms_i <- fdc_fms(obs = obs, sim = q_sim_i, lower = 0.2, upper = 0.7)
  fdc_flv_i <- fdc_flv(obs = obs, sim = q_sim_i, l = 0.9)
  fdc_fhv_i <- fdc_fhv(obs = obs, sim = q_sim_i, h = 0.02)
  metrics_i   <- as.numeric(c(KGE_i$KGE.value, KGE_i$KGE.elements, NSE_i, fdc_fms_i, fdc_flv_i, fdc_fhv_i))
  return(GR4J_i <- list(metrics = metrics_i, q_sim = q_sim_i))
  }
}

