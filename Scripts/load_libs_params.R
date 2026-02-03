library(tidyverse)
library(odbc)
library(DBI)
library(getPass)
library(keyring)
library(lifecycle)
library(data.table)
library(crabpack)
library(sf)
library(gstat)
library(rnaturalearth)
library(raster)
library(concaveman)
library(png)
library(mgcv)
library(dplyr)
library(mgcv)
library(ggplot2)
library("rnaturalearth")
library(patchwork)
library(gratia)
library(MuMIn)
library(DHARMa)
library(mgcViz)
library(akgfmaps)
library(dlm)
#library(INLA)
# library(sdmTMB)
# library(glmmTMB)
# library(broom)
# library(gstat)
# library(devtools)
# #install_github("vast-lib/tinyVAST", dependencies = TRUE)
# library(tinyVAST)
library(ecmwfr)
library(tidync)

# Set years
current.year <- 2025
years <- c(1989:2007, 2009:2011, 2013, 2015, 2017:2019, 2021:current.year)

# Specify directory
dir <- "Y:/KOD_Research/Ryznar/Crab functional maturity"

data_dir <- "Y:/KOD_Survey/EBS Shelf/Data_Processing/Data/" # for survey data

remote_dir <- "Y:/KOD_Research/Ryznar/Crab functional maturity/"

# CRS for spatial blocking
region_layers <- akgfmaps::get_base_layers("sebs")

map.crs <- region_layers$crs

in.crs = "+proj=longlat +datum=NAD83"

crs.latlon <- "epsg:4326" 

# Read in spatial layers
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

# Set coordinate reference system
ncrs <- "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# Load prediction grid
ebs_grid <- read.csv(here::here(paste0(dir, "/ebs_coarse_grid.csv")))


# FUNCTIONS ----
diagnose <- function(model){
  model.name <- deparse(substitute(model))
  acf(resid(model))
  ss <- summary(model) # model summary
  
  gam.check(model) # make sure smooth terms are appropriate
  
  plot(simulateResiduals(model)) # Checks uniformity, dispersion, outliers via DHARMa
  
  # plot facetted smooths
  sm.dat <- smooth_estimates(model) %>%
    pivot_longer(., cols = 6:ncol(.), names_to = "resp", values_to = "value")
  
  ggplot(sm.dat, aes(x = value, y = .estimate)) +
    geom_ribbon(sm.dat, mapping = aes(ymin = .estimate + 2 * .se, ymax = .estimate - 2 * .se), fill = "cadetblue", alpha = 0.25)+
    geom_line(color = "cadetblue", linewidth = 1.25) +
    facet_wrap(~ .smooth, scales = "free_x", ncol = 2) +   # facet by smooth term name
    theme_bw()+
    ggtitle(model.name)+
    ylab("Partial effect")+
    xlab("Value")+
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14)) -> sm.plot
  
  return(list(ss, sm.plot))
}

diagnose.gamm <- function(model){
  lme.part <- model$lme
  rr.lme <- resid(lme.part, type = "normalized")
  acf(rr.lme)
  
  model <- model$gam
  #model.name <- deparse(substitute(model))
  
  ss <- summary(model) # model summary
  
  gam.check(model) # make sure smooth terms are appropriate
  
  dev.explained = (model$null.deviance - model$deviance) / model$null.deviance * 100
  
  #plot(simulateResiduals(model)) # Checks uniformity, dispersion, outliers via DHARMa
  
  # plot facetted smooths
  sm.dat <- smooth_estimates(model) %>%
    dplyr::select(!c(.type, .by)) %>%
    pivot_longer(!c(1:3), names_to = "Parameter", values_to = "value") %>%
    filter(is.na(value) == FALSE) %>%
    group_by(.smooth) %>%
    mutate(
      lower = .estimate - 1.96 * .se,
      upper = .estimate + 1.96 * .se
    ) %>%
    ungroup()
  

  
  ggplot(sm.dat, aes(x = value, y = .estimate)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "cadetblue", alpha = 0.25) +
    geom_line(color = "cadetblue", linewidth = 1.25) +
    facet_wrap(~ .smooth, scales = "free_x", ncol = 2)+
    #ggtitle(model.name) +
    ylab("Partial effect") +
    theme_bw() +
    xlab("Value") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14)) -> sm.plot
  
  
  return(list(ss, sm.plot, deviance.explained = dev.explained))
}

# # DATA ----
# # Snow minima and cutline coefficients
# snow_minima <- read.csv("./Output/chela_cutline_minima.csv") %>%
#   filter(SPECIES == "SNOW") %>%
#   mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
#          BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])
# 
# snow_BETA0 <- unique(snow_minima$BETA0)
# snow_BETA1 <- unique(snow_minima$BETA1)
# 
# # Snow chela data compiled by Shannon
# snow_chela <- read.csv(paste0(data_dir, "specimen_chela.csv")) %>% # already != HT 17, only shell 2, no special projects
#   filter(SPECIES == "SNOW", HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE,
#          YEAR %in% years) %>% # filter for males, sh2, only chela msrd, not HT17
#   mutate(ratio = SIZE/CHELA_HEIGHT,
#          LN_CH = log(CHELA_HEIGHT),
#          LN_CW = log(SIZE),
#          CW = SIZE) %>%
#   filter(ratio > 2 & ratio < 35) %>% # filter extreme measurements
#   dplyr::select(!c(ratio)) %>%
#   mutate(cutoff = snow_BETA0 + snow_BETA1*LN_CW, # apply cutline model
#          MATURE = case_when((LN_CH > cutoff) ~ 1,
#                             TRUE ~ 0))
# 
# # Tanner minima data
# tanner_minima <- read.csv("./Output/chela_cutline_minima.csv") %>%
#   filter(SPECIES == "TANNER") %>%
#   mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
#          BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])
# 
# 
# tanner_BETA0 <- unique(tanner_minima$BETA0)
# tanner_BETA1 <- unique(tanner_minima$BETA1)
# 
# # Tanner chela data compiled by Shannon
# tanner_chela <- read.csv(paste0(data_dir, "specimen_chela.csv")) %>% # already != HT 17, only shell 2, no special projects
#   filter(SPECIES == "TANNER", HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE,
#          YEAR %in% years) %>% # filter for males, sh2, only chela msrd, not HT17
#   mutate(ratio = SIZE/CHELA_HEIGHT,
#          LN_CH = log(CHELA_HEIGHT),
#          LN_CW = log(SIZE),
#          CW = SIZE) %>%
#   filter(ratio > 2 & ratio < 35) %>% # filter extreme measurements
#   dplyr::select(!c(ratio)) %>%
#   mutate(cutoff = tanner_BETA0 + tanner_BETA1*LN_CW, # apply cutline model
#          MATURE = case_when((LN_CH > cutoff) ~ 1,
#                             TRUE ~ 0))
