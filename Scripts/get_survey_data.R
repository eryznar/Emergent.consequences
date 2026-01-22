# PURPOSE: 
# To get updated survey data by year for chela maturity processing

# Notes:
# 1) What is the best way to update data? Use the chela database and then specimen data each new year? Or will
# the chela database be updated in time?

# LOAD LIBS/PARAMS -----
source("./Scripts/load_libs_params.R")

# Set channel
channel <- "API"

# GET SNOW SURVEY SPECIMEN/CPUE DATA FROM CRABPACK ----
# Pull specimen data
species <- "SNOW"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             #years = years, # set in libs/params script
                                             channel = channel)

saveRDS(specimen_data, "./Data/snow_survey_specimenEBS.rda")

# Calculate CPUE for mature/immature female
snow_spec <- readRDS("./Data/snow_survey_specimenEBS.rda")

snow_cpue <- calc_cpue(crab_data = specimen_data,
                       species = species,
                       crab_category = c("mature_female", "immature_female"))

# Calculate CPUE for mature/immature male
snow_mod <- readRDS("Y:/KOD_Research/Ryznar/Crab functional maturity/SNOW/sdmTMB/sdmTMB_spVAR_SIZE_k300_BEST.rda")



# GET TANNER SURVEY SPECIMEN/CPUE/BIOABUND DATA FROM CRABPACK ----
# Pull specimen data
species <- "TANNER"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             #years = years, # set in libs/params script
                                             channel = channel)

saveRDS(specimen_data, "./Data/tanner_survey_specimenEBS.rda")

# Calculate CPUE for mature/immature male
tanner_mod <- readRDS("Y:/KOD_Research/Ryznar/Crab functional maturity/TANNER/sdmTMB/sdmTMB_spVAR_SIZE_k200_BEST.rda")
