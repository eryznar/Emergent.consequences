# PURPOSE: 
# Example workflow for querying crabpack data

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

# GET TANNER SURVEY SPECIMEN/CPUE/BIOABUND DATA FROM CRABPACK ----
# Pull specimen data
species <- "TANNER"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             #years = years, # set in libs/params script
                                             channel = channel)

saveRDS(specimen_data, "./Data/tanner_survey_specimenEBS.rda")

# Calculate CPUE for large/small males
tanner_cpue <- calc_cpue(crab_data = specimen_data,
                       species = species,
                       crab_category = c("large_female", "small_male"))