source("./Scripts/load_libs_params.R")

# 1) install crabpack
#devtools::install_github("AFSC-Shellfish-Assessment-Program/crabpack")

# 2) load crabpack
library("crabpack")
library("dplyr")
library("tidyverse")

# 3) Set channel
channel <- "API"

# 4) Set species
species <- "HYBRID"

# 5) Pull specimen data
hybrid_data <- crabpack::get_specimen_data(species = "HYBRID",
                                             region = "EBS", # can also include NBS
                                             channel = channel)

snow_data <- crabpack::get_specimen_data(species = "SNOW",
                                         region = "EBS", # can also include NBS
                                         channel = channel)

tanner_data <- crabpack::get_specimen_data(species = "TANNER",
                                         region = "EBS", # can also include NBS
                                         channel = channel)

# 6) Calculate per-station CPUE by 1mm size bin
hybrid_cpue <- crabpack::calc_cpue(crab_data = hybrid_data,
                            species = "HYBRID",
                            region = "EBS",
                            bin_1mm = TRUE,
                            size_min = NULL,
                            size_max = NULL,
                            crab_category = "all_categories")

snow_cpue <- crabpack::calc_cpue(crab_data = snow_data,
                                   species = "SNOW",
                                   region = "EBS",
                                   bin_1mm = TRUE,
                                   size_min = NULL,
                                   size_max = NULL,
                                   crab_category = "all_categories")

tanner_cpue <- crabpack::calc_cpue(crab_data = tanner_data,
                                   species = "TANNER",
                                   region = "EBS",
                                   bin_1mm = TRUE,
                                   size_min = NULL,
                                   size_max = NULL,
                                   crab_category = "all_categories")

# 7) filter by stations where we only find hybrids
positive_catch <- hybrid_cpue %>% 
                    filter(CPUE >0) %>%
                    dplyr::select(STATION_ID, YEAR) %>%
                    distinct() # this selects the unique year X station combinations for positive catch hybrids

# 8) Generature output dataframe that includes proportion hybrid by 1mm size bin

out <- rbind(tanner_cpue, snow_cpue) %>% # bind snow and tanner cpue data frames
        right_join(., positive_catch) %>% # join snow and tanner cpue data frames positive catch hybrid stations; this will only keep data that matches the year X station combination for positive catch hybrid data
        rbind(hybrid_cpue %>% filter(CPUE>0)) %>% # bind with all cpue data for positive catch hybrids
        filter(YEAR %in% c(2015:2025)) %>% # filter data to within the last ten years
          group_by(SIZE_1MM, SEX_TEXT) %>% # group by SIZE and SEX to calculate next step
          mutate(TOTAL_CPUE = sum(CPUE)) %>% # TOTAL CPUE is calculated as CPUE sum by SIZE and SEX across stations and species; mutate() adds a column to the data frame for this
          ungroup() %>% # telling R to no longer use SIZE and SEX as grouping variables
          group_by(SIZE_1MM, SPECIES, SEX_TEXT) %>% # telling R to now group by SIZE, SPECIES, and SEX for next step
          mutate(SPECIES_CPUE = sum(CPUE)) %>% # calculating species-specific CPUE by SEX and SIZE across stations
          ungroup() %>% # telling R to no longer use SIZE, SEX, and SPECIES as grouping variables
          dplyr::select(SPECIES, YEAR, SIZE_1MM, TOTAL_CPUE, SPECIES_CPUE, SEX_TEXT) %>% # only selecting columns of interest now
          distinct() %>% # distinct() removes any duplicate rows
          filter(SPECIES == "HYBRID") %>% # only use hybrid data
          group_by(SIZE_1MM, SEX_TEXT) %>% # grouping by SIZE and SEX again
          reframe(PROP_CPUE = SPECIES_CPUE/TOTAL_CPUE) # summarizing the data using reframe() to calculate proportion hybrid by size and sex

# 8) plot
ggplot(out, aes(SIZE_1MM, PROP_CPUE, color = SEX_TEXT))+
  geom_line() + 
  theme_bw() +
  geom_smooth() 

# 9) calculate biomass and abudance
hybrid_bioabund <- crabpack::calc_bioabund(crab_data = hybrid_data,
                                   species = "HYBRID",
                                   region = "EBS",
                                   size_min = 70,
                                   sex = "female")

# plot
ggplot(hybrid_bioabund, aes(YEAR, ABUNDANCE, color = CATEGORY))+
  geom_line(linewidth = 1) +
  theme_bw()
