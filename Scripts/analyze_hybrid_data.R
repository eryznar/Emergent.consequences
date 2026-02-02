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

saveRDS(hybrid_data, "./Data/hybrid_specimen.rda")

snow_data <- crabpack::get_specimen_data(species = "SNOW",
                                         region = "EBS", # can also include NBS
                                         channel = channel)

saveRDS(snow_data, "./Data/snow_specimen.rda")


tanner_data <- crabpack::get_specimen_data(species = "TANNER",
                                         region = "EBS", # can also include NBS
                                         channel = channel)

saveRDS(tanner_data, "./Data/tanner_specimen.rda")


# 6) Calculate per-station CPUE by 1mm size bin
hybrid_cpue <- crabpack::calc_cpue(crab_data = hybrid_data,
                            species = "HYBRID",
                            region = "EBS",
                            bin_1mm = TRUE,
                            size_min = NULL,
                            size_max = NULL,
                            sex = "all_categories")

snow_cpue <- crabpack::calc_cpue(crab_data = snow_data,
                                   species = "SNOW",
                                   region = "EBS",
                                   bin_1mm = TRUE,
                                   size_min = NULL,
                                   size_max = NULL,
                                   sex = "all_categories")

tanner_cpue <- crabpack::calc_cpue(crab_data = tanner_data,
                                   species = "TANNER",
                                   region = "EBS",
                                   years = c(1980:2025),
                                   bin_1mm = TRUE,
                                   size_min = NULL,
                                   size_max = NULL,
                                   sex = c("all_categories"))

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
female <- crabpack::calc_bioabund(crab_data = hybrid_data,
                                   species = "HYBRID",
                                   region = "EBS",
                                   size_min = 70,
                                   crab_category = "all_categories",
                                   sex = "female")

male <- crabpack::calc_bioabund(crab_data = hybrid_data,
                                  species = "HYBRID",
                                  region = "EBS",
                                  size_min = 80,
                                  crab_category = "all_categories",
                                  sex = "male")

write.csv(rbind(female, male), "./Data/hybrid_bioabund.csv")
          
          
# plot
ggplot(hybrid_bioabund, aes(YEAR, ABUNDANCE, color = CATEGORY))+
  geom_line(linewidth = 1) +
  theme_bw()

# Analyze size distribution
hybrid_spec <- readRDS("./Data/hybrid_specimen.rda")$specimen

pp <- hybrid_spec %>% filter(YEAR == 2025)

ggplot()+
  geom_histogram(pp %>% filter(SEX != 4), mapping = aes(SIZE)) +
  facet_wrap(~SEX, nrow = 2, scales = "free")+
  theme_bw()

# Fit models


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

hybrid_abund <- read.csv("./Data/hybrid_bioabund.csv") %>%
  filter(CATEGORY %in% c("mature_female", "large_male")) %>%
  filter(YEAR > 1989) %>%
  group_by(YEAR) %>%
  reframe(ABUNDANCE = sum(ABUNDANCE))


snow_abund <- read.csv("./Data/snow_bioabund.csv") %>%
  filter(YEAR > 1989) %>%
  group_by(YEAR) %>%
  reframe(ABUNDANCE = sum(ABUNDANCE))

tanner_abund <- read.csv("./Data/tanner_bioabund.csv") %>%
  filter(YEAR > 1989) %>%
  group_by(YEAR) %>%
  reframe(ABUNDANCE = sum(ABUNDANCE))


overlap <- rbind(read.csv("./Data/bhatt_malesnow.femTanner.csv"),
                 read.csv("./Data/bhatt_maleTanner.femsnow.csv")) %>%
  filter(YEAR > 1989) %>%
  na.omit()

ice <- read.csv("./Data/ice_means_1989-2025.csv") %>%
  group_by(year) %>%
  reframe(value = mean(value)) %>%
  filter(year> 1989 & year != 2020)

lag = 3
mod.dat <- data.frame(year = unique(hybrid_abund$YEAR), 
                      ice = lag(ice$value,6),
                      hybrid_abund = lag(hybrid_abund$ABUNDANCE/1e6),
                      snow_abund = lag(snow_abund$ABUNDANCE/1e6,3),
                      tanner_abund = lag(tanner_abund$ABUNDANCE/1e6, 3),
                      # hybrid_matfem = (hybrid_abund$ABUNDANCE[hybrid_abund$CATEGORY == "mature_female"])/1e6,
                      # hybrid_lgmale = (hybrid_abund$ABUNDANCE[hybrid_abund$CATEGORY == "large_male"])/1e6,
                      # snow_matfem = (snow_abund$ABUNDANCE[snow_abund$CATEGORY == "mature_female"])/1e6,
                      # snow_lgmale = (snow_abund$ABUNDANCE[snow_abund$CATEGORY == "large_male"])/1e6,
                      # tanner_matfem = (tanner_abund$ABUNDANCE[tanner_abund$CATEGORY == "mature_female"])/1e6,
                      # tanner_lgmale = (tanner_abund$ABUNDANCE[tanner_abund$CATEGORY == "large_male"])/1e6,
                      bhatt_snowmale_tannerfem = lag(overlap$bhatt[overlap$comparison == "snowlgmale-tannermatfem"], 0),
                      bhatt_tannermale_snowfem = lag(overlap$bhatt[overlap$comparison == "tannerlgmale-snowmatfem"]), 0)

library(corrplot)

# Select numeric columns and compute correlation
cor_vars <- mod.dat %>%
  dplyr::select(ice, snow_abund, tanner_abund, 
                bhatt_snowmale_tannerfem, bhatt_tannermale_snowfem) %>%
  na.omit()  # drop rows with any NA

cor_mat <- cor(cor_vars, use = "complete.obs")

# Plot
corrplot(cor_mat, 
         method = "circle",
         type = "upper",
         order = "hclust",
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black",
         number.cex = 0.7)



hybrid.gam <- gam(
  hybrid_abund ~  s(bhatt_snowmale_tannerfem, k=4) +
    s(bhatt_tannermale_snowfem, k = 4) + 
    s(ice, k=4) + 
    s(tanner_abund, k = 4)+
    s(snow_abund, k = 4),
  data = mod.dat,
  family = nb()
)

acf(residuals(hybrid.gam))

gam.check(hybrid.gam)
summary(hybrid.gam)

res_dharma <- simulateResiduals(hybrid.gam, n = 1000)
plot(res_dharma)

