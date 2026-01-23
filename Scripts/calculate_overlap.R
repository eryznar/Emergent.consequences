snow_matmale_cpue <- read.csv("./Data/snow_maturecpue.csv")
tanner_matmale_cpue <- read.csv("./Data/tanner_maturecpue.csv")


# females
snow_spec <- readRDS("./Data/snow_specimen.rda")
tanner_spec <- readRDS("./Data/tanner_specimen.rda")


snow_fem_cpue <- crabpack::calc_cpue(crab_data = snow_spec,
                                     region = "EBS",
                                     species = "SNOW",
                                     years = 1980:2025,
                                     crab_category = "mature_female")

tanner_fem_cpue <- crabpack::calc_cpue(crab_data = tanner_spec,
                                     region = "EBS",
                                     species = "TANNER",
                                     years = 1980:2025,
                                     crab_category = "mature_female")


# Overlap between snow mature males and tanner mature females ----
s.male_t.fem <- rbind(
  snow_matmale_cpue %>%
    filter(CATEGORY == "Mature male") %>%
    dplyr::select(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE,
                  CATEGORY, DISTRICT, STRATUM, TOTAL_AREA, CPUE, CPUE_MT, CPUE_LBS),
  tanner_fem_cpue %>%
    dplyr::select(!COUNT)) %>%
  filter(YEAR >= 1990) %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE, CATEGORY, SPECIES) %>%
  reframe(CPUE = sum(CPUE), .groups = "drop_last") %>%  # keep YEAR, STATION_ID, ... groups
  group_by(YEAR, SPECIES) %>%
  mutate(
    TOTAL_CPUE = sum(CPUE),
    P = ifelse(TOTAL_CPUE > 0, CPUE / TOTAL_CPUE, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-CPUE, -TOTAL_CPUE) %>%
  pivot_wider(
    names_from  = CATEGORY,
    values_from = P
  )
  
 
# Split by species
snow_prob <- s.male_t.fem %>%
  filter(SPECIES == "SNOW") %>%
  dplyr::select(YEAR, STATION_ID, snow_p = `Mature male`)

tanner_prob <- s.male_t.fem %>%
  filter(SPECIES == "TANNER") %>%
  dplyr::select(YEAR, STATION_ID, tanner_p = mature_female)

# Align on common YEAR × STATION_ID, fill missing probs with 0
years  <- sort(union(snow_prob$YEAR, tanner_prob$YEAR))
stns   <- sort(union(snow_prob$STATION_ID, tanner_prob$STATION_ID))

p_both <- expand_grid(YEAR = years, STATION_ID = stns) %>%
  left_join(snow_prob,   by = c("YEAR", "STATION_ID")) %>%
  left_join(tanner_prob, by = c("YEAR", "STATION_ID")) %>%
  mutate(
    snow_p   = replace_na(snow_p, 0),
    tanner_p = replace_na(tanner_p, 0)
  )

# Bhattacharyya coefficient by year
bhatt_year <- p_both %>%
  group_by(YEAR) %>%
  summarise(
    bhatt = sum(sqrt(snow_p * tanner_p))
  ) %>%
  rbind(data.frame(YEAR = 2020, bhatt = NA)) %>%
  mutate(bhatt = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
                                   TRUE ~ bhatt)) # filtering years with no chela measurments

write.csv(bhatt_year, "./Data/bhatt_matmalesnow.matfemTanner.csv")

ggplot(bhatt_year, aes(YEAR, bhatt))+
  geom_point()+
  geom_line()+
  theme_bw() +
  ggtitle("Overlap between mature male snow and mature female tanner")

# Overlap between snow mature females and tanner mature males ----
s.fem_t.male <- rbind(
  tanner_matmale_cpue %>%
    filter(CATEGORY == "Mature male") %>%
    dplyr::select(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE,
                  CATEGORY, DISTRICT, STRATUM, TOTAL_AREA, CPUE, CPUE_MT, CPUE_LBS),
  snow_fem_cpue %>%
    dplyr::select(!COUNT)) %>%
  filter(YEAR >= 1990) %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE, CATEGORY, SPECIES) %>%
  reframe(CPUE = sum(CPUE), .groups = "drop_last") %>%  # keep YEAR, STATION_ID, ... groups
  group_by(YEAR, SPECIES) %>%
  mutate(
    TOTAL_CPUE = sum(CPUE),
    P = ifelse(TOTAL_CPUE > 0, CPUE / TOTAL_CPUE, 0)
  ) %>%
  ungroup() %>%
  dplyr::select(-CPUE, -TOTAL_CPUE) %>%
  pivot_wider(
    names_from  = CATEGORY,
    values_from = P
  )


# Split by species
snow_prob <- s.fem_t.male %>%
  filter(SPECIES == "SNOW") %>%
  dplyr::select(YEAR, STATION_ID, snow_p = mature_female)

tanner_prob <- s.fem_t.male %>%
  filter(SPECIES == "TANNER") %>%
  dplyr::select(YEAR, STATION_ID, tanner_p = `Mature male`)

# Align on common YEAR × STATION_ID, fill missing probs with 0
years  <- sort(union(snow_prob$YEAR, tanner_prob$YEAR))
stns   <- sort(union(snow_prob$STATION_ID, tanner_prob$STATION_ID))

p_both <- expand_grid(YEAR = years, STATION_ID = stns) %>%
  left_join(snow_prob,   by = c("YEAR", "STATION_ID")) %>%
  left_join(tanner_prob, by = c("YEAR", "STATION_ID")) %>%
  mutate(
    snow_p   = replace_na(snow_p, 0),
    tanner_p = replace_na(tanner_p, 0)
  )

# Bhattacharyya coefficient by year
bhatt_year <- p_both %>%
  group_by(YEAR) %>%
  summarise(
    bhatt = sum(sqrt(snow_p * tanner_p))
  ) %>%
  rbind(data.frame(YEAR = 2020, bhatt = NA)) %>%
  mutate(bhatt = case_when(YEAR %in% c(2011, 2013, 2015, 2020) ~ NA,
                           TRUE ~ bhatt)) # filtering years with no chela measurments

write.csv(bhatt_year, "./Data/bhatt_matmaleTanner.matfemsnow.csv")

ggplot(bhatt_year, aes(YEAR, bhatt))+
  geom_point()+
  geom_line()+
  theme_bw() +
  ggtitle("Overlap between mature male tanner and mature female snow")

