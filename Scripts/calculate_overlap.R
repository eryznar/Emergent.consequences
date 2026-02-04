snow_matmale_cpue <- read.csv("./Data/snow_maturemale_cpue.csv")
tanner_matmale_cpue <- read.csv("./Data/tanner_maturemale_cpue.csv")


# females
snow_spec <- readRDS("./Data/snow_survey_specimenEBS.rda")
tanner_spec <- readRDS("./Data/tanner_survey_specimenEBS.rda")



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

# large and small males instead of mature/immature
snow_lgsm_cpue <- crabpack::calc_cpue(crab_data = snow_spec,
                                     region = "EBS",
                                     species = "SNOW",
                                     years = 1980:2025,
                                     crab_category = c("large_male", "small_male"))

write.csv(rbind(snow_lgsm_cpue %>% dplyr::select(!COUNT), read.csv("./Data/snow_cpue.csv") %>% dplyr::select(!X)),
          "./Data/snow_cpue.csv")

snow_lgsm_bioabund <- crabpack::calc_bioabund(crab_data = snow_spec,
                                      region = "EBS",
                                      species = "SNOW",
                                      years = 1980:2025,
                                      crab_category = c("large_male", "small_male", "mature_female", "immature_female"))

write.csv(snow_lgsm_bioabund, "./Data/snow_bioabund.csv")


tanner_lgsm_cpue <- crabpack::calc_cpue(crab_data = tanner_spec,
                                      region = "EBS",
                                      species = "TANNER",
                                      years = 1980:2025,
                                      crab_category = c("large_male", "small_male"))

write.csv(rbind(tanner_lgsm_cpue %>% dplyr::select(!COUNT), read.csv("./Data/tanner_cpue.csv") %>% dplyr::select(!X)),
          "./Data/tanner_cpue.csv")

tanner_lgsm_bioabund <- crabpack::calc_bioabund(crab_data = tanner_spec,
                                        region = "EBS",
                                        species = "TANNER",
                                        years = 1980:2025,
                                        crab_category = c("large_male", "small_male", "mature_female", "immature_female"))

write.csv(tanner_lgsm_bioabund, "./Data/tanner_bioabund.csv")


# Overlap between snow mature males and tanner mature females ----
s.male_t.fem <- rbind(
  snow_matmale_cpue %>%
    filter(CATEGORY == "Mature male") %>%
    dplyr::select(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE,
                  CATEGORY, DISTRICT, STRATUM, TOTAL_AREA, CPUE, CPUE_MT, CPUE_LBS),
  tanner_fem_cpue %>%
    dplyr::select(!COUNT)) %>%
  filter(YEAR >= 1980) %>%
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
bhatt_year1 <- p_both %>%
  group_by(YEAR) %>%
  summarise(
    bhatt = sum(sqrt(snow_p * tanner_p))
  ) %>%
  rbind(data.frame(YEAR = 2020, bhatt = NA)) %>%
  mutate(bhatt = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
                           TRUE ~ bhatt)) %>%
  mutate(comparison = "snowmatmale-tannermatfem") # filtering years with no chela measurments

#write.csv(bhatt_year, "./Data/bhatt_matmalesnow.matfemTanner.csv")

ggplot(bhatt_year1, aes(YEAR, bhatt))+
  geom_point()+
  geom_line()+
  theme_bw() +
  ggtitle("Overlap between mature male snow and mature female tanner")

# Overlap between snow large males and tanner mature females ----
s.male_t.fem <- rbind(
  snow_lgsm_cpue %>%
    filter(CATEGORY == "large_male") %>%
    dplyr::select(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE,
                  CATEGORY, DISTRICT, STRATUM, TOTAL_AREA, CPUE, CPUE_MT, CPUE_LBS),
  tanner_fem_cpue %>%
    dplyr::select(!COUNT)) %>%
  filter(YEAR >= 1980) %>%
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
  dplyr::select(YEAR, STATION_ID, snow_p = large_male)

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
bhatt_year2 <- p_both %>%
  group_by(YEAR) %>%
  summarise(
    bhatt = sum(sqrt(snow_p * tanner_p))
  ) %>%
  rbind(data.frame(YEAR = 2020, bhatt = NA)) %>%
  mutate(comparison = "snowlgmale-tannermatfem") # filtering years with no chela measurments

#write.csv(bhatt_year, "./Data/bhatt_matmalesnow.matfemTanner.csv")
bhatt_compare1 <- rbind(bhatt_year1, bhatt_year2)

ggplot(bhatt_compare1, aes(YEAR, bhatt, color = comparison))+
  geom_point()+
  annotate("text", x = 1995, y = 0.45, label = "r = 0.68") +
  geom_line()+
  theme_bw()


bhatt_year1 %>%
  na.omit() -> pp

bhatt_year2 %>%
  na.omit() %>%
  filter(YEAR %in% pp$YEAR) -> tt

cor(pp$bhatt, tt$bhatt)

write.csv(bhatt_compare1, "./Data/bhatt_malesnow.femTanner.csv")

# Overlap between snow mature females and tanner mature males ----
s.fem_t.male <- rbind(
  tanner_matmale_cpue %>%
    filter(CATEGORY == "Mature male") %>%
    dplyr::select(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE,
                  CATEGORY, DISTRICT, STRATUM, TOTAL_AREA, CPUE, CPUE_MT, CPUE_LBS),
  snow_fem_cpue %>%
    dplyr::select(!COUNT)) %>%
  filter(YEAR >= 1980) %>%
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
bhatt_year3 <- p_both %>%
  group_by(YEAR) %>%
  summarise(
    bhatt = sum(sqrt(snow_p * tanner_p))
  ) %>%
  rbind(data.frame(YEAR = 2020, bhatt = NA)) %>%
  mutate(bhatt = case_when(YEAR %in% c(2011, 2013, 2015, 2020) ~ NA,
                           TRUE ~ bhatt)) %>% # filtering years with no chela measurments
  mutate(comparison = "tannermatmale-snowmatfem")

#write.csv(bhatt_year3, "./Data/bhatt_matmaleTanner.matfemsnow.csv")

ggplot(bhatt_year3, aes(YEAR, bhatt))+
  geom_point()+
  geom_line()+
  theme_bw() +
  ggtitle("Overlap between mature male tanner and mature female snow")

# Overlap between Tanner large males and snow mature females ----
t.male_s.fem <- rbind(
  tanner_lgsm_cpue %>%
    filter(CATEGORY == "large_male") %>%
    dplyr::select(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE,
                  CATEGORY, DISTRICT, STRATUM, TOTAL_AREA, CPUE, CPUE_MT, CPUE_LBS),
  snow_fem_cpue %>%
    dplyr::select(!COUNT)) %>%
  filter(YEAR >= 1980) %>%
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
tanner_prob <- t.male_s.fem %>%
  filter(SPECIES == "TANNER") %>%
  dplyr::select(YEAR, STATION_ID, snow_p = large_male)

snow_prob <- t.male_s.fem %>%
  filter(SPECIES == "SNOW") %>%
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
bhatt_year4 <- p_both %>%
  group_by(YEAR) %>%
  summarise(
    bhatt = sum(sqrt(snow_p * tanner_p))
  ) %>%
  rbind(data.frame(YEAR = 2020, bhatt = NA)) %>%
  mutate(comparison = "tannerlgmale-snowmatfem") # filtering years with no chela measurments

#write.csv(bhatt_year, "./Data/bhatt_matmalesnow.matfemTanner.csv")
bhatt_compare2 <- rbind(bhatt_year3, bhatt_year4)

ggplot(bhatt_compare2, aes(YEAR, bhatt, color = comparison)) +
  geom_point() +
  geom_line() +
  annotate("text", x = 1995, y = 0.3, label = "r = 0.89") +
  theme_bw()


bhatt_year3 %>%
  na.omit() -> pp

bhatt_year4 %>%
  na.omit() %>%
  filter(YEAR %in% pp$YEAR) -> tt

cor(pp$bhatt, tt$bhatt)

write.csv(bhatt_compare2, "./Data/bhatt_maleTanner.femsnow.csv")

write.csv(rbind(bhatt_compare1, bhatt_compare2), "./Data/bhatt_snowTanner.csv")

ggplot(rbind(bhatt_compare1, bhatt_compare2), aes(YEAR, bhatt, color  = comparison))+
  geom_line()+
  geom_point()
