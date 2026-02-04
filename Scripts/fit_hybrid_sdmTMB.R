source("./Scripts/load_libs_params.R")

# 3) Set channel
channel <- "API"

# 4) Set species
species <- "HYBRID"

# # 5) Pull specimen data
# hybrid_data <- crabpack::get_specimen_data(species = "HYBRID",
#                                            region = "EBS", # can also include NBS
#                                            channel = channel)
# 
# saveRDS(hybrid_data, "./Data/hybrid_specimen.rda")

hybrid_data <- readRDS("./Data/hybrid_specimen.rda")


spec2 <- hybrid_data$specimen %>%
          filter(SEX == 1 & SIZE >=79| SEX == 2 & SIZE >=65) # check this!!!!

hybrid_data$specimen <- spec2


# 6) Calculate per-station CPUE by 1mm size bin
hybrid_cpue <- crabpack::calc_cpue(crab_data = hybrid_data,
                                   species = "HYBRID",
                                   region = "EBS",
                                   size_min = NULL,
                                   size_max = NULL)

db_abund <-  crabpack::calc_bioabund(crab_data = hybrid_data,
                                 species = "HYBRID",
                                 region = "EBS",
                                 size_min = NULL,
                                 size_max = NULL) %>%
                mutate(ABUNDANCE = ABUNDANCE/1e6,
                       ABUNDANCE_CI = ABUNDANCE_CI/1e6) %>%
                filter(YEAR >=1980)

mod.dat <- hybrid_cpue %>%
              group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
              reframe(CPUE = sum(CPUE)) %>%
              st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
              st_transform(., crs = "+proj=utm +zone=2") %>%
              cbind(st_coordinates(.)) %>%
              as.data.frame(.) %>%
              mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
                     LONGITUDE = X/1000,
                     YEAR_F = as.factor(YEAR)) %>%
              dplyr::select(!c(X, Y, geometry)) %>%
              filter(YEAR >=1980)


# FIT MODELS ----
# Build mesh
mesh <- make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 90, type = "kmeans")

# Fit model
# mod.1 <- sdmTMB(CPUE ~ 0 + YEAR_F, #the 0 is there so there is a factor predictor for each time slice
#               spatial = "on",
#               spatiotemporal = "iid",
#               mesh = mesh,
#               family = tweedie(link = "log"),
#               time = "YEAR",
#               extra_time = c(2020),
#               anisotropy = TRUE,
#               data = mod.dat)
# 
# saveRDS(mod.1, "./Models/hybrid_sdmTMB_tw_90.rda")

mod.2 <- sdmTMB(CPUE ~ 0 + YEAR_F, #the 0 is there so there is a factor predictor for each time slice
                spatial = "on",
                spatiotemporal = "iid",
                mesh = mesh,
                family = delta_gamma(type = "poisson-link"),
                time = "YEAR",
                extra_time = c(2020),
                anisotropy = TRUE,
                data = mod.dat)

saveRDS(mod.2, "./Models/hybrid_sdmTMB_dg_90.rda")

# Diagnostics
sanity(mod.2)

# Predict 
years <- unique(mod.dat$YEAR)
pred.grid <- read.csv("./Data/ebs_coarse_grid.csv") %>%
                dplyr::select(area_km2, X, Y) %>%
                replicate_df(., "YEAR", years) %>%
                rename(LONGITUDE = X, LATITUDE = Y) %>%
                mutate(YEAR_F = as.factor(YEAR),
                       area_nmi2 = area_km2 * 0.29155335)


mod.2 <- readRDS("./Models/hybrid_sdmTMB_dg_90.rda")
preds <- predict(mod.2, newdata= pred.grid, return_tmb_object = T)

# Get index
abund <- get_index(preds, bias_correct = TRUE, area = pred.grid$area_nmi2) 

mb_abund <- abund %>%
              mutate(ABUNDANCE = est/1e6,
                             lwr = lwr/1e6,
                             upr = upr/1e6)


# Plot index
ggplot()+
  geom_ribbon(db_abund, mapping=aes(YEAR, ymin = ABUNDANCE - ABUNDANCE_CI,
                                    ymax = ABUNDANCE + ABUNDANCE_CI, fill = as.factor(1)), alpha = 0.25)+
  geom_ribbon(mb_abund, mapping=aes(YEAR, ymin = ABUNDANCE - lwr,
                                    ymax = ABUNDANCE + upr, fill = as.factor(2)), alpha = 0.25)+
  geom_line(db_abund, mapping=aes(YEAR, ABUNDANCE, color = as.factor(1)), linewidth = 1.25)+
  geom_line(mb_abund, mapping=aes(YEAR, ABUNDANCE, color = as.factor(2)), linewidth = 1.25)+
  geom_line(tv_abund, mapping=aes(YEAR, ABUNDANCE, color = as.factor(3)), linewidth = 1.25)+
  theme_bw()+
  scale_color_manual(values = c("cadetblue", "darkgoldenrod", "salmon"), labels = c("design-based", "sdmTMB", "tinyVAST"), name = "")+
  scale_fill_manual(values = c("cadetblue", "darkgoldenrod", "salmon"), labels = c("design-based", "sdmTMB", "tinyVAST"), name = "")
  
# Plot spatial
spat.preds <- predict(mod.2, newdata= pred.grid, type = "response")


ggplot(spat.preds %>% filter(YEAR %in% c(2019:2025)), aes(LONGITUDE, LATITUDE, fill = log(est)))+
  geom_tile(width = 27, height = 27)+
  theme_bw()+
  facet_wrap(~YEAR, nrow = 3)+
  scale_fill_viridis_c(option = 'mako')


# ADD COVARIATES ----
# Ice
# Ice as sf in UTM, retaining year
ice <- read.csv("./Output/spatial_ice_means_1980-2025.csv") %>%
  group_by(year, latitude, longitude) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(crs = "+proj=utm +zone=2")

# CPUE as sf in same CRS (assumes LONGITUDE/LATITUDE are WGS84)
mod_sf <- st_as_sf(mod.dat,
                   coords = c("LONGITUDE", "LATITUDE"),
                   crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(crs = st_crs(ice))

mod_sf$ice <- NA_real_

for (y in sort(unique(mod_sf$YEAR))) {
  # CPUE points for this year
  idx_cpue <- which(mod_sf$YEAR == y)
  if (!length(idx_cpue)) next
  
  # Ice points for this year
  ice_y <- ice[ice$year == y, ]
  if (nrow(ice_y) == 0) next
  
  # Nearest ice point index for each CPUE point (within this year)
  nn_idx <- st_nearest_feature(mod_sf[idx_cpue, ], ice_y)
  
  # Assign ice value
  mod_sf$ice[idx_cpue] <- ice_y$value[nn_idx]
}

mod.dat$ICE <- mod_sf$ice[match(seq_len(nrow(mod.dat)), as.integer(st_drop_geometry(mod_sf$rowid %||% 1:nrow(mod_sf))))]

ggplot(mod.dat, aes(LONGITUDE, LATITUDE, fill = ICE)) +
  geom_tile(width = 45, height = 45) +
  theme_bw() +
  facet_wrap(~YEAR) +
  scale_fill_viridis_c()

# Snow cpue
snow_cpue <- crabpack::calc_cpue(crab_data = readRDS("./Data/snow_specimen.rda"),
                                 species = "SNOW",
                                 region = "EBS")

# Tanner cpue
tanner_cpue <- crabpack::calc_cpue(crab_data = readRDS("./Data/tanner_specimen.rda"),
                                   species = "TANNER",
                                   region = "EBS")

# Add snow and tanner cpue to mod.dat
mod.dat2 <- mod.dat %>%
  as.data.frame() %>%
  dplyr::select(!c(X, Y, geometry)) %>%
  rename(HYBRID_CPUE = CPUE) %>%
  inner_join(.,
             snow_cpue %>%
               dplyr::select(YEAR, STATION_ID, CPUE) %>%
               rename(SNOW_CPUE = CPUE) %>%
               filter(YEAR >= 1988),
             by = c("YEAR", "STATION_ID")) %>%
  inner_join(.,
             tanner_cpue %>%
               dplyr::select(YEAR, STATION_ID, CPUE) %>%
               rename(TANNER_CPUE = CPUE) %>%
               filter(YEAR >= 1988),
             by = c("YEAR", "STATION_ID"))

