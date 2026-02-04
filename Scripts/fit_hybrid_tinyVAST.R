# PURPOSE: fit tinyVAST spatiotemporal models

# LOAD LIBS/PARAMS ----
source("./Scripts/load_libs_params.R")

# LOAD/PROCESS DATA ----
hybrid_data <- readRDS("./Data/hybrid_specimen.rda")

spec2 <- hybrid_data$specimen %>%
  filter(SEX == 1 & SIZE >=79| SEX == 2 & SIZE >=65) # check this!!!!

hybrid_data$specimen <- spec2


# Calculate per-station CPUE and pop-level bioabund
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


# Prepare cpue data for modeling
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
  filter(YEAR >=1980) %>%
  droplevels()

# FIT MODELS ----
# Build mesh
mesh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 90, type = "kmeans")

tv_mesh <- mesh$mesh

mod.dat$var <- as.character("CPUE")

# Fit model
fit_tv <- tinyVAST(
  formula        = CPUE ~ 0 + YEAR_F,
  data           = mod.dat,
  family         = delta_gamma(type = "poisson-link"),
  spatial_domain = tv_mesh,
  space_columns  = c("LONGITUDE", "LATITUDE"),
  space_term     = "CPUE <-> CPUE, sd_space",
  spacetime_term = "CPUE <-> CPUE, 0, sd_spacetime",  # if you want IID S×T
  time_column    = "YEAR",
  
  control = tinyVASTcontrol(
    getsd   = TRUE,
    profile = "alpha_j"
  )
)


saveRDS(fit_tv, "./Models/hybrid_tinyVAST_db_90.rda")

# Prediction grid
years <- unique(mod.dat$YEAR)
pred.grid <- read.csv("./Data/ebs_coarse_grid.csv") %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "YEAR", years) %>%
  rename(LONGITUDE = X, LATITUDE = Y) %>%
  mutate(YEAR_F = as.factor(YEAR),
         area_nmi2 = area_km2 * 0.29155335)

years <- sort(unique(pred.grid$YEAR))

idx_by_year <- purrr::map_dfr(
  years,
  function(yy) {
    pred.grid2 <- pred.grid %>% dplyr::filter(YEAR == yy)
    
    ii <- integrate_output(
      object       = fit_tv,
      newdata      = pred.grid2,
      area         = pred.grid2$area_nmi2,
      type         = rep(1, nrow(pred.grid2)),
      getsd        = FALSE,
      bias.correct = FALSE
    )
    # integrate_output returns a numeric vector:
    # c(Estimate, Std.Error, Est.bias.correct, Std.bias.correct)
    # With getsd=FALSE, bias.correct=FALSE, only the first element is non‑NA.[web:105]
    
    est_val <- ii[1]  # first element = plug‑in estimate
    
    tibble::tibble(
      YEAR = yy,
      est  = est_val,
      se   = NA_real_
    )
  }
)


tv_abund <- idx_by_year %>%
  mutate(ABUNDANCE = est/1e6)



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
mod_sf <-  hybrid_cpue %>%
              group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
              reframe(CPUE = sum(CPUE)) %>%
              st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
              st_transform(., crs = "+proj=utm +zone=2") %>%
              cbind(st_coordinates(.)) %>%
              dplyr::select(!c(X, Y, geometry)) %>%
              filter(YEAR >=1980) %>%
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
             by = c("YEAR", "STATION_ID")) %>%
  mutate(HYBRID_CPUE = scale(log(HYBRID_CPUE + 10))[,1],
         SNOW_CPUE = scale(log(SNOW_CPUE +10))[,1],
         TANNER_CPUE = scale(log(TANNER_CPUE + 10))[,1]) # log transform and scale

# Make into long format where hybrid and tanner cpue are the variables
long_dat <- mod.dat2 %>%
              pivot_longer(., cols = c("HYBRID_CPUE", "TANNER_CPUE"), names_to = "resp_name", values_to = "resp_value")

# FIT DSEM MODELS ----
# Fit model
fit_tv <- tinyVAST(
  formula        = CPUE ~ 0 + YEAR_F,
  data           = mod.dat,
  family         = delta_gamma(type = "poisson-link"),
  spatial_domain = tv_mesh,
  space_columns  = c("LONGITUDE", "LATITUDE"),
  space_term     = "CPUE <-> CPUE, sd_space",
  spacetime_term = "CPUE <-> CPUE, 0, sd_spacetime",  # if you want IID S×T
  time_column    = "YEAR",
  
  control = tinyVASTcontrol(
    getsd   = TRUE,
    profile = "alpha_j"
  )
)
