source("./Scripts/load_libs_params.R")


hybrid_data <- readRDS("./Data/hybrid_specimen.rda")


spec2 <- hybrid_data$specimen %>%
          filter(SEX == 1 & SIZE >=79| SEX == 2 & SIZE >=65) # check this!!!!

hybrid_data$specimen <- spec2


# Calculate per-station CPUE by 1mm size bin
hybrid_cpue <- crabpack::calc_cpue(crab_data = hybrid_data,
                                   species = "HYBRID",
                                   region = "EBS",
                                   size_min = NULL,
                                   size_max = NULL)


# transform model dat to sf
mod.dat <- hybrid_cpue %>%
              group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
              reframe(CPUE = sum(CPUE)) %>%
              st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
              st_transform(., crs = "+proj=utm +zone=2")

# Add sediment
f <- "Y:/KOD_Survey/EBS Shelf/Spatial crab/Environmental layers/EBS_phi_1km.grd"
sed <- raster(f)
sed.df <- as.data.frame(sed, xy = TRUE)  # xy = coordinates of cell centers
tt <- crs(sed)
sed.sf <- sed.df %>%
  st_as_sf(coords = c("x", "y"),
           crs = tt) %>%
  st_transform(crs = "+proj=utm +zone=2")

mod.dat <- st_join(mod.dat, sed.sf, join = st_nearest_feature) %>%
            rename(sed = X3.pred)

# Add ice by year
ice <- read.csv("./Output/spatial_ice_means_1980-2025.csv") %>%
  group_by(year, latitude, longitude) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(crs = st_crs(mod.dat))   # ensure same CRS

mod.dat$ice <- NA_real_

for (y in sort(unique(mod.dat$YEAR))) {
  idx_cpue <- which(mod.dat$YEAR == y)
  if (!length(idx_cpue)) next
  
  ice_y <- ice[ice$year == y, ]
  if (nrow(ice_y) == 0) next
  
  # nearest ice cell for that year
  nn_idx <- st_nearest_feature(mod.dat[idx_cpue, ], ice_y)
  
  # assign ice value
  mod.dat$ice[idx_cpue] <- ice_y$value[nn_idx]
}

# Transform to data frame
mod.df <- mod.dat %>%
            cbind(st_coordinates(.)) %>%
            as.data.frame(.) %>%
            mutate(LONGITUDE = X/1000,
                   LATITUDE = Y/1000) %>%
            dplyr::select(!c(X, Y, geometry)) %>%
            rename(SED = sed, ICE = ice)

# Add depth and bottom temp
mod.df2 <- hybrid_data$haul %>%
            dplyr::select(YEAR, STATION_ID, BOTTOM_DEPTH, GEAR_TEMPERATURE) %>%
            right_join(mod.df, .) %>%
            filter(YEAR >=1988) %>%
            rename(DEPTH = BOTTOM_DEPTH, BTEMP = GEAR_TEMPERATURE) %>%
            mutate(YEAR_F = as.factor(YEAR),
                   YEAR_SCALED = scale(YEAR),
                   ICE_SCALED = scale(ICE),
                   DEPTH_SCALED = scale(DEPTH),
                   BTEMP_SCALED = scale(BTEMP),
                   SED_SCALED = scale(SED),
                DEPTH_SCALED = as.numeric(DEPTH_SCALED),
                BTEMP_SCALED = as.numeric(BTEMP_SCALED),
                ICE_SCALED   = as.numeric(ICE_SCALED),
                SED_SCALED   = as.numeric(SED_SCALED)) %>%
            filter(is.finite(DEPTH_SCALED),
                   is.finite(BTEMP_SCALED),
                   is.finite(ICE_SCALED),
                   is.finite(SED_SCALED))
          

# Plot
ggplot(mod.df2, aes(LONGITUDE, LATITUDE, fill = ICE)) +
  geom_tile(width = 45, height = 45) +
  theme_bw() +
  facet_wrap(~YEAR) +
  scale_fill_viridis_c()

ggplot(mod.df2, aes(LONGITUDE, LATITUDE, fill = SED)) +
  geom_tile(width = 45, height = 45) +
  theme_bw() +
  facet_wrap(~YEAR) +
  scale_fill_viridis_c()

ggplot(mod.df2, aes(LONGITUDE, LATITUDE, fill = DEPTH)) +
  geom_tile(width = 45, height = 45) +
  theme_bw() +
  facet_wrap(~YEAR) +
  scale_fill_viridis_c()

ggplot(mod.df2, aes(LONGITUDE, LATITUDE, fill = log(CPUE+10))) +
  geom_tile(width = 45, height = 45) +
  theme_bw() +
  facet_wrap(~YEAR) +
  scale_fill_viridis_c()


# FIT MODELS WITH COVARIATEES ----
# Build mesh
mesh <- make_mesh(mod.df2, c("LONGITUDE","LATITUDE"), n_knots = 90, type = "kmeans")

# Fit model
mod <- sdmTMB(
        CPUE ~ 0 + s(DEPTH_SCALED, k=4) + s(BTEMP_SCALED, k=4) + s(ICE_SCALED, k=4) + s(SED_SCALED, k=4),
        time_varying = ~ 0 + DEPTH_SCALED + ICE_SCALED + BTEMP_SCALED + SED_SCALED,
        time_varying_type = "ar1",           # or "ar1"
        mesh = mesh,
        extra_time = c(2020),
        family = tweedie(link = "log"),
        time = "YEAR",
        spatial = "on",
        data = mod.df2,
        silent = FALSE
      )

# Diagnostics
sanity(mod)

# Plot effects
cov_means <- mod.df2 %>%
  summarise(
    DEPTH_SCALED = mean(DEPTH_SCALED, na.rm = TRUE),
    BTEMP_SCALED = mean(BTEMP_SCALED, na.rm = TRUE),
    ICE_SCALED   = mean(ICE_SCALED,   na.rm = TRUE),
    SED_SCALED   = mean(SED_SCALED,   na.rm = TRUE)
  )

years_use <- sort(unique(mod.df2$YEAR))


make_tv_effect <- function(mod, data, focal, n,
                           years,
                           cov_means) {
  
  rng <- range(data[[focal]], na.rm = TRUE)
  x   <- seq(rng[1], rng[2], length.out = n)
  
  nd <- expand.grid(
    focal_val = x,
    YEAR      = years
  ) |>
    as_tibble() |>
    mutate(
      DEPTH_SCALED = cov_means$DEPTH_SCALED,
      BTEMP_SCALED = cov_means$BTEMP_SCALED,
      ICE_SCALED   = cov_means$ICE_SCALED,
      SED_SCALED   = cov_means$SED_SCALED
    ) |>
    mutate(!!focal := focal_val, .keep = "unused")
  
  # predictions on link scale
  pred_raw <- predict(
    mod,
    newdata = nd,
    re_form = ~ 0,
    type    = "link",   # IMPORTANT: link scale
    se_fit  = TRUE
  )
  
  # keep only prediction columns
  pred <- pred_raw |>
    dplyr::select(est, est_se)
  
  out <- bind_cols(nd, pred)
  
  # back-transform Tweedie/log link to response scale
  out <- out |>
    mutate(
      est_resp = exp(est))
  
  out
}

# call function (no duplicate-column messages now)
pred_depth <- make_tv_effect(
  mod,
  data      = mod.df2,
  focal     = "DEPTH_SCALED",
  n         = 100,
  years     = years_use,
  cov_means = cov_means
)
options(bitmapType = "cairo") 
ggplot() +
  geom_line(pred_depth, mapping = aes(DEPTH_SCALED, est_resp, color = YEAR, group = YEAR), size = 1)+
  theme_bw()+
  scale_color_viridis_c()+
  geom_ribbon(pred_depth, mapping = aes(DEPTH_SCALED, ymin = est_resp - est_se,
                                        ymax = est_resp + est_se, fill = YEAR), alpha = 0.5)


# FIT MODELS WITHOUT COVARIATES ----
# Build mesh
mesh <- make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 90, type = "kmeans")

# Fit model
mod <- sdmTMB(CPUE ~ 0 + YEAR_F, #the 0 is there so there is a factor predictor for each time slice
                spatial = "on",
                spatiotemporal = "iid",
                mesh = mesh,
                family = delta_gamma(type = "poisson-link"),
                time = "YEAR",
                extra_time = c(2020),
                anisotropy = TRUE,
                data = mod.df2)

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
