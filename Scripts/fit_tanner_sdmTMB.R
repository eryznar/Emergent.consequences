source("./Scripts/load_libs_params.R")
source("./Scripts/make_helper_functions.R")

tanner_data <- readRDS("./Data/tanner_specimen.rda")


spec2 <- tanner_data$specimen 
tanner_data$specimen <- spec2

# Calculate per-station CPUE by 1mm size bin
tanner_cpue <- crabpack::calc_cpue(crab_data = tanner_data,
                                 species = "TANNER",
                                 region = "EBS",
                                 size_min = 50,
                                 size_max = NULL,
                                shell_condition = c("soft_molting", "new_hardshell"))


# transform model dat to sf
mod.dat <- tanner_cpue %>%
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
  filter(name == "Mar-Apr ice") %>%
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
mod.df2 <- tanner_data$haul %>%
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
         is.finite(SED_SCALED)) %>%
  mutate(
    DEPTH_SCALED2 = DEPTH_SCALED^2,
    SED_SCALED2   = SED_SCALED^2
  )


# # Plot
# ggplot(mod.df2 %>% filter(YEAR %in% 2016:2025), aes(LONGITUDE, LATITUDE, fill = ICE)) +
#   geom_tile(width = 45, height = 45) +
#   theme_bw() +
#   facet_wrap(~YEAR) +
#   scale_fill_viridis_c()
# 
# ggplot(mod.df2, aes(LONGITUDE, LATITUDE, fill = SED)) +
#   geom_tile(width = 45, height = 45) +
#   theme_bw() +
#   facet_wrap(~YEAR) +
#   scale_fill_viridis_c()
# 
# ggplot(mod.df2, aes(LONGITUDE, LATITUDE, fill = DEPTH)) +
#   geom_tile(width = 45, height = 45) +
#   theme_bw() +
#   facet_wrap(~YEAR) +
#   scale_fill_viridis_c()
# 
# ggplot(mod.df2%>% filter(YEAR %in% 2016:2025), aes(LONGITUDE, LATITUDE, fill = log(CPUE+10))) +
#   geom_tile(width = 45, height = 45) +
#   theme_bw() +
#   facet_wrap(~YEAR) +
#   scale_fill_viridis_c()
# 
# ggplot(mod.df2 %>% filter(YEAR %in% 2016:2025), aes(LONGITUDE, LATITUDE, fill = BTEMP)) +
#   geom_tile(width = 45, height = 45) +
#   theme_bw() +
#   facet_wrap(~YEAR) +
#   scale_fill_viridis_c()
# 
# 
# cor_dat <- mod.df2 %>%
#   dplyr::select(where(is.numeric)) %>%
#   dplyr::select(ICE_SCALED, DEPTH_SCALED, BTEMP_SCALED, SED_SCALED)
# 
# # compute correlation matrix
# cor_mat <- cor(cor_dat, use = "pairwise.complete.obs")
# 
# # plot
# corrplot::corrplot(
#   cor_mat,
#   method      = "color",
#   type        = "upper",
#   tl.col      = "black",
#   tl.srt      = 45,
#   addCoef.col = "black",   # show correlation values
#   number.cex  = 0.7        # adjust as needed
# )

# FIT MODELS WITH COVARIATEES ----
# Build mesh
mesh <- make_mesh(mod.df2, c("LONGITUDE","LATITUDE"), n_knots = 90, type = "kmeans")

# Fit model
mod <- sdmTMB(
  CPUE ~ 1 +
    DEPTH_SCALED + DEPTH_SCALED2 +      # unimodal depth
    SED_SCALED   + SED_SCALED2 +        # unimodal sediment
    BTEMP_SCALED +                      # baseline linear temp
    ICE_SCALED,                         # baseline linear ice
  time_varying = ~ 1 +
    DEPTH_SCALED + DEPTH_SCALED2 +      # depth niche shifts through time
    SED_SCALED   + SED_SCALED2 +        # sediment niche shifts
    BTEMP_SCALED +                      # temp effect can change by year
    ICE_SCALED,                         # ice effect can change by year
  time_varying_type = "ar1",
  mesh      = mesh,
  extra_time = 2020,
  family    = tweedie(link = "log"),
  time      = "YEAR",
  spatial   = "on",
  data      = mod.df2,
  silent    = FALSE
)
