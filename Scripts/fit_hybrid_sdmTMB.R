source("./Scripts/load_libs_params.R")
source("./Scripts/make_helper_functions.R")



hybrid_data <- readRDS("./Data/hybrid_specimen.rda")


spec2 <- hybrid_data$specimen %>%
          filter(SIZE >=48 & SIZE <=103) # central 90% of abundance for hybrids in 2025 across sexes

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
                   is.finite(SED_SCALED)) %>%
            mutate(
              DEPTH_SCALED2 = DEPTH_SCALED^2,
              SED_SCALED2   = SED_SCALED^2
            )
          

# Plot
ggplot(mod.df2 %>% filter(YEAR %in% 2016:2025), aes(LONGITUDE, LATITUDE, fill = ICE)) +
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

ggplot(mod.df2%>% filter(YEAR %in% 2016:2025), aes(LONGITUDE, LATITUDE, fill = log(CPUE+10))) +
  geom_tile(width = 45, height = 45) +
  theme_bw() +
  facet_wrap(~YEAR) +
  scale_fill_viridis_c()

ggplot(mod.df2 %>% filter(YEAR %in% 2016:2025), aes(LONGITUDE, LATITUDE, fill = BTEMP)) +
  geom_tile(width = 45, height = 45) +
  theme_bw() +
  facet_wrap(~YEAR) +
  scale_fill_viridis_c()


cor_dat <- mod.df2 %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(ICE_SCALED, DEPTH_SCALED, BTEMP_SCALED, SED_SCALED)

# compute correlation matrix
cor_mat <- cor(cor_dat, use = "pairwise.complete.obs")

# plot
corrplot::corrplot(
  cor_mat,
  method      = "color",
  type        = "upper",
  tl.col      = "black",
  tl.srt      = 45,
  addCoef.col = "black",   # show correlation values
  number.cex  = 0.7        # adjust as needed
)

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

# Diagnostics
sanity(mod)

# covariate means (for non-focal covariates)
cov_means <- mod.df2 %>%
  summarise(
    DEPTH_SCALED = mean(DEPTH_SCALED, na.rm = TRUE),
    SED_SCALED   = mean(SED_SCALED,   na.rm = TRUE),
    BTEMP_SCALED = mean(BTEMP_SCALED, na.rm = TRUE),
    ICE_SCALED   = mean(ICE_SCALED,   na.rm = TRUE)
  )

# Focal years
years_use <- c(1990:2019, 2021:2025)

# Depth
pred_depth <- make_tv_effect_quad(
  mod,
  data      = mod.df2,
  focal     = "DEPTH_SCALED",
  n         = 100,
  years     = years_use,
  cov_means = cov_means
)

depth_mu_sd <- pred_depth %>%
  group_by(YEAR) %>%
  reframe(get_mu_sigma(DEPTH_SCALED, est_resp)) %>%
  mutate(period = case_when(YEAR %in% 2024:2025 ~ "2024:2025",
                            TRUE ~ "<2024"))

niche_n_mu_sigma_test(pred_depth,
                      focal = c("DEPTH_SCALED"),
                      break_year = 2024,
                      nperm = 999)

ggplot(depth_mu_sd, aes(x = sigma, y = mu, color = period)) +
  scale_color_manual(values = c("cadetblue", "gold"))+
  geom_point(size = 2) +
  geom_text(aes(label = YEAR),
            nudge_y = 0.03, size = 4)+
  theme_bw()+
  ylab("Optimum (μ)")+
  xlab("Breadth (σ)")+
  #ggtitle("Depth")+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("./Figures/hybrid_depth_musigma.png", width = 7, height = 5)
  

ggplot() +
  geom_line(pred_depth, mapping = aes(DEPTH_SCALED, est_resp, color = period, group = YEAR), size = 1)+
  theme_bw()+
  ylab("Hybrid CPUE")+
  #ggtitle("Depth")+
  xlab("Depth scaled")+
  scale_color_manual(values = c("cadetblue", "gold"), name = "")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/hybrid_depth_predcurves.png", width = 7, height = 5)


# Sediment
pred_sed <- make_tv_effect_quad(
  mod,
  data      = mod.df2,
  focal     = "SED_SCALED",
  n         = 100,
  years     = years_use,
  cov_means = cov_means
)


niche_n_mu_sigma_test(pred_sed,
                      focal = c("SED_SCALED"),
                      break_year = 2024,
                      nperm = 999)

sed_mu_sd <- pred_sed %>%
  group_by(YEAR) %>%
  reframe(get_mu_sigma(SED_SCALED, est_resp)) %>%
  mutate(period = case_when(YEAR %in% 2024:2025 ~ "2024:2025",
                            TRUE ~ "<2024"))

ggplot(sed_mu_sd, aes(x = sigma, y = mu, color = period)) +
  scale_color_manual(values = c("cadetblue", "gold"))+
  geom_point() +
  geom_text(aes(label = YEAR),
            nudge_y = 0.03, size = 4)+
  theme_bw()+
  ylab("Optimum (μ)")+
  xlab("Breadth (σ)")+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("./Figures/hybrid_sed_musigma.png", width = 7, height = 5)

ggplot() +
  geom_line(pred_sed, mapping = aes(SED_SCALED, est_resp, color = period, group = YEAR), size = 1)+
  theme_bw()+
  ylab("Hybrid CPUE")+
  #ggtitle("Depth")+
  xlab("Sediment scaled")+
  scale_color_manual(values = c("cadetblue", "gold"), name = "")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/hybrid_sed_predcurves.png", width = 7, height = 5)

# ICE
pred_ice <- make_tv_effect_quad(
  mod,
  data      = mod.df2,
  focal     = "ICE_SCALED",
  n         = 100,
  years     = years_use,
  cov_means = cov_means
)

# choose a small step for numerical derivative
h <- 0.1

# for each year, compute derivative of log-mean at ICE_SCALED = 0
ice_slopes <- pred_ice %>%
  group_by(YEAR, period) %>%
  summarise(
    # nearest points to -h, 0, +h
    ice_minus = ICE_SCALED[which.min(abs(ICE_SCALED + h))],
    ice_zero  = ICE_SCALED[which.min(abs(ICE_SCALED))],
    ice_plus  = ICE_SCALED[which.min(abs(ICE_SCALED - h))],
    eta_minus = est[which.min(abs(ICE_SCALED + h))],
    eta_zero  = est[which.min(abs(ICE_SCALED))],
    eta_plus  = est[which.min(abs(ICE_SCALED - h))],
    slope = (eta_plus - eta_minus) / (ice_plus - ice_minus)
  ) %>%
  full_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), period = unique(.$period)))

niche_slope_test(ice_slopes, slope_col = "slope", break_year = 2024, nperm = 999)

ggplot(ice_slopes, aes(x = YEAR, y = slope, colour = period)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  #ggtitle("Ice")+
  ylab("slope")+
  xlab("Year")+
  theme_bw()+
  scale_color_manual(values = c("cadetblue", "gold"), name = "")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/hybrid_ice_slope.png", width = 7, height = 5)

ggplot() +
  geom_line(pred_ice, mapping = aes(ICE_SCALED, est_resp, color = period, group = YEAR), size = 1)+
  theme_bw()+
  ylab("Hybrid CPUE")+
  #ggtitle("Depth")+
  xlab("Ice scaled")+
  scale_color_manual(values = c("cadetblue", "gold"), name = "")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/hybrid_ice_predcurves.png", width = 7, height = 5)


# BTEMP
pred_bt <- make_tv_effect_quad(
  mod,
  data      = mod.df2,
  focal     = "BTEMP_SCALED",
  n         = 100,
  years     = years_use,
  cov_means = cov_means
)

# 1. Compute slope of log-mean vs BTEMP_SCALED at BTEMP ≈ 0 for each year
btemp_slopes <- pred_bt %>%
  group_by(YEAR, period) %>%
  summarise(
    # nearest points to -h, 0, +h
    btemp_minus = BTEMP_SCALED[which.min(abs(BTEMP_SCALED + h))],
    btemp_zero  = BTEMP_SCALED[which.min(abs(BTEMP_SCALED))],
    btemp_plus  = BTEMP_SCALED[which.min(abs(BTEMP_SCALED - h))],
    eta_minus   = est[which.min(abs(BTEMP_SCALED + h))],
    eta_zero    = est[which.min(abs(BTEMP_SCALED))],
    eta_plus    = est[which.min(abs(BTEMP_SCALED - h))],
    slope  = (eta_plus - eta_minus) / (btemp_plus - btemp_minus)
  ) %>%
  full_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), period = unique(.$period)))

niche_slope_test(btemp_slopes, slope_col = "slope", break_year = 2024, nperm = 999)

ggplot(btemp_slopes, aes(x = YEAR, y = slope, colour = period)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2) +
  #ggtitle("Ice")+
  ylab("slope")+
  xlab("Year")+
  theme_bw()+
  scale_color_manual(values = c("cadetblue", "gold"), name = "")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/hybrid_bt_slope.png", width = 7, height = 5)

ggplot() +
  geom_line(pred_bt, mapping = aes(BTEMP_SCALED, est_resp, color = period, group = YEAR), size = 1)+
  theme_bw()+
  ylab("Hybrid CPUE")+
  #ggtitle("Depth")+
  xlab("Bottom temperature scaled")+
  scale_color_manual(values = c("cadetblue", "gold"), name = "")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/hybrid_bt_predcurves.png", width = 7, height = 5)

