source("./Scripts/make_helper_functions.R")

years_use = c(1988:2019, 2021:2025)
break_year = 2024
h = 0.1

hybrid_metrics <- run_niche_metrics(hybrid.mod, hybrid.mod.df2, years_use, break_year, h)

hybrid_plots <- plot_niche_results(hybrid_metrics, "hybrid", fig_dir = "./Figures")

snow_metrics <- run_niche_metrics(snow.mod, snow.mod.df2, years_use, break_year, h)

snow_plots <- plot_niche_results(snow_metrics, "snow", fig_dir = "./Figures")

tanner_metrics <- run_niche_metrics(tanner.mod, tanner.mod.df2, years_use, break_year, h)

tanner_plots <- plot_niche_results(tanner_metrics, "tanner", fig_dir = "./Figures")

# Plot combined depth ----
depth_pred_all <- dplyr::bind_rows(
  snow_metrics$depth_pred   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$depth_pred %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$depth_pred %>% dplyr::mutate(species = "Hybrid"))%>%
  dplyr::mutate(
    period = factor(period,
                    levels = c("Pre-heatwave", "Heatwave", "Post-heatwave")))

depth_mu_all <- dplyr::bind_rows(
  snow_metrics$depth_mu   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$depth_mu %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$depth_mu %>% dplyr::mutate(species = "Hybrid"))

p_depth_mu_tbl <- dplyr::bind_rows(
  hybrid_metrics$depth_mu_test %>% dplyr::mutate(species = "Hybrid"),
  snow_metrics$depth_mu_test   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$depth_mu_test %>% dplyr::mutate(species = "Tanner")) %>%
  dplyr::mutate(
    p_label = case_when(p_value < 0.05 & p_value>0.001 ~ "p<0.05*",
                        p_value < 0.001 ~ "p<0.001*",
                        TRUE ~  paste0("p=", signif(p_value, 2)))) %>%
  mutate(x = c(2003, 2003, 2003),
         y = c(1.45, 1.8, 2.25))


p_depth_resp <- ggplot(depth_pred_all,
                       aes(DEPTH_SCALED, est_resp, colour = period, group = YEAR)) +
  theme_bw() +
  geom_line()+
  labs(x = "Depth scaled", y = "CPUE") +
  scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                "Heatwave"     = "darkred",
                                "Post-heatwave"= "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free") +      # rows = species, 1st column
  theme(legend.position = "bottom", legend.direction = "horizontal")

p_depth_mu <- ggplot(depth_mu_all,
                     aes(YEAR, mu, colour = period)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Year", y = "Optimum (μ)") +
  scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                "Heatwave"     = "darkred",
                                "Post-heatwave"= "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free")+
  geom_text(
    data = p_depth_mu_tbl,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = p_label),
    hjust = 0, vjust = 0, size = 3.5)+
  guides(color = "none")

depth_fig <- p_depth_resp + p_depth_mu +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 11))

ggsave("./Figures/depth_combined_fig.png", width = 7, height = 6)

# Plot combined sed ----
sed_pred_all <- dplyr::bind_rows(
  snow_metrics$sed_pred   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$sed_pred %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$sed_pred %>% dplyr::mutate(species = "Hybrid")) %>%
  dplyr::mutate(
    period = factor(period,
                    levels = c("Pre-heatwave", "Heatwave", "Post-heatwave")))

sed_mu_all <- dplyr::bind_rows(
  snow_metrics$sed_mu   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$sed_mu %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$sed_mu %>% dplyr::mutate(species = "Hybrid"))

p_sed_mu_tbl <- dplyr::bind_rows(
  hybrid_metrics$sed_mu_test %>% dplyr::mutate(species = "Hybrid"),
  snow_metrics$sed_mu_test   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$sed_mu_test %>% dplyr::mutate(species = "Tanner")) %>%
  dplyr::mutate(
    p_label = case_when(p_value < 0.05 & p_value>0.001 ~ "p<0.05*",
                        p_value < 0.001 ~ "p<0.001*",
                        TRUE ~  paste0("p=", signif(p_value, 2)))) %>%
  mutate(x = c(2003, 2003, 2003),
         y = c(0.35, 0.9, 0.1))


p_sed_resp <- ggplot(sed_pred_all,
                       aes(SED_SCALED, est_resp, colour =  period, group = YEAR)) +
  theme_bw() +
  geom_line()+
  labs(x = "Sediment scaled", y = "CPUE") +
  scale_color_manual(values = c("cadetblue",
                                "darkred",
                                "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free") +      # rows = species, 1st column
  theme(legend.position = "bottom", legend.direction = "horizontal")

p_sed_mu <- ggplot(sed_mu_all,
                     aes(YEAR, mu, colour = period)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Year", y = "Optimum (μ)") +
  scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                "Heatwave"     = "darkred",
                                "Post-heatwave"= "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free")+
  geom_text(
    data = p_sed_mu_tbl,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = p_label),
    hjust = 0, vjust = 0, size = 3.5)+
  guides(color = "none")

sed_fig <- p_sed_resp + p_sed_mu +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 11))

ggsave("./Figures/sed_combined_fig.png", width = 7, height = 6)


# Plot combined ice ----
ice_pred_all <- dplyr::bind_rows(
  snow_metrics$ice_pred   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$ice_pred %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$ice_pred %>% dplyr::mutate(species = "Hybrid")) %>%
  dplyr::mutate(
    period = factor(period,
                    levels = c("Pre-heatwave", "Heatwave", "Post-heatwave")))

ice_slope_all <- dplyr::bind_rows(
  snow_metrics$ice_slopes   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$ice_slopes %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$ice_slopes %>% dplyr::mutate(species = "Hybrid"))

p_ice_slope_tbl <- dplyr::bind_rows(
  hybrid_metrics$ice_slope_test %>% dplyr::mutate(species = "Hybrid"),
  snow_metrics$ice_slope_test   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$ice_slope_test %>% dplyr::mutate(species = "Tanner")) %>%
  dplyr::mutate(
    p_label = case_when(p_value < 0.05 & p_value>0.001 ~ "p<0.05*",
                        p_value < 0.001 ~ "p<0.001*",
                        TRUE ~  paste0("p=", signif(p_value, 2)))) %>%
  mutate(x = c(2015, 2015, 2015),
         y = c(0.3, 0.18, -0.1))


p_ice_resp <- ggplot(ice_pred_all,
                     aes(ICE_SCALED, est_resp, colour =  period, group = YEAR)) +
  theme_bw() +
  geom_line()+
  labs(x = "Ice scaled", y = "CPUE") +
  scale_color_manual(values = c("cadetblue",
                                "darkred",
                                "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free") +      # rows = species, 1st column
  theme(legend.position = "bottom", legend.direction = "horizontal")

p_ice_slope <- ggplot(ice_slope_all,
                   aes(YEAR, slope, colour = period)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Year", y = "Slope") +
  scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                "Heatwave"     = "darkred",
                                "Post-heatwave"= "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free")+
  geom_text(
    data = p_ice_slope_tbl,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = p_label),
    hjust = 0, vjust = 0, size = 3.5)+
  guides(color = "none")

ice_fig <- p_ice_resp + p_ice_slope +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 11))

ggsave("./Figures/ice_combined_fig.png", width = 7, height = 6)


# Plot combined btemp ----
btemp_pred_all <- dplyr::bind_rows(
  snow_metrics$btemp_pred   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$btemp_pred %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$btemp_pred %>% dplyr::mutate(species = "Hybrid")) %>%
  dplyr::mutate(period = factor(period,
                    levels = c("Pre-heatwave", "Heatwave", "Post-heatwave")))

btemp_slope_all <- dplyr::bind_rows(
  snow_metrics$btemp_slopes   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$btemp_slopes %>% dplyr::mutate(species = "Tanner"),
  hybrid_metrics$btemp_slopes %>% dplyr::mutate(species = "Hybrid"))

p_btemp_slope_tbl <- dplyr::bind_rows(
  hybrid_metrics$btemp_slope_test %>% dplyr::mutate(species = "Hybrid"),
  snow_metrics$btemp_slope_test   %>% dplyr::mutate(species = "Snow"),
  tanner_metrics$btemp_slope_test %>% dplyr::mutate(species = "Tanner")) %>%
  dplyr::mutate(
    p_label = case_when(p_value < 0.05 & p_value>0.001 ~ "p<0.05*",
                        p_value < 0.001 ~ "p<0.001*",
                        TRUE ~  paste0("p=", signif(p_value, 2)))) %>%
  mutate(x = c(1995, 1995, 1995),
         y = c(-0.7, -0.8, -0.48))


p_btemp_resp <- ggplot(btemp_pred_all,
                     aes(BTEMP_SCALED, est_resp, colour =  period, group = YEAR)) +
  theme_bw() +
  geom_line()+
  labs(x = "Bottom temperature scaled", y = "CPUE") +
  scale_color_manual(values = c("cadetblue",
                                "darkred",
                                "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free") +      # rows = species, 1st column
  theme(legend.position = "bottom", legend.direction = "horizontal")

p_btemp_slope <- ggplot(btemp_slope_all,
                      aes(YEAR, slope, colour = period)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(x = "Year", y = "Slope") +
  scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                "Heatwave"     = "darkred",
                                "Post-heatwave"= "gold"),
                     name = "") +
  facet_grid(rows = vars(species), scales = "free")+
  geom_text(
    data = p_btemp_slope_tbl,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = p_label),
    hjust = 0, vjust = 0, size = 3.5)+
  guides(color = "none")

btemp_fig <- p_btemp_resp + p_btemp_slope +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center",
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 11))

ggsave("./Figures/btemp_combined_fig.png", width = 7, height = 6)
