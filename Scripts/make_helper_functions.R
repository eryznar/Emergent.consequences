make_tv_effect_quad <- function(mod, data,
                                focal,
                                n,
                                years,
                                cov_means) {
  
  # variables actually used in the fitted model matrix
  mm_names <- colnames(stats::model.matrix(mod))  # robust to interactions, etc.
  
  # use full observed range for the focal covariate
  rng <- range(data[[focal]], na.rm = TRUE)
  x   <- seq(rng[1], rng[2], length.out = n)
  
  nd <- expand.grid(
    focal_val = x,
    YEAR      = years
  ) |>
    dplyr::as_tibble() |>
    dplyr::mutate(
      DEPTH_SCALED = cov_means$DEPTH_SCALED,
      SED_SCALED   = cov_means$SED_SCALED,
      BTEMP_SCALED = cov_means$BTEMP_SCALED,
      ICE_SCALED   = cov_means$ICE_SCALED
    ) |>
    # overwrite focal covariate with grid values
    dplyr::mutate(!!focal := focal_val, .keep = "unused")
  
  ## --- quadratic terms: only add if the model matrix has them ---
  if ("DEPTH_SCALED2" %in% mm_names) {
    nd$DEPTH_SCALED2 <- nd$DEPTH_SCALED^2
  }
  if ("SED_SCALED2" %in% mm_names) {
    nd$SED_SCALED2 <- nd$SED_SCALED^2
  }
  
  ## --- family-specific prediction handling ---
  fam <- mod$family
  is_delta <- inherits(fam, "sdmTMB_delta_family")
  
  if (is_delta) {
    # delta_* family: combined mean on response scale
    pred_raw <- predict(
      mod,
      newdata = nd,
      re_form = ~0,
      type    = "response",
      se_fit  = FALSE
    )
    est_lin  <- pred_raw$est_lin
    est_resp <- pred_raw$est
  } else {
    # e.g. Tweedie(log)
    pred_raw <- predict(
      mod,
      newdata = nd,
      re_form = ~0,
      type    = "link",
      se_fit  = FALSE
    )
    est_lin  <- pred_raw$est
    est_resp <- exp(pred_raw$est)
  }
  
  dplyr::bind_cols(
    nd,
    tibble::tibble(
      est      = est_lin,
      est_resp = est_resp
    )
  ) |>
    dplyr::mutate(
      period = dplyr::case_when(
        YEAR < 2018         ~ "Pre-heatwave",
        YEAR %in% 2018:2019 ~ "Heatwave",
        TRUE                ~ "Post-heatwave"
      )
    )
}

get_mu <- function(x, w) {
  w <- pmax(w, 0)
  w <- w / sum(w)
  mu <- sum(x * w)
  tibble(mu = mu)
}

niche_mu_test <- function(pred_dat,
                          focal = c("DEPTH_SCALED", "SED_SCALED",
                                    "ICE_SCALED", "BTEMP_SCALED"),
                          break_year = 2024,
                          nperm = 999) {
  focal <- match.arg(focal)
  
  mu_df <- pred_dat %>%
    dplyr::group_by(YEAR) %>%
    dplyr::reframe(get_mu(.data[[focal]], est_resp)) %>%
    dplyr::mutate(period = if_else(YEAR >= break_year, "post", "pre"))
  
  obs_stats <- mu_df %>%
    dplyr::group_by(period) %>%
    dplyr::summarise(mu_mean = mean(mu), .groups = "drop")
  
  obs_T_mu <- obs_stats$mu_mean[obs_stats$period == "post"] -
    obs_stats$mu_mean[obs_stats$period == "pre"]
  
  perm_T_mu <- numeric(nperm)
  
  for (i in seq_len(nperm)) {
    perm_period <- sample(mu_df$period)
    perm_stats <- mu_df %>%
      dplyr::mutate(perm_period = perm_period) %>%
      dplyr::group_by(perm_period) %>%
      dplyr::summarise(mu_mean = mean(mu), .groups = "drop")
    
    perm_T_mu[i] <- perm_stats$mu_mean[perm_stats$perm_period == "post"] -
      perm_stats$mu_mean[perm_stats$perm_period == "pre"]
  }
  
  p_mu <- mean(abs(perm_T_mu) >= abs(obs_T_mu))
  
  tibble(
    focal    = focal,
    obs_T_mu = obs_T_mu,
    p_mu     = p_mu
  )
}

niche_slope_test <- function(slope_dat,
                             slope_col = "slope",
                             break_year = 2024,
                             nperm = 999) {
  
  library(dplyr)
  slope_sym <- rlang::sym(slope_col)
  
  dat <- slope_dat %>%
    ungroup() %>%
    filter(!is.na(!!slope_sym)) %>%        # remove NA slopes
    mutate(period = if_else(YEAR >= break_year, "post", "pre"))
  
  # observed stats
  obs_stats <- dat %>%
    group_by(period) %>%
    summarise(mean_slope = mean(!!slope_sym), .groups = "drop")
  
  obs_T <- obs_stats$mean_slope[obs_stats$period == "post"] -
    obs_stats$mean_slope[obs_stats$period == "pre"]
  
  nperm <- 999
  perm_T <- numeric(nperm)
  
  for (i in seq_len(nperm)) {
    perm_period <- sample(dat$period)
    perm_stats <- dat %>%
      mutate(perm_period = perm_period) %>%
      group_by(perm_period) %>%
      summarise(mean_slope = mean(!!slope_sym), .groups = "drop")
    perm_T[i] <- perm_stats$mean_slope[perm_stats$perm_period == "post"] -
      perm_stats$mean_slope[perm_stats$perm_period == "pre"]
  }
  
  p_val <- mean(abs(perm_T) >= abs(obs_T))
  
  tibble(obs_T = obs_T, p_value = p_val)
}

test_mu_break <- function(mu_df, break_year = 2018) {
  dat <- mu_df %>%
    dplyr::mutate(
      period_bin = ifelse(YEAR >= break_year, 1L, 0L)
    )
  
  fit <- lm(mu ~ period_bin, data = dat)
  
  broom::tidy(fit) %>%
    dplyr::filter(term == "period_bin") %>%
    dplyr::rename(
      estimate_diff = estimate,   # mean(post) - mean(pre)
      p_value       = p.value
    )
}

test_slope_break <- function(slope_df, break_year) {

  dat <- slope_df %>%
    dplyr::mutate(
      period_bin = ifelse(YEAR >= break_year, 1L, 0L)
    )
  
  fit <- lm(slope ~ period_bin, data = dat)
  
  broom::tidy(fit) %>%
    dplyr::filter(term == "period_bin") %>%
    dplyr::rename(
      estimate_diff = estimate,   # mean(post) - mean(pre)
      p_value       = p.value
    )
}

get_sdmTMB_tv_coef <- function(mod, tv_term) {
  # tidy() with effects = "fixed" and "ran_pars" doesn't give tv coefs,
  # but effects = "ran_vals" does, and includes time-varying betas
  
  tv <- tidy(mod, effects = "ran_vals") %>%
    dplyr::filter(grepl(tv_term, .data$term))
  
 
  tv %>%
    dplyr::transmute(
      YEAR  = as.integer(sub(".*:(\\d+)$", "\\1", term)),
      slope = estimate
    ) %>%
    dplyr::arrange(YEAR)
}

run_niche_metrics <- function(mod,
                              dat,
                              years_use = c(1988:2019, 2021:2025),
                              break_year = 2024,
                              h = 0.1) {
  
  # 1. Covariate means (for non-focal covariates)
  cov_means <- dat %>%
    dplyr::summarise(
      DEPTH_SCALED = mean(DEPTH_SCALED, na.rm = TRUE),
      SED_SCALED   = mean(SED_SCALED,   na.rm = TRUE),
      BTEMP_SCALED = mean(BTEMP_SCALED, na.rm = TRUE),
      ICE_SCALED   = mean(ICE_SCALED,   na.rm = TRUE)
    )
  
  add_period <- function(df) {
    df %>%
      dplyr::mutate(
        period = dplyr::case_when(
          YEAR < 2018 ~ "Pre-heatwave",
          YEAR %in% 2018:2019 ~ "Heatwave",
          TRUE ~ "Post-heatwave"
        )
      )
  }
  
  ## ---------- DEPTH (μ) ----------
  message("Depth")
  pred_depth <- make_tv_effect_quad(
    mod       = mod,
    data      = dat,
    focal     = "DEPTH_SCALED",
    n         = 100,
    years     = years_use,
    cov_means = cov_means
  )
  
  depth_mu <- pred_depth %>%
    dplyr::group_by(YEAR) %>%
    dplyr::reframe(get_mu(DEPTH_SCALED, est_resp)) %>%
    add_period()
  
  depth_mu_test <- test_mu_break(depth_mu, break_year = break_year)
  
  ## ---------- SEDIMENT (μ) ----------
  message("Sediment")
  pred_sed <- make_tv_effect_quad(
    mod       = mod,
    data      = dat,
    focal     = "SED_SCALED",
    n         = 100,
    years     = years_use,
    cov_means = cov_means
  )
  
  sed_mu <- pred_sed %>%
    dplyr::group_by(YEAR) %>%
    dplyr::reframe(get_mu(SED_SCALED, est_resp)) %>%
    add_period()
  
  sed_mu_test <- test_mu_break(sed_mu, break_year = break_year)
  
  ## ---------- ICE (slopes) ----------
  message("Ice")
  pred_ice <- make_tv_effect_quad(
    mod       = mod,
    data      = dat,
    focal     = "ICE_SCALED",
    n         = 100,
    years     = years_use,
    cov_means = cov_means
  ) %>%
    add_period()
  
  ice_slopes <- pred_ice %>%
    dplyr::group_by(YEAR, period) %>%
    dplyr::summarise(
      ice_minus = ICE_SCALED[which.min(abs(ICE_SCALED + h))],
      ice_zero  = ICE_SCALED[which.min(abs(ICE_SCALED))],
      ice_plus  = ICE_SCALED[which.min(abs(ICE_SCALED - h))],
      eta_minus = est[which.min(abs(ICE_SCALED + h))],
      eta_zero  = est[which.min(abs(ICE_SCALED))],
      eta_plus  = est[which.min(abs(ICE_SCALED - h))],
      slope     = (eta_plus - eta_minus) / (ice_plus - ice_minus),
      .groups = "drop"
    ) %>%
    dplyr::full_join(
      expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)),
                  period = unique(.$period)),
      by = c("YEAR", "period")
    )
  
  ice_slope_test <- test_slope_break(
    ice_slopes, break_year = break_year
  )
  
  ## ---------- BTEMP (slopes) ----------
  message("Bottom temperature")
  pred_bt <- make_tv_effect_quad(
    mod       = mod,
    data      = dat,
    focal     = "BTEMP_SCALED",
    n         = 100,
    years     = years_use,
    cov_means = cov_means
  ) %>%
    add_period()
  
  btemp_slopes <- pred_bt %>%
    dplyr::group_by(YEAR, period) %>%
    dplyr::summarise(
      btemp_minus = BTEMP_SCALED[which.min(abs(BTEMP_SCALED + h))],
      btemp_zero  = BTEMP_SCALED[which.min(abs(BTEMP_SCALED))],
      btemp_plus  = BTEMP_SCALED[which.min(abs(BTEMP_SCALED - h))],
      eta_minus   = est[which.min(abs(BTEMP_SCALED + h))],
      eta_zero    = est[which.min(abs(BTEMP_SCALED))],
      eta_plus    = est[which.min(abs(BTEMP_SCALED - h))],
      slope       = (eta_plus - eta_minus) / (btemp_plus - btemp_minus),
      .groups = "drop"
    ) %>%
    dplyr::full_join(
      expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)),
                  period = unique(.$period)),
      by = c("YEAR", "period")
    )
  
  btemp_slope_test <- test_slope_break(
    btemp_slopes, break_year = break_year
  )
  
  # Return everything needed for plotting / further analysis
  list(
    cov_means        = cov_means,
    depth_pred       = pred_depth,
    depth_mu         = depth_mu,
    depth_mu_test    = depth_mu_test,
    sed_pred         = pred_sed,
    sed_mu           = sed_mu,
    sed_mu_test      = sed_mu_test,
    ice_pred         = pred_ice,
    ice_slopes       = ice_slopes,
    ice_slope_test   = ice_slope_test,
    btemp_pred       = pred_bt,
    btemp_slopes     = btemp_slopes,
    btemp_slope_test = btemp_slope_test
  )
}


plot_niche_results <- function(metrics,
                               species_name,
                               fig_dir = "./Figures") {
  
  species_title <- tools::toTitleCase(species_name)
  
  # helper: ordered period factor for plots
  order_period <- function(df) {
    df %>%
      dplyr::mutate(
        period = factor(
          period,
          levels = c("Pre-heatwave", "Heatwave", "Post-heatwave")
        )
      )
  }
  
  ## ---------- Depth μ and curves ----------
  
  depth_mu   <- order_period(metrics$depth_mu)
  pred_depth <- order_period(metrics$depth_pred)
  
  p_depth_mu <- ggplot(depth_mu,
                       aes(x = YEAR, y = mu, colour = period)) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      x = "Year",
      y = "Optimum (μ)",
      title = paste0(species_title, " depth optimum (μ)")
    ) +
    theme(legend.position = "bottom",
          axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  p_depth_curves <- ggplot() +
    geom_line(pred_depth,
              mapping = aes(DEPTH_SCALED, est_resp,
                            colour = period, group = YEAR),
              size = 1) +
    theme_bw() +
    labs(
      x = "Depth scaled",
      y = "CPUE",
      title = paste0(species_title, " depth response")
    ) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    theme(axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ## ---------- Sediment μ and curves ----------
  
  sed_mu   <- order_period(metrics$sed_mu)
  pred_sed <- order_period(metrics$sed_pred)
  
  p_sed_mu <- ggplot(sed_mu,
                     aes(x = YEAR, y = mu, colour = period)) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      x = "Year",
      y = "Optimum (μ)",
      title = paste0(species_title, " sediment optimum (μ)")
    ) +
    theme(legend.position = "bottom",
          axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  p_sed_curves <- ggplot() +
    geom_line(pred_sed,
              mapping = aes(SED_SCALED, est_resp,
                            colour = period, group = YEAR),
              size = 1) +
    theme_bw() +
    labs(
      x = "Sediment scaled",
      y = "CPUE",
      title = paste0(species_title, " sediment response")
    ) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    theme(axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ## ---------- ICE slopes and curves ----------
  
  ice_slopes <- order_period(metrics$ice_slopes)
  pred_ice   <- order_period(metrics$ice_pred)
  
  p_ice_slope <- ggplot(ice_slopes,
                        aes(x = YEAR, y = slope, colour = period)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      x = "Year",
      y = "Slope",
      title = paste0(species_title, " ice slope")
    ) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    theme(axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  p_ice_curves <- ggplot() +
    geom_line(pred_ice,
              mapping = aes(ICE_SCALED, est_resp,
                            colour = period, group = YEAR),
              size = 1) +
    theme_bw() +
    labs(
      x = "Ice scaled",
      y = "CPUE",
      title = paste0(species_title, " ice response")
    ) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    theme(axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  ## ---------- BTEMP slopes and curves ----------
  
  btemp_slopes <- order_period(metrics$btemp_slopes)
  pred_bt      <- order_period(metrics$btemp_pred)
  
  p_bt_slope <- ggplot(btemp_slopes,
                       aes(x = YEAR, y = slope, colour = period)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      x = "Year",
      y = "Slope",
      title = paste0(species_title, " bottom temperature slope")
    ) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    theme(axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  p_bt_curves <- ggplot() +
    geom_line(pred_bt,
              mapping = aes(BTEMP_SCALED, est_resp,
                            colour = period, group = YEAR),
              size = 1) +
    theme_bw() +
    labs(
      x = "Bottom temperature scaled",
      y = "CPUE",
      title = paste0(species_title, " bottom temperature response")
    ) +
    scale_color_manual(values = c("Pre-heatwave" = "cadetblue",
                                  "Heatwave"     = "darkred",
                                  "Post-heatwave"= "gold"),
                       name = "") +
    theme(axis.text  = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  # return plots (no ggsave, nothing written to disk)
  list(
    bt_curves  = p_bt_curves,
    bt_slope   = p_bt_slope,
    ice_curves = p_ice_curves,
    ice_slope  = p_ice_slope,
    sed_curves = p_sed_curves,
    sed_mu     = p_sed_mu,
    depth_curves = p_depth_curves,
    depth_mu     = p_depth_mu
  )
}


# DIAGOSTIC FUNCTION ----
plot.resids <- function(model, species, model_name){
  resids <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, return_DHARMa = TRUE)
  
  dat <- cbind(model$data, DHARMa_resid = resids$scaledResiduals)
  
  rr_yr  <- dat %>%
    group_by(YEAR) %>%
    arrange(DHARMa_resid, .by_group = TRUE) %>%
    mutate(
      n = n(),
      expected = ppoints(n),         # uniform quantiles
      observed = sort(DHARMa_resid)  # sort residuals for QQ
    ) %>%
    ungroup() %>%
    mutate(model = model_name)
  
  #  QQ plot with ggplot2
  ggplot()+
    theme_bw()+
    geom_point(rr_yr, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("observed")+
    xlab("expected")+
    facet_wrap(~YEAR)+
    scale_x_continuous(breaks = c(0, 0.5, 1))+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    ggtitle(paste0(species, " ", model_name)) -> by_yr
  
  
  dat2 <- dat %>%
    group_by(STATION_ID) %>%
    mutate(LONGITUDE = mean(LONGITUDE), LATITUDE = mean(LATITUDE)) %>%
    ungroup()
  
  ggplot(dat2, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
    geom_point(shape = 21, size = 1.75, stroke = NA)+
    facet_wrap(~YEAR)+
    scale_fill_gradient2(midpoint = 0.5)+
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text = element_text(size = 10)) +
    ggtitle(ggtitle(paste0(species, " ", model_name))) -> by_yr_sp
  
  return(list(by_yr = by_yr, by_yr_sp = by_yr_sp,
              rr_yr = rr_yr))
}
