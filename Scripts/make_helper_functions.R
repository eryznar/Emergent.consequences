make_tv_effect_quad <- function(mod, data,
                                focal,
                                n,
                                years,
                                cov_means) {
  
  # variables in fixed effects (for detecting *_SCALED2)
  mod_vars <- all.vars(formula(mod))
  
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
  
  # quadratic terms (fixed effects only; no time-varying quadratics)
  if ("DEPTH_SCALED2" %in% mod_vars)
    nd$DEPTH_SCALED2 <- nd$DEPTH_SCALED^2
  if ("SED_SCALED2" %in% mod_vars)
    nd$SED_SCALED2   <- nd$SED_SCALED^2
  
  pred_raw <- predict(
    mod,
    newdata = nd,
    re_form = ~0,          # fixed effects (incl. quadratics) only
    type    = "link",
    se_fit  = FALSE
  )
  
  dplyr::bind_cols(
    nd,
    dplyr::select(pred_raw, est)
  ) |>
    dplyr::mutate(
      est_resp = exp(est),
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
