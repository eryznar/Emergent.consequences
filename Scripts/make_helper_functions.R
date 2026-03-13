
# Helper functions
make_tv_effect_quad <- function(mod, data,
                                focal,          # "DEPTH_SCALED" or "SED_SCALED"
                                n,
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
      SED_SCALED   = cov_means$SED_SCALED,
      BTEMP_SCALED = cov_means$BTEMP_SCALED,
      ICE_SCALED   = cov_means$ICE_SCALED
    ) |>
    mutate(!!focal := focal_val, .keep = "unused")
  
  # ALWAYS create both quadratic cols expected by the model
  nd$DEPTH_SCALED2 <- nd$DEPTH_SCALED^2
  nd$SED_SCALED2   <- nd$SED_SCALED^2
  
  pred_raw <- predict(
    mod,
    newdata = nd,
    re_form = ~ 0,
    type    = "link",
    se_fit  = TRUE
  )
  
  out <- bind_cols(
    nd,
    dplyr::select(pred_raw, est, est_se)
  ) |>
    mutate(est_resp = exp(est),
           period = case_when((YEAR <2018) ~ "Pre-heatwave",
                                       (YEAR %in% 2018:2019) ~ "Heatwave",
                                       TRUE ~ "Post-heatwave"))
  
  out
}


get_mu_sigma <- function(x, w) {
  w <- pmax(w, 0)
  w <- w / sum(w)
  mu <- sum(x * w)
  sigma <- sqrt(sum((x - mu)^2 * w))
  tibble(mu = mu, sigma = sigma)   # <- data frame, not vector
}

niche_n_mu_sigma_test <- function(pred_dat,
                                  focal = c("DEPTH_SCALED", "SED_SCALED"),
                                  break_year = 2024,
                                  nperm = 999) {
  focal <- match.arg(focal)
  
  # summarise to mu, sigma per year
  mu_sd <- pred_dat %>%
    group_by(YEAR) %>%
    reframe(get_mu_sigma(.data[[focal]], est_resp)) %>%
    mutate(period = if_else(YEAR >= break_year, "post", "pre"))
  
  # observed stats
  obs_stats <- mu_sd %>%
    group_by(period) %>%
    summarise(
      mu_mean    = mean(mu),
      sigma_mean = mean(sigma),
      .groups = "drop"
    )
  
  obs_T_mu    <- obs_stats$mu_mean[obs_stats$period == "post"] -
    obs_stats$mu_mean[obs_stats$period == "pre"]
  obs_T_sigma <- obs_stats$sigma_mean[obs_stats$period == "post"] -
    obs_stats$sigma_mean[obs_stats$period == "pre"]
  
  # permutation
  perm_T_mu    <- numeric(nperm)
  perm_T_sigma <- numeric(nperm)
  
  for (i in seq_len(nperm)) {
    perm_period <- sample(mu_sd$period)
    
    perm_stats <- mu_sd %>%
      mutate(perm_period = perm_period) %>%
      group_by(perm_period) %>%
      summarise(
        mu_mean    = mean(mu),
        sigma_mean = mean(sigma),
        .groups = "drop"
      )
    
    perm_T_mu[i] <- perm_stats$mu_mean[perm_stats$perm_period == "post"] -
      perm_stats$mu_mean[perm_stats$perm_period == "pre"]
    perm_T_sigma[i] <- perm_stats$sigma_mean[perm_stats$perm_period == "post"] -
      perm_stats$sigma_mean[perm_stats$perm_period == "pre"]
  }
  
  p_mu    <- mean(abs(perm_T_mu)    >= abs(obs_T_mu))
  p_sigma <- mean(abs(perm_T_sigma) >= abs(obs_T_sigma))
  
  tibble(
    focal       = focal,
    obs_T_mu    = obs_T_mu,
    p_mu        = p_mu,
    obs_T_sigma = obs_T_sigma,
    p_sigma     = p_sigma
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
