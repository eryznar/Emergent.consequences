
# Hybrid abundance
hybrid_data <- readRDS("./Data/hybrid_specimen.rda")


spec2 <- hybrid_data$specimen %>%
  filter(SEX == 1 & SIZE >=79| SEX == 2 & SIZE >=65) # check this!!!!

hybrid_data$specimen <- spec2


hybrid_abund <- crabpack::calc_bioabund(crab_data = hybrid_data,
                                   species = "HYBRID",
                                   region = "EBS",
                                   size_min = NULL,
                                   size_max = NULL) %>%
  dplyr::select(YEAR, ABUNDANCE) %>%
  rename(HYBRID_ABUND = ABUNDANCE)


# Snow abundance
snow_abund <- read.csv("./Data/snow_bioabund.csv") %>%
  filter(YEAR >= 1980) %>%
  group_by(YEAR) %>%
  reframe(SNOW_ABUND = sum(ABUNDANCE)/1e6)

# Tanner abundance
tanner_abund <- read.csv("./Data/tanner_bioabund.csv") %>%
  filter(YEAR >= 1980) %>%
  group_by(YEAR) %>%
  reframe(TANNER_ABUND = sum(ABUNDANCE)/1e6)

# Ice
ice <- read.csv("./Output/ice_means_1980-2025.csv") %>%
  group_by(year) %>%
  reframe(ICE = mean(value)) %>%
  filter(year>= 1980 & year != 2020)

# Overlap
st_overlap <- read.csv("./Data/bhatt_snowTanner.csv") %>%
  filter(comparison == "snowlgmale-tannermatfem", YEAR !=2020)

# Join
lagdat <- cbind(hybrid_abund, SNOW_ABUND = snow_abund$SNOW_ABUND, TANNER_ABUND = tanner_abund$TANNER_ABUND, ICE = ice$ICE,
                BHATT = st_overlap$bhatt) %>%
          full_join(., data.frame(YEAR = 2020)) %>%
          arrange(., YEAR) %>%
          mutate(SNOW_ABUND_avg2 = zoo::rollmean(SNOW_ABUND, 2, fill = NA, align = "right"),
                 SNOW_ABUND_avg3 = zoo::rollmean(SNOW_ABUND, 3, fill = NA, align = "right"),
                 TANNER_ABUND_avg2 = zoo::rollmean(TANNER_ABUND, 2, fill = NA, align = "right"),
                 TANNER_ABUND_avg3 = zoo::rollmean(TANNER_ABUND, 3, fill = NA, align = "right"),
                 ICE_avg2 = zoo::rollmean(ICE, 2, fill = NA, align = "right"),
                 ICE_avg3 = zoo::rollmean(ICE, 3, fill = NA, align = "right"),
                 BHATT_avg2 = zoo::rollmean(BHATT, 2, fill = NA, align = "right"),
                 BHATT_avg3 = zoo::rollmean(BHATT, 3, fill = NA, align = "right"))

# Run ccf
vars <- c("SNOW_ABUND", "TANNER_ABUND", "ICE", "ICE_avg3", "ICE_avg2",
          "SNOW_ABUND_avg2", "SNOW_ABUND_avg3",
          "TANNER_ABUND_avg2", "TANNER_ABUND_avg3", "BHATT", "BHATT_avg2", "BHATT_avg3")

# candidate lags (in years)
lag <- 6   # adjust as needed

ccf_list <- list()

for (v in vars) {
  x <- lagdat$HYBRID_ABUND
  y <- lagdat[[v]]
  
  # remove rows with any NA in the pair
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  
  cc <- ccf(x2, y2, lag.max = lag,
            plot = FALSE, na.action = na.omit)
  
  ccf_df <- data.frame(
    var  = v,
    lag  = cc$lag,
    ccf  = cc$acf
  )
  
  ccf_list[[v]] <- ccf_df
}

ccf_all <-  bind_rows(ccf_list) %>%
                mutate(
                  smooth = case_when(
                    grepl("avg2", var) ~ "2-year",
                    grepl("avg3", var) ~ "3-year",
                    TRUE               ~ "none"
                  ),
                  short_var = case_when(
                    grepl("avg2", var) ~ gsub("_avg2", "", var),
                    grepl("avg3", var) ~ gsub("_avg3", "", var),
                    TRUE               ~ var
                  )
                )

# Plot
ggplot(ccf_all, aes(lag, ccf, fill = smooth))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~short_var)+
  scale_x_continuous(breaks = seq(-lag, lag, by = 1))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())


# Select top lag/smooths
best_lags <- ccf_all %>%
              group_by(short_var) %>%
              slice_max(abs(ccf), n=2)

# Add top lags to model df
lagdat2 <- lagdat %>%
           mutate(ICE_lag6 = lag(ICE, 6),
                  ICE_avg2lag4 = lag(ICE_avg2, 4),
                  SNOW_ABUND_lag4 = lag(SNOW_ABUND, 4)) %>%
           dplyr::select(YEAR, HYBRID_ABUND,
                         BHATT_avg2, BHATT_avg3,
                         ICE_lag6, ICE_avg2lag4,
                         SNOW_ABUND_lag4, SNOW_ABUND,
                         TANNER_ABUND_avg3, TANNER_ABUND_avg2)

# CV function
k_folds <- 5

cv_rmse_gamm <- function(fml, data, k_folds = 5) {
  data <- data %>% arrange(YEAR)
  n    <- nrow(data)
  folds <- cut(seq_len(n), breaks = k_folds, labels = FALSE)
  
  errs <- numeric(k_folds)
  
  for (k in seq_len(k_folds)) {
    test_idx  <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    train_dat <- data[train_idx, , drop = FALSE]
    test_dat  <- data[test_idx,  , drop = FALSE]
    
    # fit GAMM with AR(1) on YEAR
    fit_k <- gamm(
      formula     = fml,
      data        = train_dat,
      family      = nb(),
      correlation = corAR1(),
      method      = "REML"
    )
    
    # use the GAM component for prediction
    pred <- predict(fit_k$gam, newdata = test_dat, type = "response")
    
    errs[k] <- sqrt(mean((test_dat$HYBRID_ABUND - pred)^2, na.rm = TRUE))
  }
  
  mean(errs)
}

## 2. Parameter grid (unchanged) ----

snow.pars   <- c(NA, names(lagdat2)[grep("SNOW_ABUND",   names(lagdat2))])
tanner.pars <- c(NA, names(lagdat2)[grep("TANNER_ABUND", names(lagdat2))])
bhatt.pars  <- c(NA, names(lagdat2)[grep("BHATT",        names(lagdat2))])
ice.pars    <- c(NA, names(lagdat2)[grep("ICE",          names(lagdat2))])

combos <- tidyr::expand_grid(
  ice    = ice.pars,
  snow   = snow.pars,
  tanner = tanner.pars,
  bhatt  = bhatt.pars
) %>%
  dplyr::filter(!(is.na(bhatt) & is.na(snow) & is.na(ice) & is.na(tanner)))

## 3. Safe wrapper for gamm() ----

safe_gamm <- purrr::safely(gamm)

response <- "HYBRID_ABUND"
k_folds  <- 5

## 4. Fit all combinations with GAMM + AR1 ----

fits <- purrr::pmap_dfr(
  combos,
  function(snow, tanner, ice, bhatt) {
    terms <- c(
      if (!is.na(snow))   paste0("s(", snow,   ", k = 4)") else NULL,
      if (!is.na(tanner)) paste0("s(", tanner, ", k = 4)") else NULL,
      if (!is.na(bhatt))  paste0("s(", bhatt,  ", k = 4)") else NULL,
      if (!is.na(ice))    paste0("s(", ice,    ", k = 4)") else NULL
    )
    
    fml <- as.formula(paste(response, "~", paste(terms, collapse = " + ")))
    
    fit <- safe_gamm(
      formula     = fml,
      data        = lagdat2,
      family      = nb(),
      correlation = corAR1(),
      method      = "REML"
    )
    
    if (!is.null(fit$error)) {
      return(tibble(
        snow_term   = snow,
        tanner_term = tanner,
        ice_term    = ice,
        bhatt_term  = bhatt,
        k_terms     = length(terms),
        AIC         = NA_real_,
        GCV         = NA_real_,
        cv_rmse     = NA_real_,
        edf_total   = NA_real_,
        error       = conditionMessage(fit$error)
      ))
    }
    
    # k-fold CV RMSE using the GAMM
    cv_err <- cv_rmse_gamm(fml, data = lagdat2, k_folds = k_folds)
    
    tibble(
      snow_term   = snow,
      tanner_term = tanner,
      ice_term    = ice,
      bhatt_term  = bhatt,
      k_terms     = length(terms),
      AIC         = AIC(fit$result$lme),             # AIC from lme component
      GCV         = fit$result$gam$gcv.ubre %||% NA, # mgcv GAM GCV/UBRE if present
      cv_rmse     = cv_err,
      edf_total   = sum(fit$result$gam$edf),
      error       = NA_character_
    )
  }
)

fits %>% arrange(cv_rmse, AIC)

              
mod <- gam(HYBRID_ABUND ~ 
            s(TANNER_ABUND_avg2, k=4)+
            s(ICE_avg2lag4, k=4) +
            s(BHATT_avg2, k=4),
            family = nb(),
           data = lagdat2)
