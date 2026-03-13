

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

# HYBRID RESIDUALS ----
hybrid <- readRDS("./Models/hybrid_sdmTMB_tw_90.rda")

plot.resids(hybrid, "Hybrid", "all_sizes") -> out.hybrid

