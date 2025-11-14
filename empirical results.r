# spatial plot already loaded in the rdata
source("packages.r")
load("empirical results.RData")

# metrics for combined income
direct_income <- MO_county_sf %>%
  mutate(
    type     = "Direct estimate",
    MSE      = mse_dir_g,
    Bias     = bias_dir_g,
    coverage = NA_real_,
    IS       = NA_real_
  )

unis_income <- MO_county_sf %>%
  mutate(
    type     = "Univariate model",
    MSE      = mse_ind_g,
    Bias     = bias_ind_g,
    coverage = cov_ind_g,
    IS       = IS_ind_g
  )

muts_income <- MO_county_sf %>%
  mutate(
    type     = "Multi-type model",
    MSE      = mse_gibbs_g,
    Bias     = bias_gibbs_g,
    coverage = cov_gibbs_g,
    IS       = IS_gibbs_g
  )

# metrics for combined poverty
direct_pov <- povrate21 %>%
  mutate(
    type     = "Direct estimate",
    MSE      = mse_dir_b,
    Bias     = bias_dir_b,
    coverage = NA_real_,
    IS       = NA_real_
  )

unis_pov <- povrate21 %>%
  mutate(
    type     = "Univariate model",
    MSE      = mse_ind_b,
    Bias     = bias_ind_b,
    coverage = cov_ind_b,
    IS       = IS_ind_b
  )

muts_pov <- povrate21 %>%
  mutate(
    type     = "Multi-type model",
    MSE      = mse_gibbs_b,
    Bias     = bias_gibbs_b,
    coverage = cov_gibbs_b,
    IS       = IS_gibbs_b
  )

# plot function
plot_spatial <- function(df,
                         est,
                         type = "type",
                         title = NULL,
                         legend_lab = NULL,
                         prob = 0.85) {

  est_sym  <- ensym(est)
  type_sym <- ensym(type)
  est_name <- as_string(est_sym)

  vals  <- df %>% pull(!!est_sym)
  vmin  <- min(vals, na.rm = TRUE)
  vmax  <- quantile(vals, prob, na.rm = TRUE)

  if (is.null(title)) {
    title <- paste0(est_name, " (", prob*100, "% truncation)")
  }
  if (is.null(legend_lab)) {
    legend_lab <- est_name
  }

  ggplot(df) +
    geom_sf(aes(fill = !!est_sym), colour = NA) +
    facet_wrap(vars(!!type_sym)) +
    scale_fill_gradientn(
      colours = choose_pal,
      name    = legend_lab,
      limits  = c(vmin, vmax),
      oob     = scales::squish
    ) +
    labs(
      title   = title,
      caption = "Data source: 2017–2021 5-year ACS, US Census Bureau"
    ) +
    theme_minimal()
}

# direc mse comparison
gaus_dm <- bind_rows(
  direct_income %>% mutate(type = "Direct estimate"),
  muts_income %>% mutate(type = "Multi-type model")
)
bion_dm <- bind_rows(
  direct_pov %>% mutate(type = "Direct estimate"),
  muts_pov %>% mutate(type = "Multi-type model")
)
p_gaus_dm <- plot_spatial(gaus_dm, est = MSE,
title = "Gaussian MSE", legend_lab = "MSE", prob = 0.9)
p_bion_dm <- plot_spatial(bion_dm, est = MSE,
title = "Binomial MSE", legend_lab = "MSE", prob = 0.9)
p_gaus_dm
p_bion_dm
# gaus mse comparison 
gaus_df <- bind_rows(
  unis_income %>% mutate(type = "Univariate model"),
  muts_income %>% mutate(type = "Multi-type model")
)
p_gaus <- plot_spatial(gaus_df, est = MSE,
title = "Gaussian MSE", legend_lab = "MSE")

# binom mse comparison
bio_df <- bind_rows(
  unis_pov %>% mutate(type = "Univariate model"),
  muts_pov %>% mutate(type = "Multi-type model")
)
p_bio <- plot_spatial(bio_df, est = MSE,
title = "Binomial MSE", legend_lab = "MSE", prob = 0.9)
p_bio

# gasu interval score comparison
p_gaus_is <- plot_spatial(gaus_df, est = IS,
title = "Gaussian Interval Score", legend_lab = "Interval Score")

# binom interval score comparison
p_bio_is <- plot_spatial(bio_df, est = IS,
title = "Binomial Interval Score", legend_lab = "Interval Score", prob = 0.92)
p_bio_is
# output the results of mse comparison and interval score comparison

## Gaussian
tab_gaus <- tibble(
  Type = c("Direct estimate", "UNIS model", "MUTS model"),
  MSE  = c(mean(mse_dir_g),
           mean(mse_ind_g),
           mean(mse_gibbs_g)),
  Coverage = c(NA,
               mean(cov_ind_g),
               mean(cov_gibbs_g)),
  IS   = c(NA,
           mean(IS_ind_g),
           mean(IS_gibbs_g))
) %>%
  mutate(
    `MSE Red (%)` = c(
      NA,
      (1 - MSE[2] / MSE[1]) * 100,
      (1 - MSE[3] / MSE[1]) * 100
    )
  )

tab_gaus_fmt <- tab_gaus %>%
  mutate(
    MSE        = sprintf("%.4f", MSE),
    Coverage   = ifelse(is.na(Coverage), "–",
                        sprintf("%.1f%%", Coverage * 100)),
    IS         = ifelse(is.na(IS), "–",
                        sprintf("%.3f", IS)),
    `MSE Red (%)` = ifelse(is.na(`MSE Red (%)`), "–",
                           sprintf("%.2f%%", `MSE Red (%)`))
  )

## Binomial
tab_binom <- tibble(
  Type = c("Direct estimate", "UNIS model", "MUTS model"),
  MSE  = c(mean(mse_dir_b),
           mean(mse_ind_b),
           mean(mse_gibbs_b)),
  Coverage = c(NA,
               mean(cov_ind_b),
               mean(cov_gibbs_b)),
  IS   = c(NA,
           mean(IS_ind_b),
           mean(IS_gibbs_b))
) %>%
  mutate(
    `MSE Red (%)` = c(
      NA,
      (1 - MSE[2] / MSE[1]) * 100,
      (1 - MSE[3] / MSE[1]) * 100
    )
  )


tab_binom_fmt <- tab_binom %>%
  mutate(
    MSE        = sprintf("%.5f", MSE),
    Coverage   = ifelse(is.na(Coverage), "–",
                        sprintf("%.1f%%", Coverage * 100)),
    IS         = ifelse(is.na(IS), "–",
                        sprintf("%.3f", IS)),
    `MSE Red (%)` = ifelse(is.na(`MSE Red (%)`), "–",
                           sprintf("%.2f%%", `MSE Red (%)`))
  )

tab_gaus_fmt
tab_binom_fmt

for (p in c("p_gaus_dm","p_bion_dm","p_gaus","p_bio","p_gaus_is","p_bio_is"))
  ggsave(paste0(p, ".png"), plot = get(p), width = 12, height = 4, dpi = 300)
