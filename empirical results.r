source("packages.r")
load("SD cleaned.RData")
load("empirical results 0.5.RData")
cor()
if (!exists("choose_pal")) {
  choose_pal <- hcl.colors(9, "YlOrRd", rev = TRUE)
}

make_metric_sf <- function(geo_sf, type, mse, coverage = NULL, is = NULL) {
  geo_sf %>%
    mutate(
      type = type,
      MSE = mse,
      coverage = coverage,
      IS = is
    )
}

plot_spatial <- function(
  df,
  est,
  type = "type",
  title = NULL,
  legend_lab = NULL,
  prob = 0.9,
  ncol = NULL
) {
  est_sym <- rlang::ensym(est)
  type_sym <- rlang::ensym(type)
  est_name <- rlang::as_string(est_sym)

  vals <- dplyr::pull(df, !!est_sym)
  vals <- vals[is.finite(vals)]

  vmin <- min(vals, na.rm = TRUE)
  vmax <- stats::quantile(vals, prob, na.rm = TRUE)

  if (is.null(title)) {
    title <- paste0(est_name, " (", prob * 100, "% truncation)")
  }
  if (is.null(legend_lab)) {
    legend_lab <- est_name
  }

  ggplot(df) +
    geom_sf(aes(fill = !!est_sym), colour = NA) +
    facet_wrap(vars(!!type_sym), ncol = ncol) +
    scale_fill_gradientn(
      colours = choose_pal,
      name = legend_lab,
      limits = c(vmin, vmax),
      oob = scales::squish
    ) +
    labs(
      title = title,
      caption = "Data source: 2017–2021 5-year ACS, U.S. Census Bureau"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}

save_plot <- function(plot_obj, filename, width = 12, height = 4, dpi = 300) {
  ggsave(filename, plot = plot_obj, width = width, height = height, dpi = dpi)
}

income_direct <- make_metric_sf(
  tract_sf,
  type = "Direct estimate",
  mse = mse_dir_g,
  coverage = NA_real_,
  is = NA_real_
)

income_uni <- make_metric_sf(
  tract_sf,
  type = "Univariate model",
  mse = mse_ind_g,
  coverage = cov_ind_g,
  is = IS_ind_g
)

income_multi <- make_metric_sf(
  tract_sf,
  type = "Multi-type model",
  mse = mse_gibbs_g,
  coverage = cov_gibbs_g,
  is = IS_gibbs_g
)

pov_direct <- make_metric_sf(
  tract_sf,
  type = "Direct estimate",
  mse = mse_dir_b,
  coverage = NA_real_,
  is = NA_real_
)

pov_uni <- make_metric_sf(
  tract_sf,
  type = "Univariate model",
  mse = mse_ind_b,
  coverage = cov_ind_b,
  is = IS_ind_b
)

pov_multi <- make_metric_sf(
  tract_sf,
  type = "Multi-type model",
  mse = mse_gibbs_b,
  coverage = cov_gibbs_b,
  is = IS_gibbs_b
)

income_dm <- bind_rows(
  income_direct,
  income_multi
)

pov_dm <- bind_rows(
  pov_direct,
  pov_multi
)

income_model_cmp <- bind_rows(
  income_uni,
  income_multi
)

pov_model_cmp <- bind_rows(
  pov_uni,
  pov_multi
)

p_income_dm <- plot_spatial(
  income_dm,
  est = MSE,
  title = "Gaussian Response: MSE",
  legend_lab = "MSE",
  prob = 0.9
)

p_pov_dm <- plot_spatial(
  pov_dm,
  est = MSE,
  title = "Binomial Response: MSE",
  legend_lab = "MSE",
  prob = 0.9
)

p_income_mse <- plot_spatial(
  income_model_cmp,
  est = MSE,
  title = "Gaussian Response: MSE",
  legend_lab = "MSE",
  prob = 0.9
)

p_pov_mse <- plot_spatial(
  pov_model_cmp,
  est = MSE,
  title = "Binomial Response: MSE",
  legend_lab = "MSE",
  prob = 0.9
)

p_income_is <- plot_spatial(
  income_model_cmp,
  est = IS,
  title = "Gaussian Response: Interval Score",
  legend_lab = "Interval Score",
  prob = 0.9
)

p_pov_is <- plot_spatial(
  pov_model_cmp,
  est = IS,
  title = "Binomial Response: Interval Score",
  legend_lab = "Interval Score",
  prob = 0.92
)

p_income_dm
p_pov_dm
p_income_mse
p_pov_mse
p_income_is
p_pov_is

tab_gaus <- tibble(
  Type = c("Direct estimate", "Univariate model", "Multi-type model"),
  MSE = c(mean(mse_dir_g), mean(mse_ind_g), mean(mse_gibbs_g)),
  Coverage = c(NA_real_, mean(cov_ind_g), mean(cov_gibbs_g)),
  IS = c(NA_real_, mean(IS_ind_g), mean(IS_gibbs_g))
) %>%
  mutate(
    `MSE Red (%)` = c(
      NA_real_,
      (1 - MSE[2] / MSE[1]) * 100,
      (1 - MSE[3] / MSE[1]) * 100
    )
  )

tab_gaus_fmt <- tab_gaus %>%
  mutate(
    MSE = sprintf("%.4f", MSE),
    Coverage = ifelse(is.na(Coverage), "–", sprintf("%.1f%%", 100 * Coverage)),
    IS = ifelse(is.na(IS), "–", sprintf("%.3f", IS)),
    `MSE Red (%)` = ifelse(is.na(`MSE Red (%)`), "–", sprintf("%.2f%%", `MSE Red (%)`))
  )

tab_binom <- tibble(
  Type = c("Direct estimate", "Univariate model", "Multi-type model"),
  MSE = c(mean(mse_dir_b), mean(mse_ind_b), mean(mse_gibbs_b)),
  Coverage = c(NA_real_, mean(cov_ind_b), mean(cov_gibbs_b)),
  IS = c(NA_real_, mean(IS_ind_b), mean(IS_gibbs_b))
) %>%
  mutate(
    `MSE Red (%)` = c(
      NA_real_,
      (1 - MSE[2] / MSE[1]) * 100,
      (1 - MSE[3] / MSE[1]) * 100
    )
  )

tab_binom_fmt <- tab_binom %>%
  mutate(
    MSE = sprintf("%.5f", MSE),
    Coverage = ifelse(is.na(Coverage), "–", sprintf("%.1f%%", 100 * Coverage)),
    IS = ifelse(is.na(IS), "–", sprintf("%.3f", IS)),
    `MSE Red (%)` = ifelse(is.na(`MSE Red (%)`), "–", sprintf("%.2f%%", `MSE Red (%)`))
  )

tab_gaus_fmt
tab_binom_fmt

save_plot(p_income_dm, "p_income_dm.png")
save_plot(p_pov_dm, "p_pov_dm.png")
save_plot(p_income_mse, "p_income_mse.png")
save_plot(p_pov_mse, "p_pov_mse.png")
save_plot(p_income_is, "p_income_is.png")
save_plot(p_pov_is, "p_pov_is.png")