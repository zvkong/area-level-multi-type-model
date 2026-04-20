source("packages.r")
load("SD cleaned.RData")
load("empirical results.RData")

library(patchwork)

out_dir <- "main_text_figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_plot <- function(plot_obj, filename, width = 8, height = 5, dpi = 300) {
  ggsave(
    file.path(out_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
}

theme_main <- theme_minimal(base_size = 11, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14, family = "Arial"),
    plot.subtitle = element_text(size = 10, family = "Arial"),
    axis.title = element_text(size = 11, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 10, family = "Arial"),
    legend.text = element_text(size = 9, family = "Arial"),
    strip.text = element_text(face = "bold", size = 11, family = "Arial"),
    plot.margin = margin(6, 8, 6, 8)
  )

map_theme <- theme_main +
  theme(
    panel.grid.major = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 9, family = "Arial")
  )

resp_cols <- c(
  "Binomial response" = "#E6A04B",
  "Gaussian response" = "#6F8FB8"
)

make_sf_metric <- function(sf_obj, values, value_name) {
  stopifnot(nrow(sf_obj) == length(values))
  out <- sf_obj
  out[[value_name]] <- as.numeric(values)
  out
}

plot_reduction_map <- function(sf_obj, values, value_name, title, prob = 0.98) {
  plot_df <- make_sf_metric(sf_obj, values, value_name)
  vals <- plot_df[[value_name]]
  vals <- vals[is.finite(vals)]
  vmax <- stats::quantile(abs(vals), prob, na.rm = TRUE)
  vmax <- max(vmax, 1e-8)

  ggplot(plot_df) +
    geom_sf(aes(fill = .data[[value_name]]), colour = NA) +
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "grey96",
      high = "#2A9D8F",
      midpoint = 0,
      limits = c(-vmax, vmax),
      oob = scales::squish,
      na.value = "grey90",
      name = "Reduction (%)"
    ) +
    labs(
      title = title,
      subtitle = "Positive values favor the Multi-type model."
    ) +
    map_theme
}

plot_violin_reduction <- function(df, yvar, title, ylab) {
  ggplot(df, aes(x = Response, y = .data[[yvar]], fill = Response)) +
    geom_violin(alpha = 0.9, trim = FALSE, linewidth = 0.4) +
    geom_boxplot(width = 0.14, outlier.size = 0.5, alpha = 0.95) +
    geom_hline(yintercept = 0, linetype = 2, colour = "red", linewidth = 0.8) +
    scale_fill_manual(values = resp_cols) +
    labs(
      title = title,
      subtitle = "Positive values favor the Multi-type model.",
      x = NULL,
      y = ylab
    ) +
    theme_main +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
}

mse_red_pct_g <- 100 * (1 - mse_gibbs_g / mse_ind_g)
mse_red_pct_b <- 100 * (1 - mse_gibbs_b / mse_ind_b)

is_red_pct_g <- 100 * (1 - IS_gibbs_g / IS_ind_g)
is_red_pct_b <- 100 * (1 - IS_gibbs_b / IS_ind_b)

violin_df <- dplyr::bind_rows(
  tibble(
    Response = "Binomial response",
    mse_reduction = mse_red_pct_b,
    is_reduction = is_red_pct_b
  ),
  tibble(
    Response = "Gaussian response",
    mse_reduction = mse_red_pct_g,
    is_reduction = is_red_pct_g
  )
)

violin_df$Response <- factor(
  violin_df$Response,
  levels = c("Binomial response", "Gaussian response")
)

p1 <- plot_violin_reduction(
  violin_df,
  yvar = "mse_reduction",
  title = "A. Tract-level MSE reduction",
  ylab = "MSE reduction (%)"
)

p2 <- plot_violin_reduction(
  violin_df,
  yvar = "is_reduction",
  title = "B. Tract-level interval score reduction",
  ylab = "Interval score reduction (%)"
)

p3 <- plot_reduction_map(
  tract_sf,
  mse_red_pct_g,
  "mse_red_pct_g",
  "C. Gaussian response: tract-level MSE reduction"
)

p4 <- plot_reduction_map(
  tract_sf,
  mse_red_pct_b,
  "mse_red_pct_b",
  "D. Binomial response: tract-level MSE reduction"
)

final_fig <- (p1 | p2) / (p3 | p4) +
  plot_layout(heights = c(1, 1.1))

print(final_fig)

save_plot(p1, "panel_A_violin_mse_reduction.png", width = 7.2, height = 4.8)
save_plot(p2, "panel_B_violin_is_reduction.png", width = 7.2, height = 4.8)
save_plot(p3, "panel_C_map_mse_reduction_gaussian.png", width = 7.8, height = 5.2)
save_plot(p4, "panel_D_map_mse_reduction_binomial.png", width = 7.8, height = 5.2)

ggsave(
  file.path(out_dir, "main_text_four_panel_figure.png"),
  plot = final_fig,
  width = 14,
  height = 11,
  dpi = 300
)

ggsave(
  file.path(out_dir, "main_text_four_panel_figure.pdf"),
  plot = final_fig,
  width = 14,
  height = 11,
  device = grDevices::cairo_pdf
)