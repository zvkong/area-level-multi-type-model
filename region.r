source("packages.r")
source("functions.r")
load("westnorthc.Rdata")

choose_pal <- rev(RColorBrewer::brewer.pal(9, "RdBu"))

eps <- 1e-6
p_direct <- pmin(pmax(povrate21$pov_rate, eps), 1 - eps)
mu1_est <- income21.ed$lg_income

Z_1 <- as.numeric(scale(mu1_est))
D <- income21.ed$lg_var / stats::var(mu1_est)
m <- p_direct * (1 - p_direct) / pmax(povrate21$var_pov, 1e-12)
m <- pmax(m, 1e-8)
Z_2 <- m * p_direct

nburn <- 1000
nsim <- 1000
nthin <- 1

fit_shared <- block_sampler_share_RE(
  X_1, X_2, Z_1, Z_2, D, m, S,
  nburn = nburn,
  nsim = nsim,
  nthin = nthin,
  tau_1_init = 1,
  tau_2_init = 1,
  tau_3_init = 1
)

fit_binom_uni <- SRE_binomial(X_2, Z_2, m, S, nburn, nsim, nthin)
fit_gauss_uni <- SRE_sampler(X_1, Z_1, D, S, nburn, nsim, nthin)

est <- extract_estimates(
  fit_shared,
  fit_gauss_uni,
  fit_binom_uni,
  mu1_est = mu1_est
)

plot_facets <- function(sf_obj, fill_vec_list, fill_name, title) {
  labs <- names(fill_vec_list)

  make_one <- function(v, lab) {
    x <- sf_obj
    x[[fill_name]] <- as.numeric(v)
    x$type <- lab
    x
  }

  df_all <- do.call(rbind, Map(make_one, fill_vec_list, labs))

  ggplot(df_all, aes(fill = .data[[fill_name]])) +
    geom_sf(colour = NA) +
    facet_wrap(~type) +
    labs(title = title, fill = fill_name) +
    scale_fill_gradientn(colours = choose_pal, na.value = NA) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold")
    )
}

mu1_shared_mean <- est$mu1_gibbs
mu1_uni_mean <- est$mu1_ind
p_shared_mean <- est$p_gibbs
p_uni_mean <- est$p_ind

p_income <- plot_facets(
  sf_obj = MO_county_sf,
  fill_vec_list = list(
    "Direct estimate" = mu1_est,
    "Univariate model" = mu1_uni_mean,
    "Multi-type model" = mu1_shared_mean
  ),
  fill_name = "log_median_income",
  title = "Log-median income (2017-2021)"
)
print(p_income)

p_pov <- plot_facets(
  sf_obj = povrate21,
  fill_vec_list = list(
    "Direct estimate" = p_direct,
    "Univariate model" = p_uni_mean,
    "Multi-type model" = p_shared_mean
  ),
  fill_name = "pov_rate",
  title = "Poverty rate (2017-2021)"
)
print(p_pov)

mu1_scale <- stats::sd(mu1_est)
var_mu1_shared <- apply(fit_shared$Mu_1.chain, 1, var) * (mu1_scale^2)
var_mu1_uni <- apply(fit_gauss_uni$Mu.chain, 1, var) * (mu1_scale^2)
var_p_shared <- apply(plogis(fit_shared$Mu_2.chain), 1, var)
var_p_uni <- apply(plogis(fit_binom_uni$Mu.chain), 1, var)

p_var_g <- plot_facets(
  sf_obj = MO_county_sf,
  fill_vec_list = list(
    "Direct estimate" = income21.ed$lg_var,
    "Univariate model" = var_mu1_uni,
    "Multi-type model" = var_mu1_shared
  ),
  fill_name = "var_log_income",
  title = "Posterior variance: Gaussian response"
)
print(p_var_g)

p_var_b <- plot_facets(
  sf_obj = povrate21,
  fill_vec_list = list(
    "Direct estimate" = povrate21$var_pov,
    "Univariate model" = var_p_uni,
    "Multi-type model" = var_p_shared
  ),
  fill_name = "var_pov_rate",
  title = "Posterior variance: Binomial response"
)
print(p_var_b)

plot(var_mu1_shared, var_mu1_uni,
     xlab = "Multi-type model variance",
     ylab = "Univariate model variance")
abline(0, 1, col = "red", lty = 2, lwd = 2)

boxplot(
  list(
    "UNIS / Direct" = var_mu1_uni / income21.ed$lg_var,
    "MUTS / Direct" = var_mu1_shared / income21.ed$lg_var,
    "MUTS / UNIS" = var_mu1_shared / var_mu1_uni
  ),
  main = "Variance ratios: Gaussian"
)
abline(h = 1, col = "red", lwd = 2, lty = 2)

plot(var_p_shared, var_p_uni,
     xlab = "Multi-type model variance",
     ylab = "Univariate model variance")
abline(0, 1, col = "red", lty = 2, lwd = 2)

boxplot(
  list(
    "UNIS / Direct" = var_p_uni / povrate21$var_pov,
    "MUTS / Direct" = var_p_shared / povrate21$var_pov,
    "MUTS / UNIS" = var_p_shared / var_p_uni
  ),
  main = "Variance ratios: Binomial"
)
abline(h = 1, col = "red", lwd = 2, lty = 2)

save(
  fit_shared,
  fit_gauss_uni,
  fit_binom_uni,
  mu1_shared_mean,
  mu1_uni_mean,
  p_shared_mean,
  p_uni_mean,
  var_mu1_shared,
  var_mu1_uni,
  var_p_shared,
  var_p_uni,
  file = "region results.RData"
)
