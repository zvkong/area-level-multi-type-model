load("./westnorthc.Rdata")           # expects MO_county_sf, povrate21, income21.ed, X_1, X_2, S
source("packages.r")     # region_analysis_clean.R
source("functions.r")

choose_pal <- rev(RColorBrewer::brewer.pal(9, "RdBu"))

# true values and design
p_true    <- povrate21$pov_rate
mu1_true  <- income21.ed$lg_income
mu2_true  <- qlogis(p_true)
n_area    <- length(mu1_true)

var_D  <- income21.ed$lg_var
v_rate <- povrate21$var_pov
v      <- (1 / (p_true * (1 - p_true)))^2 * v_rate

Z_1 <- as.numeric(scale(mu1_true))
D   <- var_D / var(mu1_true)

pi_hat   <- plogis(mu2_true)
v_pi_hat <- (pi_hat * (1 - pi_hat))^2 * v
m        <- pi_hat * (1 - pi_hat) / pmax(v_pi_hat, 1e-12)
Z_2      <- m * pi_hat

# MCMC
nburn <- 1000
nsim  <- 1000
nthin <- 1

# fits (standardized names)
fit_shared     <- block_sampler_share_RE(X_1, X_2, Z_1, Z_2, D, m, S,
                                         nburn = nburn, nsim = nsim, nthin = nthin,
                                         tau_1_init = 1, tau_2_init = 1, tau_3_init = 1)
fit_binom_uni  <- SRE_binomial(X_2, Z_2, m, S, nburn, nsim, nthin)
fit_gauss_uni  <- SRE_sampler (X_1, Z_1, D, S, nburn, nsim, nthin)

# plotting helper
plot_facets <- function(sf_obj, fill_vec_list, fill_name, title) {
  stopifnot(inherits(sf_obj, "sf"))
  labs <- names(fill_vec_list)
  stopifnot(all(vapply(fill_vec_list, length, 0L) == nrow(sf_obj)))
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
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"))
}

# means
mu1_shared_mean <- rowMeans(fit_shared$Mu_1.chain)
mu1_uni_mean    <- rowMeans(fit_gauss_uni$Mu.chain)
p_shared_mean   <- plogis(rowMeans(fit_shared$Mu_2.chain))
p_uni_mean      <- plogis(rowMeans(fit_binom_uni$Mu.chain))

p_income <- plot_facets(
  sf_obj = MO_county_sf,
  fill_vec_list = list(
    "Direct Estimate"           = income21.ed$income_scaled,
    "Predict By Shared Random"  = mu1_shared_mean,
    "Predict By Spatial Random" = mu1_uni_mean
  ),
  fill_name = "median_income",
  title = "Median Income (2017–2021)"
)
print(p_income)

p_pov <- plot_facets(
  sf_obj = povrate21,
  fill_vec_list = list(
    "Direct Estimate"           = p_true,
    "Predict By Shared Random"  = p_shared_mean,
    "Predict By Spatial Random" = p_uni_mean
  ),
  fill_name = "pov_rate",
  title = "Poverty Rate (2017–2021)"
)
print(p_pov)

# variances
var_mu1_shared <- apply(fit_shared$Mu_1.chain, 1, var)
var_mu1_uni    <- apply(fit_gauss_uni$Mu.chain, 1, var)
var_p_shared   <- apply(plogis(fit_shared$Mu_2.chain), 1, var)
var_p_uni      <- apply(plogis(fit_binom_uni$Mu.chain), 1, var)

p_var_g <- plot_facets(
  sf_obj = MO_county_sf,
  fill_vec_list = list(
    "Direct Estimate"           = D,
    "Predict By Shared Random"  = var_mu1_shared,
    "Predict By Spatial Random" = var_mu1_uni
  ),
  fill_name = "var_mu1",
  title = "Posterior Variance: Gaussian"
)
print(p_var_g)

p_var_b <- plot_facets(
  sf_obj = povrate21,
  fill_vec_list = list(
    "Direct Estimate"           = povrate21$var_pov,
    "Predict By Shared Random"  = var_p_shared,
    "Predict By Spatial Random" = var_p_uni
  ),
  fill_name = "var_pov",
  title = "Posterior Variance: Binomial"
)
print(p_var_b)

# quick numeric checks
plot(var_mu1_shared, var_mu1_uni,
     xlab = "Var (Shared Random)", ylab = "Var (Univariate)")
abline(0, 1, col = "red", lty = 2, lwd = 2)

boxplot(var_mu1_shared / var_mu1_uni,
        main = "Var ratio (Shared / Univariate) — Gaussian")
abline(h = 1, col = 'red', lwd = 2, lty = 2)

boxplot(var_mu1_shared / D,
        main = "Var ratio (Shared / Direct) — Gaussian")
abline(h = 1, col = 'red', lwd = 2, lty = 2)

plot(var_p_shared, var_p_uni,
     xlab = "Var (Shared Random)", ylab = "Var (Univariate)")
abline(0, 1, col = "red", lty = 2, lwd = 2)

boxplot(var_p_shared / var_p_uni,
        main = "Var ratio (Shared / Univariate) — Binomial")
abline(h = 1, col = 'red', lwd = 2, lty = 2)

boxplot(var_p_shared / povrate21$var_pov,
        main = "Var ratio (Shared / Direct) — Binomial")
abline(h = 1, col = 'red', lwd = 2, lty = 2)
