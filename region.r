source("packages.r")
source("functions.r")
load("westnorth_region_clean.RData")

choose_pal <- hcl.colors(9, "YlOrRd", rev = TRUE)

eps <- 1e-6
nburn <- 1000
nsim <- 1000
nthin <- 1
beta_prec <- 0.001

a_tau1 <- 1
a_tau2 <- -1
a_tau3 <- 0

mu1_direct <- income_df$log_income
p_direct <- pmin(pmax(poverty_df$poverty_rate, eps), 1 - eps)
mu2_direct <- qlogis(p_direct)

var_income_direct <- income_df$log_income_variance
var_p_direct <- poverty_df$poverty_rate_variance
var_logit_direct <- var_p_direct / pmax((p_direct * (1 - p_direct))^2, 1e-12)

design_obj <- prep_design(
  mu1 = mu1_direct,
  mu2 = mu2_direct,
  var_D = var_income_direct,
  v = var_logit_direct
)

fits <- run_models(
  X_1 = X_income,
  X_2 = X_poverty,
  Z1 = design_obj$Z1,
  Z2 = design_obj$Z2,
  D = design_obj$D,
  m = design_obj$m,
  S = basis_matrix,
  nburn = nburn,
  nsim = nsim,
  nthin = nthin,
  tau1 = a_tau1,
  tau2 = a_tau2,
  tau3 = a_tau3,
  a_zeta_1 = 2,
  b_zeta_1 = 0.5,
  a_zeta_2 = 2,
  b_zeta_2 = 0.5,
  beta_prec = beta_prec
)

fit_shared <- fits$shared
fit_gauss_uni <- fits$gauss_uni
fit_binom_uni <- fits$binom_uni

est <- extract_estimates(
  fit_shared,
  fit_gauss_uni,
  fit_binom_uni,
  mu1_mean = design_obj$mu1_mean,
  mu1_sd = design_obj$mu1_sd
)

mu1_shared_mean <- est$mu1_gibbs
mu1_uni_mean <- est$mu1_ind
p_shared_mean <- est$p_gibbs
p_uni_mean <- est$p_ind

p_income <- plot_facets(
  sf_obj = region_tract_sf,
  value_list = list(
    "Direct estimate" = mu1_direct,
    "Univariate model" = mu1_uni_mean,
    "Multi-type model" = mu1_shared_mean
  ),
  value_name = "log_median_income",
  title = "Log median income",
  choose_pal = choose_pal
)

p_poverty <- plot_facets(
  sf_obj = region_tract_sf,
  value_list = list(
    "Direct estimate" = p_direct,
    "Univariate model" = p_uni_mean,
    "Multi-type model" = p_shared_mean
  ),
  value_name = "poverty_rate",
  title = "Poverty rate",
  choose_pal = choose_pal
)

print(p_income)
print(p_poverty)

var_mu1_shared <- apply(fit_shared$Mu_1.chain, 1, stats::var) * (design_obj$mu1_sd^2)
var_mu1_uni <- apply(fit_gauss_uni$Mu.chain, 1, stats::var) * (design_obj$mu1_sd^2)
var_p_shared <- apply(plogis(fit_shared$Mu_2.chain), 1, stats::var)
var_p_uni <- apply(plogis(fit_binom_uni$Mu.chain), 1, stats::var)

p_var_income <- plot_facets(
  sf_obj = region_tract_sf,
  value_list = list(
    "Direct estimate" = var_income_direct,
    "Univariate model" = var_mu1_uni,
    "Multi-type model" = var_mu1_shared
  ),
  value_name = "var_log_income",
  title = "Posterior variance: Gaussian response",
  choose_pal = choose_pal
)

p_var_poverty <- plot_facets(
  sf_obj = region_tract_sf,
  value_list = list(
    "Direct estimate" = var_p_direct,
    "Univariate model" = var_p_uni,
    "Multi-type model" = var_p_shared
  ),
  value_name = "var_poverty_rate",
  title = "Posterior variance: Binomial response",
  choose_pal = choose_pal
)

print(p_var_income)
print(p_var_poverty)

plot(
  var_mu1_shared,
  var_mu1_uni,
  xlab = "Multi-type model variance",
  ylab = "Univariate model variance",
  main = "Gaussian response"
)
abline(0, 1, col = "red", lty = 2, lwd = 2)

boxplot(
  list(
    "UNIS / Direct" = var_mu1_uni / var_income_direct,
    "MUTS / Direct" = var_mu1_shared / var_income_direct,
    "MUTS / UNIS" = var_mu1_shared / var_mu1_uni
  ),
  main = "Variance ratios: Gaussian"
)
abline(h = 1, col = "red", lwd = 2, lty = 2)

plot(
  var_p_shared,
  var_p_uni,
  xlab = "Multi-type model variance",
  ylab = "Univariate model variance",
  main = "Binomial response"
)
abline(0, 1, col = "red", lty = 2, lwd = 2)

boxplot(
  list(
    "UNIS / Direct" = var_p_uni / var_p_direct,
    "MUTS / Direct" = var_p_shared / var_p_direct,
    "MUTS / UNIS" = var_p_shared / var_p_uni
  ),
  main = "Variance ratios: Binomial"
)
abline(h = 1, col = "red", lwd = 2, lty = 2)

save(
  fit_shared,
  fit_gauss_uni,
  fit_binom_uni,
  design_obj,
  mu1_direct,
  p_direct,
  mu1_shared_mean,
  mu1_uni_mean,
  p_shared_mean,
  p_uni_mean,
  var_income_direct,
  var_p_direct,
  var_mu1_shared,
  var_mu1_uni,
  var_p_shared,
  var_p_uni,
  basis_q,
  beta_prec,
  p_income,
  p_poverty,
  p_var_income,
  p_var_poverty,
  file = "region_results.RData"
)
