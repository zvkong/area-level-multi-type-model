suppressPackageStartupMessages({
  library(spdep)
  library(sf)
  library(dplyr)
})

source("packages.r")
source("functions.r")
load("SD cleaned.RData")

tract_nb <- poly2nb(tract_sf, queen = TRUE)
W_list <- nb2listw(tract_nb, style = "B", zero.policy = TRUE)
B <- as.matrix(as_dgRMatrix_listw(W_list))

eig <- eigen(B, symmetric = TRUE)
idx <- large_abs(eig$values, .25)
S <- eig$vectors[, idx, drop = FALSE]
X_1 <- as.matrix(cbind(1, X_white))
X_2 <- as.matrix(cbind(1, X_white))

eps <- 1e-6
Mu1_true <- income21.ed$lg_income
p_true <- pmin(pmax(povrate21$pov_rate, eps), 1 - eps)
Mu2_true <- qlogis(p_true)
N <- length(Mu1_true)
inflation_factor <- 1
var_D <- income21.ed$lg_var * inflation_factor
v_rate <- povrate21$var_pov
v <- ((1 / (p_true * (1 - p_true)))^2 * v_rate) * inflation_factor

n_data <- 100
nburn <- 1000
nsim <- 1000
nthin <- 1

Mu1_dir <- matrix(0, N, n_data)
Mu2_dir <- matrix(0, N, n_data)
Z1_set <- matrix(0, N, n_data)
Z2_set <- matrix(0, N, n_data)
Pi_set <- matrix(0, N, n_data)

G_gibbs <- matrix(0, N, n_data)
G_ind <- matrix(0, N, n_data)
P_gibbs <- matrix(0, N, n_data)
P_ind <- matrix(0, N, n_data)

GG_l <- matrix(0, N, n_data)
GG_u <- matrix(0, N, n_data)
GI_l <- matrix(0, N, n_data)
GI_u <- matrix(0, N, n_data)
GB_l <- matrix(0, N, n_data)
GB_u <- matrix(0, N, n_data)
IB_l <- matrix(0, N, n_data)
IB_u <- matrix(0, N, n_data)

for (i in seq_len(n_data)) {
  cat("\n==================================================\n")
  cat(sprintf("   Starting Dataset %d / %d\n", i, n_data))
  cat("==================================================\n")

  set.seed(i)

  # drw <- draw_direct_independent(Mu1_true, Mu2_true, var_D, v)
  drw <- draw_direct_correlated(
    Mu1_true,
    Mu2_true,
    var_D,
    v,
    rho = 0.9)
  des <- prep_design(drw$mu1, drw$mu2, var_D, v)

  Mu1_dir[, i] <- drw$mu1
  Mu2_dir[, i] <- drw$mu2
  Z1_set[, i] <- des$Z1
  Z2_set[, i] <- des$Z2
  Pi_set[, i] <- des$pi_hat

  fits <- run_models(
    X_1 = X_1,
    X_2 = X_2,
    Z1 = des$Z1,
    Z2 = des$Z2,
    D = des$D,
    m = des$m,
    S = S,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    tau1 = 1,
    tau2 = -1,
    tau3 = 0,
    a_zeta_1 = 3,
    b_zeta_1 = 0.5,
    a_zeta_2 = 3,
    b_zeta_2 = 0.5,
    beta_prec = 0.01
  )

  est <- extract_estimates(
    fits$shared,
    fits$gauss_uni,
    fits$binom_uni,
    mu1_mean = des$mu1_mean,
    mu1_sd = des$mu1_sd
  )

  GG_l[, i] <- est$gibbs_g_l
  GG_u[, i] <- est$gibbs_g_u
  GI_l[, i] <- est$ind_g_l
  GI_u[, i] <- est$ind_g_u
  GB_l[, i] <- est$gibbs_b_l
  GB_u[, i] <- est$gibbs_b_u
  IB_l[, i] <- est$ind_b_l
  IB_u[, i] <- est$ind_b_u

  G_gibbs[, i] <- est$mu1_gibbs
  G_ind[, i] <- est$mu1_ind
  P_gibbs[, i] <- est$p_gibbs
  P_ind[, i] <- est$p_ind
}

cat("\n\nAll datasets completed. Calculating evaluation metrics...\n")
length(Mu1_true)
mse_dir_g <- vapply(seq_len(N), function(i) MSE(Mu1_true[i], Mu1_dir[i, ]), 0.0)
mse_gibbs_g <- vapply(seq_len(N), function(i) MSE(Mu1_true[i], G_gibbs[i, ]), 0.0)
mse_ind_g <- vapply(seq_len(N), function(i) MSE(Mu1_true[i], G_ind[i, ]), 0.0)

mse_dir_b <- vapply(seq_len(N), function(i) MSE(p_true[i], Pi_set[i, ]), 0.0)
mse_gibbs_b <- vapply(seq_len(N), function(i) MSE(p_true[i], P_gibbs[i, ]), 0.0)
mse_ind_b <- vapply(seq_len(N), function(i) MSE(p_true[i], P_ind[i, ]), 0.0)

Mu1_mat <- matrix(Mu1_true, nrow = N, ncol = n_data)
p_mat <- matrix(p_true, nrow = N, ncol = n_data)

cov_gibbs_g <- coverage_count(GG_l, GG_u, Mu1_mat) / n_data
cov_ind_g <- coverage_count(GI_l, GI_u, Mu1_mat) / n_data
cov_gibbs_b <- coverage_count(GB_l, GB_u, p_mat) / n_data
cov_ind_b <- coverage_count(IB_l, IB_u, p_mat) / n_data

IS_gibbs_g <- rowMeans(interval_score(GG_l, GG_u, Mu1_mat))
IS_ind_g <- rowMeans(interval_score(GI_l, GI_u, Mu1_mat))
IS_gibbs_b <- rowMeans(interval_score(GB_l, GB_u, p_mat))
IS_ind_b <- rowMeans(interval_score(IB_l, IB_u, p_mat))

tab_gaus_fmt <- data.frame(
  Type = c("Direct estimate", "Univariate model", "Multi-type model"),
  MSE = c(mean(mse_dir_g), mean(mse_ind_g), mean(mse_gibbs_g)) * 1e3,
  Coverage = c(NA_real_, mean(cov_ind_g), mean(cov_gibbs_g)),
  IS = c(NA_real_, mean(IS_ind_g), mean(IS_gibbs_g)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tab_gaus_fmt$`MSE Red (%)` <- c(
  NA_real_,
  (1 - tab_gaus_fmt$MSE[2] / tab_gaus_fmt$MSE[1]) * 100,
  (1 - tab_gaus_fmt$MSE[3] / tab_gaus_fmt$MSE[1]) * 100
)

tab_gaus_fmt$MSE <- sprintf("%.1f", tab_gaus_fmt$MSE)
tab_gaus_fmt$Coverage <- ifelse(
  is.na(tab_gaus_fmt$Coverage),
  "-",
  sprintf("%.1f%%", 100 * as.numeric(tab_gaus_fmt$Coverage))
)
tab_gaus_fmt$IS <- ifelse(
  is.na(tab_gaus_fmt$IS),
  "-",
  sprintf("%.3f", as.numeric(tab_gaus_fmt$IS))
)
tab_gaus_fmt$`MSE Red (%)` <- ifelse(
  is.na(tab_gaus_fmt$`MSE Red (%)`),
  "-",
  sprintf("%.2f%%", as.numeric(tab_gaus_fmt$`MSE Red (%)`))
)

tab_binom_fmt <- data.frame(
  Type = c("Direct estimate", "Univariate model", "Multi-type model"),
  MSE = c(mean(mse_dir_b), mean(mse_ind_b), mean(mse_gibbs_b)) * 1e3,
  Coverage = c(NA_real_, mean(cov_ind_b), mean(cov_gibbs_b)),
  IS = c(NA_real_, mean(IS_ind_b), mean(IS_gibbs_b)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tab_binom_fmt$`MSE Red (%)` <- c(
  NA_real_,
  (1 - tab_binom_fmt$MSE[2] / tab_binom_fmt$MSE[1]) * 100,
  (1 - tab_binom_fmt$MSE[3] / tab_binom_fmt$MSE[1]) * 100
)

tab_binom_fmt$MSE <- sprintf("%.2f", tab_binom_fmt$MSE)
tab_binom_fmt$Coverage <- ifelse(
  is.na(tab_binom_fmt$Coverage),
  "-",
  sprintf("%.1f%%", 100 * as.numeric(tab_binom_fmt$Coverage))
)
tab_binom_fmt$IS <- ifelse(
  is.na(tab_binom_fmt$IS),
  "-",
  sprintf("%.3f", as.numeric(tab_binom_fmt$IS))
)
tab_binom_fmt$`MSE Red (%)` <- ifelse(
  is.na(tab_binom_fmt$`MSE Red (%)`),
  "-",
  sprintf("%.2f%%", as.numeric(tab_binom_fmt$`MSE Red (%)`))
)

cat("\n--- Gaussian Response (Log Income) ---\n")
print(tab_gaus_fmt)

cat("\n--- Binomial Response (Poverty Rate) ---\n")
print(tab_binom_fmt)

save(
  Mu1_true,
  p_true,
  Mu2_true,
  var_D,
  v,
  mse_dir_g,
  mse_gibbs_g,
  mse_ind_g,
  mse_dir_b,
  mse_gibbs_b,
  mse_ind_b,
  cov_gibbs_g,
  cov_ind_g,
  cov_gibbs_b,
  cov_ind_b,
  IS_gibbs_g,
  IS_ind_g,
  IS_gibbs_b,
  IS_ind_b,
  tab_gaus_fmt,
  tab_binom_fmt,
  file = "empirical results 3.RData"
)