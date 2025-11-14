# simulation_clean.R

suppressPackageStartupMessages({ library(spdep) })

source("packages.r")
load("SD cleaned.RData")          # expects MO_county_sf, income21.ed, povrate21, X_white, etc.
source("functions.r")

# basis
MO_nb  <- poly2nb(MO_county_sf, queen = TRUE)
W_list <- nb2listw(MO_nb, style = "B", zero.policy = TRUE)
B      <- as.matrix(as_dgRMatrix_listw(W_list))
eig    <- eigen(B, symmetric = TRUE)
idx    <- large_abs(eig$values, 0.25)
S      <- eig$vectors[, idx, drop = FALSE]

# fixed quantities
X_1 <- as.matrix(cbind(1, X_white))
X_2 <- as.matrix(cbind(1, X_white))
Mu1_true <- income21.ed$lg_income
p_true   <- povrate21$pov_rate
Mu2_true <- qlogis(p_true)
N        <- length(Mu1_true)
var_D    <- income21.ed$lg_var
v_rate   <- povrate21$var_pov
v        <- (1/(p_true*(1 - p_true)))^2 * v_rate

# mcmc config
n_data <- 100
nburn  <- 1000
nsim   <- 1000
nthin  <- 1

# containers
Mu1_dir <- Mu2_dir <- Z1_set <- Z2_set <- Pi_set <- matrix(0, N, n_data)
G_gibbs <- G_ind   <- P_gibbs <- P_ind  <- matrix(0, N, n_data)
GG_l <- GG_u <- GI_l <- GI_u <- GB_l <- GB_u <- IB_l <- IB_u <- matrix(0, N, n_data)

# simulation loop
for (i in 1:n_data) {
  set.seed(i)
  drw <- draw_direct(N, Mu1_true, Mu2_true, var_D, v)
  des <- prep_design(drw$mu1, drw$mu2, var_D, v)

  Mu1_dir[, i] <- drw$mu1
  Mu2_dir[, i] <- drw$mu2
  Z1_set[, i]  <- des$Z1
  Z2_set[, i]  <- des$Z2
  Pi_set[, i]  <- des$pi_hat

  fits <- run_models(X_1, X_2, des$Z1, des$Z2, des$D, des$m, S,
                     nburn, nsim, nthin, tau1 = 2, tau2 = -2, tau3 = 0)

  est <- extract_estimates(fits$shared, fits$gaus, fits$bio, drw$mu1)

  GG_l[, i] <- est$gibbs_g_l; GG_u[, i] <- est$gibbs_g_u
  GI_l[, i] <- est$ind_g_l;   GI_u[, i] <- est$ind_g_u
  GB_l[, i] <- est$gibbs_b_l; GB_u[, i] <- est$gibbs_b_u
  IB_l[, i] <- est$ind_b_l;   IB_u[, i] <- est$ind_b_u

  G_gibbs[, i] <- est$mu1_gibbs; G_ind[, i] <- est$mu1_ind
  P_gibbs[, i] <- est$p_gibbs;   P_ind[, i] <- est$p_ind
}

# metrics
mse_dir_g   <- vapply(1:N, function(i) MSE(Mu1_true[i], Mu1_dir[i,]), 0.0)
mse_gibbs_g <- vapply(1:N, function(i) MSE(Mu1_true[i], G_gibbs[i,]), 0.0)
mse_ind_g   <- vapply(1:N, function(i) MSE(Mu1_true[i], G_ind[i,]),   0.0)

mse_dir_b   <- vapply(1:N, function(i) MSE(p_true[i],   Pi_set[i,]),  0.0)
mse_gibbs_b <- vapply(1:N, function(i) MSE(p_true[i],   P_gibbs[i,]), 0.0)
mse_ind_b   <- vapply(1:N, function(i) MSE(p_true[i],   P_ind[i,]),   0.0)

bias_dir_g   <- abs(rowMeans(Mu1_dir) - Mu1_true)
bias_gibbs_g <- abs(rowMeans(G_gibbs) - Mu1_true)
bias_ind_g   <- abs(rowMeans(G_ind)   - Mu1_true)

bias_dir_b   <- abs(rowMeans(Pi_set)  - p_true)
bias_gibbs_b <- abs(rowMeans(P_gibbs) - p_true)
bias_ind_b   <- abs(rowMeans(P_ind)   - p_true)

Mu1_mat <- matrix(Mu1_true, nrow = N, ncol = n_data)
p_mat   <- matrix(p_true,   nrow = N, ncol = n_data)

cov_gibbs_g <- coverage_count(GG_l, GG_u, Mu1_mat) / n_data
cov_ind_g   <- coverage_count(GI_l, GI_u, Mu1_mat) / n_data
cov_gibbs_b <- coverage_count(GB_l, GB_u, p_mat)   / n_data
cov_ind_b   <- coverage_count(IB_l, IB_u, p_mat)   / n_data

IS_gibbs_g  <- rowMeans(interval_score(GG_l, GG_u, Mu1_mat))
IS_ind_g    <- rowMeans(interval_score(GI_l, GI_u, Mu1_mat))
IS_gibbs_b  <- rowMeans(interval_score(GB_l, GB_u, p_mat))
IS_ind_b    <- rowMeans(interval_score(IB_l, IB_u, p_mat))

# summary
cat("\nGaussian\n")
cat("mean(MSE) Direct / Shared / Univariate: ",
    mean(mse_dir_g), " / ", mean(mse_gibbs_g), " / ", mean(mse_ind_g), "\n", sep = "")
cat("mean(Bias) Direct / Shared / Univariate: ",
    mean(bias_dir_g), " / ", mean(bias_gibbs_g), " / ", mean(bias_ind_g), "\n", sep = "")
cat("mean(Coverage) Shared / Univariate: ",
    mean(cov_gibbs_g), " / ", mean(cov_ind_g), "\n", sep = "")
cat("mean(Interval Score) Shared / Univariate: ",
    mean(IS_gibbs_g), " / ", mean(IS_ind_g), "\n", sep = "")

cat("\nBinomial\n")
cat("mean(MSE) Direct / Shared / Univariate: ",
    mean(mse_dir_b), " / ", mean(mse_gibbs_b), " / ", mean(mse_ind_b), "\n", sep = "")
cat("mean(Bias) Direct / Shared / Univariate: ",
    mean(bias_dir_b), " / ", mean(bias_gibbs_b), " / ", mean(bias_ind_b), "\n", sep = "")
cat("mean(Coverage) Shared / Univariate: ",
    mean(cov_gibbs_b), " / ", mean(cov_ind_b), "\n", sep = "")
cat("mean(Interval Score) Shared / Univariate: ",
    mean(IS_gibbs_b), " / ", mean(IS_ind_b), "\n", sep = "")
