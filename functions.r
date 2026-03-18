large_abs <- function(v, prop = 0.25) {
  n_keep <- max(1L, ceiling(prop * length(v)))
  order(abs(v), decreasing = TRUE)[seq_len(n_keep)]
}

MSE <- function(real, prediction) {
  mean((real - prediction)^2)
}

rmvn_once <- function(mu, Sigma) {
  Sigma <- Matrix::forceSymmetric(Sigma)
  L <- Matrix::chol(Sigma)
  z <- matrix(stats::rnorm(length(mu)), ncol = 1)
  as.vector(mu + L %*% z)
}

block_sampler_share_RE <- function(
  X_1, X_2, Z_1, Z_2, D, m, S,
  nburn, nsim, nthin,
  tau_1_init = 1,
  tau_2_init = 1,
  tau_3_init = 1,
  Beta_1_init = rep(0, ncol(X_1)),
  Beta_2_init = rep(0, ncol(X_2)),
  kappa_init = rep(0, ncol(S)),
  eta_init = rep(0, ncol(S)),
  zeta_1_init = rep(0, nrow(X_1)),
  zeta_2_init = rep(0, nrow(X_2)),
  a_zeta_1 = 1,
  b_zeta_1 = 1,
  a_zeta_2 = 1,
  b_zeta_2 = 1
) {
  N <- length(Z_1)
  p_1 <- ncol(X_1)
  p_2 <- ncol(X_2)
  r <- ncol(S)

  prior_prec_beta <- 1e-4
  prior_prec_tau <- 1 / 25
  m_safe <- pmax(m, 1e-8)

  d_inv <- Matrix::Diagonal(x = 1 / D)
  tS <- t(S)
  SD <- tS %*% d_inv
  M <- Matrix::forceSymmetric(tS %*% d_inv %*% S)
  tX1 <- t(X_1)
  tX2 <- t(X_2)
  XD1 <- Matrix::forceSymmetric(tX1 %*% d_inv %*% X_1)
  Z_m <- Z_2 - m_safe / 2

  tau_1 <- tau_1_init
  tau_2 <- tau_2_init
  tau_3 <- tau_3_init
  Beta_1 <- Beta_1_init
  Beta_2 <- Beta_2_init
  kappa <- kappa_init
  eta <- eta_init
  zeta_1 <- zeta_1_init
  zeta_2 <- zeta_2_init
  sigma2_zeta_1 <- 1
  sigma2_zeta_2 <- 1

  Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + zeta_1)
  Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2)

  keep <- nsim / nthin
  Beta_1_chain <- matrix(0, p_1, keep)
  Beta_2_chain <- matrix(0, p_2, keep)
  eta_chain <- matrix(0, r, keep)
  kappa_chain <- matrix(0, r, keep)
  zeta_1_chain <- matrix(0, N, keep)
  zeta_2_chain <- matrix(0, N, keep)
  Mu_1_chain <- matrix(0, N, keep)
  Mu_2_chain <- matrix(0, N, keep)
  tau_1_chain <- numeric(keep)
  tau_2_chain <- numeric(keep)
  tau_3_chain <- numeric(keep)
  sigma2_zeta_1_chain <- numeric(keep)
  sigma2_zeta_2_chain <- numeric(keep)

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)
  on.exit(close(pb), add = TRUE)

  for (iter in seq_len(nsim + nburn)) {
    w <- BayesLogit::rpg(N, m_safe, Mu_2)
    Omega <- Matrix::Diagonal(x = w)
    SO <- tS %*% Omega
    OM <- Matrix::forceSymmetric(SO %*% S)

    eta_prec <- Matrix::forceSymmetric(tau_1^2 * M + tau_2^2 * OM + Matrix::Diagonal(r))
    eta_rhs <-
      tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1 - zeta_1) +
      tau_2 * tS %*% (Z_m - Omega %*% (X_2 %*% Beta_2 + tau_3 * S %*% kappa + zeta_2))
    eta_var <- Matrix::solve(eta_prec)
    eta_mean <- Matrix::solve(eta_prec, eta_rhs)
    eta <- rmvn_once(as.vector(eta_mean), eta_var)

    beta1_prec <- Matrix::forceSymmetric(XD1 + prior_prec_beta * Matrix::Diagonal(p_1))
    beta1_rhs <- tX1 %*% d_inv %*% (Z_1 - tau_1 * S %*% eta - zeta_1)
    beta1_var <- Matrix::solve(beta1_prec)
    beta1_mean <- Matrix::solve(beta1_prec, beta1_rhs)
    Beta_1 <- rmvn_once(as.vector(beta1_mean), beta1_var)

    tau1_var <- 1 / as.numeric(t(eta) %*% M %*% eta + prior_prec_tau)
    tau1_mean <- tau1_var * as.numeric(t(eta) %*% SD %*% (Z_1 - X_1 %*% Beta_1 - zeta_1))
    tau_1 <- stats::rnorm(1, tau1_mean, sqrt(tau1_var))

    sigma2_zeta_1 <- 1 / stats::rgamma(
      1,
      shape = a_zeta_1 + N / 2,
      rate = b_zeta_1 + 0.5 * as.numeric(crossprod(zeta_1))
    )

    zeta1_prec_diag <- as.numeric(1 / D + 1 / sigma2_zeta_1)
    zeta1_var_diag <- 1 / zeta1_prec_diag
    zeta1_mean <- zeta1_var_diag * (1 / D) * (Z_1 - X_1 %*% Beta_1 - tau_1 * S %*% eta)
    zeta_1 <- stats::rnorm(N, zeta1_mean, sqrt(zeta1_var_diag))

    kappa_prec <- Matrix::forceSymmetric(tau_3^2 * OM + Matrix::Diagonal(r))
    kappa_rhs <- tau_3 * tS %*% (Z_m - Omega %*% (X_2 %*% Beta_2 + tau_2 * S %*% eta + zeta_2))
    kappa_var <- Matrix::solve(kappa_prec)
    kappa_mean <- Matrix::solve(kappa_prec, kappa_rhs)
    kappa <- rmvn_once(as.vector(kappa_mean), kappa_var)

    beta2_prec <- Matrix::forceSymmetric(tX2 %*% Omega %*% X_2 + prior_prec_beta * Matrix::Diagonal(p_2))
    beta2_rhs <- tX2 %*% (Z_m - Omega %*% (tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2))
    beta2_var <- Matrix::solve(beta2_prec)
    beta2_mean <- Matrix::solve(beta2_prec, beta2_rhs)
    Beta_2 <- rmvn_once(as.vector(beta2_mean), beta2_var)

    tau2_var <- 1 / as.numeric(t(eta) %*% OM %*% eta + prior_prec_tau)
    tau2_mean <- tau2_var * as.numeric(
      t(eta) %*% tS %*% (Z_m - Omega %*% (X_2 %*% Beta_2 + tau_3 * S %*% kappa + zeta_2))
    )
    tau_2 <- stats::rnorm(1, tau2_mean, sqrt(tau2_var))

    tau3_var <- 1 / as.numeric(t(kappa) %*% OM %*% kappa + prior_prec_tau)
    tau3_mean <- tau3_var * as.numeric(
      t(kappa) %*% tS %*% (Z_m - Omega %*% (X_2 %*% Beta_2 + tau_2 * S %*% eta + zeta_2))
    )
    tau_3 <- stats::rnorm(1, tau3_mean, sqrt(tau3_var))

    sigma2_zeta_2 <- 1 / stats::rgamma(
      1,
      shape = a_zeta_2 + N / 2,
      rate = b_zeta_2 + 0.5 * as.numeric(crossprod(zeta_2))
    )

    zeta2_prec_diag <- w + 1 / sigma2_zeta_2
    zeta2_var_diag <- 1 / zeta2_prec_diag
    zeta2_mean <- zeta2_var_diag * (Z_m - w * (X_2 %*% Beta_2 + tau_2 * S %*% eta + tau_3 * S %*% kappa))
    zeta_2 <- stats::rnorm(N, zeta2_mean, sqrt(zeta2_var_diag))

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + zeta_1)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2)

    utils::setTxtProgressBar(pb, iter)

    if (iter > nburn && (iter - nburn) %% nthin == 0) {
      pos <- (iter - nburn) / nthin
      Beta_1_chain[, pos] <- Beta_1
      Beta_2_chain[, pos] <- Beta_2
      eta_chain[, pos] <- eta
      kappa_chain[, pos] <- kappa
      zeta_1_chain[, pos] <- zeta_1
      zeta_2_chain[, pos] <- zeta_2
      Mu_1_chain[, pos] <- Mu_1
      Mu_2_chain[, pos] <- Mu_2
      tau_1_chain[pos] <- tau_1
      tau_2_chain[pos] <- tau_2
      tau_3_chain[pos] <- tau_3
      sigma2_zeta_1_chain[pos] <- sigma2_zeta_1
      sigma2_zeta_2_chain[pos] <- sigma2_zeta_2
    }
  }

  list(
    Beta_1.chain = Beta_1_chain,
    Beta_2.chain = Beta_2_chain,
    eta.chain = eta_chain,
    kappa.chain = kappa_chain,
    zeta_1.chain = zeta_1_chain,
    zeta_2.chain = zeta_2_chain,
    Mu_1.chain = Mu_1_chain,
    Mu_2.chain = Mu_2_chain,
    tau_1.chain = tau_1_chain,
    tau_2.chain = tau_2_chain,
    tau_3.chain = tau_3_chain,
    sigma2_zeta_1.chain = sigma2_zeta_1_chain,
    sigma2_zeta_2.chain = sigma2_zeta_2_chain
  )
}

SRE_sampler <- function(
  X, Y, D, S, nburn, nsim, nthin,
  zeta_init = rep(0, nrow(X)),
  a_zeta = 1,
  b_zeta = 1
) {
  N <- length(Y)
  p <- ncol(X)
  r <- ncol(S)

  prior_prec_beta <- 1e-4

  d_inv <- Matrix::Diagonal(x = 1 / D)
  tS <- t(S)
  SD <- tS %*% d_inv
  M <- Matrix::forceSymmetric(tS %*% d_inv %*% S)
  tX <- t(X)
  XD <- Matrix::forceSymmetric(tX %*% d_inv %*% X)

  Beta <- rep(0, p)
  U <- rep(0, r)
  zeta <- zeta_init
  sigma2_zeta <- 1
  Mu <- as.vector(X %*% Beta + S %*% U + zeta)

  keep <- nsim / nthin
  Beta_chain <- matrix(0, p, keep)
  U_chain <- matrix(0, r, keep)
  zeta_chain <- matrix(0, N, keep)
  Mu_chain <- matrix(0, N, keep)
  sigma2_zeta_chain <- numeric(keep)

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)
  on.exit(close(pb), add = TRUE)

  for (iter in seq_len(nsim + nburn)) {
    beta_prec <- Matrix::forceSymmetric(XD + prior_prec_beta * Matrix::Diagonal(p))
    beta_rhs <- tX %*% d_inv %*% (Y - S %*% U - zeta)
    beta_var <- Matrix::solve(beta_prec)
    beta_mean <- Matrix::solve(beta_prec, beta_rhs)
    Beta <- rmvn_once(as.vector(beta_mean), beta_var)

    u_prec <- Matrix::forceSymmetric(M + Matrix::Diagonal(r))
    u_rhs <- SD %*% (Y - X %*% Beta - zeta)
    u_var <- Matrix::solve(u_prec)
    u_mean <- Matrix::solve(u_prec, u_rhs)
    U <- rmvn_once(as.vector(u_mean), u_var)

    sigma2_zeta <- 1 / stats::rgamma(
      1,
      shape = a_zeta + N / 2,
      rate = b_zeta + 0.5 * as.numeric(crossprod(zeta))
    )

    zeta_prec_diag <- as.numeric(1 / D + 1 / sigma2_zeta)
    zeta_var_diag <- 1 / zeta_prec_diag
    zeta_mean <- zeta_var_diag * (1 / D) * (Y - X %*% Beta - S %*% U)
    zeta <- stats::rnorm(N, zeta_mean, sqrt(zeta_var_diag))

    Mu <- as.vector(X %*% Beta + S %*% U + zeta)

    utils::setTxtProgressBar(pb, iter)

    if (iter > nburn && (iter - nburn) %% nthin == 0) {
      pos <- (iter - nburn) / nthin
      Beta_chain[, pos] <- Beta
      U_chain[, pos] <- U
      zeta_chain[, pos] <- zeta
      Mu_chain[, pos] <- Mu
      sigma2_zeta_chain[pos] <- sigma2_zeta
    }
  }

  list(
    Beta.chain = Beta_chain,
    U.chain = U_chain,
    zeta.chain = zeta_chain,
    Mu.chain = Mu_chain,
    Sigma2_zeta.chain = sigma2_zeta_chain
  )
}

SRE_binomial <- function(
  X, Y, m, S, nburn, nsim, nthin,
  zeta_init = rep(0, nrow(X)),
  a_zeta = 1,
  b_zeta = 1
) {
  N <- length(Y)
  p <- ncol(X)
  r <- ncol(S)

  prior_prec_beta <- 1e-4
  m_safe <- pmax(m, 1e-8)
  Z_m <- Y - m_safe / 2
  tX <- t(X)
  tS <- t(S)

  Beta <- rep(0, p)
  U <- rep(0, r)
  zeta <- zeta_init
  sigma2_zeta <- 1
  Mu <- as.vector(X %*% Beta + S %*% U + zeta)

  keep <- nsim / nthin
  Beta_chain <- matrix(0, p, keep)
  U_chain <- matrix(0, r, keep)
  Mu_chain <- matrix(0, N, keep)
  zeta_chain <- matrix(0, N, keep)
  sigma2_zeta_chain <- numeric(keep)

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)
  on.exit(close(pb), add = TRUE)

  for (iter in seq_len(nsim + nburn)) {
    w <- BayesLogit::rpg(N, m_safe, Mu)
    Omega <- Matrix::Diagonal(x = w)
    SO <- tS %*% Omega
    OM <- Matrix::forceSymmetric(SO %*% S)

    beta_prec <- Matrix::forceSymmetric(tX %*% Omega %*% X + prior_prec_beta * Matrix::Diagonal(p))
    beta_rhs <- tX %*% (Z_m - Omega %*% (S %*% U + zeta))
    beta_var <- Matrix::solve(beta_prec)
    beta_mean <- Matrix::solve(beta_prec, beta_rhs)
    Beta <- rmvn_once(as.vector(beta_mean), beta_var)

    u_prec <- Matrix::forceSymmetric(OM + Matrix::Diagonal(r))
    u_rhs <- tS %*% (Z_m - Omega %*% (X %*% Beta + zeta))
    u_var <- Matrix::solve(u_prec)
    u_mean <- Matrix::solve(u_prec, u_rhs)
    U <- rmvn_once(as.vector(u_mean), u_var)

    sigma2_zeta <- 1 / stats::rgamma(
      1,
      shape = a_zeta + N / 2,
      rate = b_zeta + 0.5 * as.numeric(crossprod(zeta))
    )

    zeta_prec_diag <- w + 1 / sigma2_zeta
    zeta_var_diag <- 1 / zeta_prec_diag
    zeta_mean <- zeta_var_diag * (Z_m - w * (X %*% Beta + S %*% U))
    zeta <- stats::rnorm(N, zeta_mean, sqrt(zeta_var_diag))

    Mu <- as.vector(X %*% Beta + S %*% U + zeta)

    utils::setTxtProgressBar(pb, iter)

    if (iter > nburn && (iter - nburn) %% nthin == 0) {
      pos <- (iter - nburn) / nthin
      Beta_chain[, pos] <- Beta
      U_chain[, pos] <- U
      Mu_chain[, pos] <- Mu
      zeta_chain[, pos] <- zeta
      sigma2_zeta_chain[pos] <- sigma2_zeta
    }
  }

  list(
    Beta.chain = Beta_chain,
    U.chain = U_chain,
    Mu.chain = Mu_chain,
    zeta.chain = zeta_chain,
    Sigma2_zeta.chain = sigma2_zeta_chain
  )
}

draw_direct <- function(N, Mu1_true, Mu2_true, var_D, v) {
  mu1 <- stats::rnorm(N, Mu1_true, sqrt(var_D))
  mu2 <- stats::rnorm(N, Mu2_true, sqrt(v))
  list(mu1 = mu1, mu2 = mu2)
}

prep_design <- function(mu1, mu2, var_D, v) {
  eps_p <- 1e-6
  eps_v <- 1e-12
  eps_m <- 1e-8

  pi_hat <- plogis(mu2)
  pi_hat <- pmin(pmax(pi_hat, eps_p), 1 - eps_p)

  v_pi_hat <- (pi_hat * (1 - pi_hat))^2 * v
  v_pi_hat[!is.finite(v_pi_hat) | v_pi_hat <= 0] <- eps_v

  m <- pi_hat * (1 - pi_hat) / v_pi_hat
  m[!is.finite(m) | m <= 0] <- eps_m
  Z2 <- m * pi_hat

  mu1_mean <- mean(mu1)
  mu1_sd <- stats::sd(mu1)
  if (!is.finite(mu1_sd) || mu1_sd <= 0) {
    mu1_sd <- 1
  }

  Z1 <- (mu1 - mu1_mean) / mu1_sd
  D <- var_D / (mu1_sd^2)
  D[!is.finite(D) | D <= 0] <- eps_v

  list(
    Z1 = as.numeric(Z1),
    Z2 = as.numeric(Z2),
    D = as.numeric(D),
    m = as.numeric(m),
    pi_hat = as.numeric(pi_hat),
    mu1_mean = mu1_mean,
    mu1_sd = mu1_sd
  )
}

run_models <- function(X_1, X_2, Z1, Z2, D, m, S, nburn, nsim, nthin,
                       tau1 = 2, tau2 = -2, tau3 = 1) {
  fit_shared <- block_sampler_share_RE(
    X_1, X_2, Z1, Z2, D, m, S,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    tau_1_init = tau1,
    tau_2_init = tau2,
    tau_3_init = tau3
  )
  fit_gauss_uni <- SRE_sampler(X_1, Z1, D, S, nburn, nsim, nthin)
  fit_binom_uni <- SRE_binomial(X_2, Z2, m, S, nburn, nsim, nthin)

  list(
    shared = fit_shared,
    gauss_uni = fit_gauss_uni,
    binom_uni = fit_binom_uni
  )
}

extract_estimates <- function(fit_shared, fit_gauss_uni, fit_binom_uni, mu1_est = NULL) {
  q_gibbs_g <- apply(fit_shared$Mu_1.chain, 1, stats::quantile, c(0.025, 0.975))
  q_ind_g <- apply(fit_gauss_uni$Mu.chain, 1, stats::quantile, c(0.025, 0.975))
  q_gibbs_b <- apply(plogis(fit_shared$Mu_2.chain), 1, stats::quantile, c(0.025, 0.975))
  q_ind_b <- apply(plogis(fit_binom_uni$Mu.chain), 1, stats::quantile, c(0.025, 0.975))

  mu1_gibbs <- rowMeans(fit_shared$Mu_1.chain)
  mu1_ind <- rowMeans(fit_gauss_uni$Mu.chain)
  gibbs_g_l <- q_gibbs_g[1, ]
  gibbs_g_u <- q_gibbs_g[2, ]
  ind_g_l <- q_ind_g[1, ]
  ind_g_u <- q_ind_g[2, ]

  if (!is.null(mu1_est)) {
    sd1 <- stats::sd(mu1_est)
    mn1 <- mean(mu1_est)
    gibbs_g_l <- gibbs_g_l * sd1 + mn1
    gibbs_g_u <- gibbs_g_u * sd1 + mn1
    ind_g_l <- ind_g_l * sd1 + mn1
    ind_g_u <- ind_g_u * sd1 + mn1
    mu1_gibbs <- mu1_gibbs * sd1 + mn1
    mu1_ind <- mu1_ind * sd1 + mn1
  }

  list(
    gibbs_g_l = gibbs_g_l,
    gibbs_g_u = gibbs_g_u,
    ind_g_l = ind_g_l,
    ind_g_u = ind_g_u,
    gibbs_b_l = q_gibbs_b[1, ],
    gibbs_b_u = q_gibbs_b[2, ],
    ind_b_l = q_ind_b[1, ],
    ind_b_u = q_ind_b[2, ],
    mu1_gibbs = mu1_gibbs,
    mu1_ind = mu1_ind,
    p_gibbs = rowMeans(plogis(fit_shared$Mu_2.chain)),
    p_ind = rowMeans(plogis(fit_binom_uni$Mu.chain))
  )
}

interval_score <- function(lower, upper, truth, alpha = 0.05) {
  (upper - lower) +
    (2 / alpha) * pmax(lower - truth, 0) +
    (2 / alpha) * pmax(truth - upper, 0)
}

coverage_count <- function(lower, upper, truth) {
  rowSums((truth >= lower) & (truth <= upper))
}
