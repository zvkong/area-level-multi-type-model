# large_abs <- function(v, q) {
#   if (!is.numeric(q) || length(q) != 1L || q <= 0 || q > 1) {
#     stop("q must be a single number in (0, 1].")
#   }
#   n_keep <- max(1L, ceiling(q * length(v)))
#   order(abs(v), decreasing = TRUE)[seq_len(n_keep)]
# }
large_abs <- function(v, q) {
  pos_idx <- which(v > 0)

  if (length(pos_idx) == 0L) {
    stop("No positive eigenvalues found.")
  }

  pos_idx
}

cs_index <- function(eig_values, q) {
  if (!is.numeric(q) || length(q) != 1L || q <= 0 || q > 1) {
    stop("q must be a single number in (0, 1].")
  }
  score <- abs(eig_values)
  if (sum(score) <= 0) {
    return(list(r = 1L, index = 1L))
  }
  dt <- cbind(score, seq_along(eig_values), 0)
  dt <- dt[order(dt[, 1], decreasing = TRUE), , drop = FALSE]
  dt[, 3] <- cumsum(dt[, 1]) / sum(dt[, 1])
  k <- max(1L, sum(dt[, 3] <= q))
  index <- dt[seq_len(k), 2]
  list(r = length(index), index = index)
}

cs_index_p <- function(eig_values, q) {
  if (!is.numeric(q) || length(q) != 1L || q <= 0 || q > 1) {
    stop("q must be a single number in (0, 1].")
  }

  pos <- which(eig_values > 0)
  neg <- which(eig_values < 0)

  pos_v <- cbind(eig_values[pos]^2, pos, 0)
  neg_v <- cbind(eig_values[neg]^2, neg, 0)

  if (length(pos) > 0) {
    pos_v <- pos_v[order(pos_v[, 1], decreasing = TRUE), , drop = FALSE]
  }
  if (length(neg) > 0) {
    neg_v <- neg_v[order(neg_v[, 1], decreasing = TRUE), , drop = FALSE]
  }

  dt <- rbind(pos_v, neg_v)
  if (nrow(dt) == 0 || sum(dt[, 1]) <= 0) {
    return(list(r = 1L, index = 1L))
  }

  dt[, 3] <- cumsum(dt[, 1]) / sum(dt[, 1])
  k <- max(1L, sum(dt[, 3] <= q))
  index <- dt[seq_len(k), 2]
  list(r = length(index), index = index)
}

MSE <- function(real, prediction) {
  mean((real - prediction)^2)
}

MSE_lu <- function(real, prediction) {
  mean((exp(real) - exp(prediction))^2)
}

interval_score <- function(lower, upper, truth, alpha = 0.05) {
  (upper - lower) +
    (2 / alpha) * pmax(lower - truth, 0) +
    (2 / alpha) * pmax(truth - upper, 0)
}

coverage_count <- function(lower, upper, truth) {
  rowSums((truth >= lower) & (truth <= upper))
}

sample_from_prec <- function(prec, rhs) {
  prec <- Matrix::forceSymmetric(prec)
  U <- Matrix::chol(prec)
  mu <- Matrix::solve(U, Matrix::solve(t(U), rhs))
  z <- matrix(stats::rnorm(length(rhs)), ncol = 1)
  y <- Matrix::solve(U, z)
  as.vector(mu + y)
}

block_sampler_share_RE <- function(
  X_1, X_2, Z_1, Z_2, D, m, S, nburn, nsim, nthin,
  tau_1_init = 1, tau_2_init = 1, tau_3_init = 1,
  Beta_1_init = NULL, Beta_2_init = NULL,
  kappa_init = NULL, eta_init = NULL,
  zeta_1_init = NULL, zeta_2_init = NULL,
  a1 = 0.5, b1 = 0.5, a2 = 0.1, b2 = 0.1,
  beta_prec = 0.2
) {
  if (nsim %% nthin != 0) {
    stop("nsim must be divisible by nthin.")
  }

  N <- length(Z_1)
  p_1 <- ncol(X_1)
  p_2 <- ncol(X_2)
  r <- ncol(S)

  if (nrow(X_1) != N || nrow(X_2) != N || nrow(S) != N || length(D) != N || length(m) != N || length(Z_2) != N) {
    stop("Input dimensions do not match.")
  }

  if (is.null(Beta_1_init)) Beta_1_init <- rep(1, p_1)
  if (is.null(Beta_2_init)) Beta_2_init <- rep(1, p_2)
  if (is.null(kappa_init)) kappa_init <- rep(0, r)
  if (is.null(eta_init)) eta_init <- rep(0, r)
  if (is.null(zeta_1_init)) zeta_1_init <- rep(0, N)
  if (is.null(zeta_2_init)) zeta_2_init <- rep(0, N)

  # inv_sigma_tau <- 0.25
  inv_sigma_tau <- 0.01
  I_r <- Matrix::Diagonal(r)
  I_p1 <- Matrix::Diagonal(p_1)
  I_p2 <- Matrix::Diagonal(p_2)

  D_inv <- 1 / D
  M <- crossprod(S, D_inv * S)
  XD <- crossprod(X_1, D_inv * X_1)
  Z_m <- Z_2 - m / 2

  tau_1 <- tau_1_init
  tau_2 <- tau_2_init
  tau_3 <- tau_3_init
  Beta_1 <- Beta_1_init
  Beta_2 <- Beta_2_init
  kappa <- kappa_init
  eta <- eta_init
  zeta_1 <- zeta_1_init
  zeta_2 <- zeta_2_init

  sigma_zeta_1 <- 1
  sigma_zeta_2 <- 1

  Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + zeta_1)
  Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2)

  keep <- nsim %/% nthin
  tau_1.chain <- numeric(keep)
  tau_2.chain <- numeric(keep)
  tau_3.chain <- numeric(keep)
  Beta_1.chain <- matrix(0, nrow = p_1, ncol = keep)
  Beta_2.chain <- matrix(0, nrow = p_2, ncol = keep)
  Sigma_zeta_1.chain <- numeric(keep)
  Sigma_zeta_2.chain <- numeric(keep)
  kappa.chain <- matrix(0, nrow = r, ncol = keep)
  eta.chain <- matrix(0, nrow = r, ncol = keep)
  zeta_1.chain <- matrix(0, nrow = N, ncol = keep)
  zeta_2.chain <- matrix(0, nrow = N, ncol = keep)
  Mu_1.chain <- matrix(0, nrow = N, ncol = keep)
  Mu_2.chain <- matrix(0, nrow = N, ncol = keep)

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    w <- BayesLogit::rpg(N, pmax(m, 1e-8), Mu_2)
    OM <- crossprod(S, w * S)

    resid_1_no_eta <- Z_1 - X_1 %*% Beta_1 - zeta_1
    resid_2_no_eta <- Z_m - w * as.numeric(X_2 %*% Beta_2 + tau_3 * S %*% kappa + zeta_2)

    prec_eta <- Matrix::forceSymmetric(tau_1^2 * M + I_r + tau_2^2 * OM)
    rhs_eta <- tau_1 * crossprod(S, D_inv * resid_1_no_eta) +
      tau_2 * crossprod(S, resid_2_no_eta)
    eta <- sample_from_prec(prec_eta, rhs_eta)

    resid_1_no_beta <- Z_1 - tau_1 * S %*% eta - zeta_1
    prec_Beta_1 <- Matrix::forceSymmetric(XD + beta_prec * I_p1)
    rhs_Beta_1 <- crossprod(X_1, D_inv * resid_1_no_beta)
    Beta_1 <- sample_from_prec(prec_Beta_1, rhs_Beta_1)

    resid_1_no_tau <- Z_1 - X_1 %*% Beta_1 - zeta_1
    var_tau_1 <- as.numeric(1 / (crossprod(eta, M %*% eta) + inv_sigma_tau))
    mean_tau_1 <- as.numeric(var_tau_1 * crossprod(eta, crossprod(S, D_inv * resid_1_no_tau)))
    tau_1 <- stats::rnorm(1, mean_tau_1, sqrt(var_tau_1))

    prec_diag_z1 <- D_inv + 1 / sigma_zeta_1
    var_vec_z1 <- 1 / prec_diag_z1
    mean_vec_z1 <- var_vec_z1 * D_inv * as.numeric(Z_1 - X_1 %*% Beta_1 - tau_1 * S %*% eta)
    zeta_1 <- stats::rnorm(N, mean_vec_z1, sqrt(var_vec_z1))

    resid_2_no_tau2 <- Z_m - w * as.numeric(X_2 %*% Beta_2 + tau_3 * S %*% kappa + zeta_2)
    var_tau_2 <- as.numeric(1 / (crossprod(eta, OM %*% eta) + inv_sigma_tau))
    mean_tau_2 <- as.numeric(var_tau_2 * crossprod(eta, crossprod(S, resid_2_no_tau2)))
    tau_2 <- stats::rnorm(1, mean_tau_2, sqrt(var_tau_2))

    resid_2_no_tau3 <- Z_m - w * as.numeric(X_2 %*% Beta_2 + tau_2 * S %*% eta + zeta_2)
    var_tau_3 <- as.numeric(1 / (crossprod(kappa, OM %*% kappa) + inv_sigma_tau))
    mean_tau_3 <- as.numeric(var_tau_3 * crossprod(kappa, crossprod(S, resid_2_no_tau3)))
    tau_3 <- stats::rnorm(1, mean_tau_3, sqrt(var_tau_3))

    resid_2_no_kappa <- Z_m - w * as.numeric(X_2 %*% Beta_2 + tau_2 * S %*% eta + zeta_2)
    prec_kappa <- Matrix::forceSymmetric(tau_3^2 * OM + I_r)
    rhs_kappa <- tau_3 * crossprod(S, resid_2_no_kappa)
    kappa <- sample_from_prec(prec_kappa, rhs_kappa)

    resid_2_no_beta <- Z_m - w * as.numeric(tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2)
    prec_Beta_2 <- Matrix::forceSymmetric(crossprod(X_2, w * X_2) + beta_prec * I_p2)
    rhs_Beta_2 <- crossprod(X_2, resid_2_no_beta)
    Beta_2 <- sample_from_prec(prec_Beta_2, rhs_Beta_2)

    prec_diag_z2 <- w + 1 / sigma_zeta_2
    var_vec_z2 <- 1 / prec_diag_z2
    mean_vec_z2 <- var_vec_z2 * (
      Z_m - w * as.numeric(X_2 %*% Beta_2 + tau_2 * S %*% eta + tau_3 * S %*% kappa)
    )
    zeta_2 <- stats::rnorm(N, mean_vec_z2, sqrt(var_vec_z2))

    a_gaus <- a1 + N / 2
    b_gaus <- b1 + 0.5 * as.numeric(crossprod(zeta_1))
    sigma_zeta_1 <- 1 / stats::rgamma(1, shape = a_gaus, rate = b_gaus)

    a_bio <- a2 + N / 2
    b_bio <- b2 + 0.5 * as.numeric(crossprod(zeta_2))
    sigma_zeta_2 <- 1 / stats::rgamma(1, shape = a_bio, rate = b_bio)

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + zeta_1)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) %/% nthin
      tau_1.chain[pos] <- tau_1
      tau_2.chain[pos] <- tau_2
      tau_3.chain[pos] <- tau_3
      Beta_1.chain[, pos] <- Beta_1
      Beta_2.chain[, pos] <- Beta_2
      kappa.chain[, pos] <- kappa
      eta.chain[, pos] <- eta
      zeta_1.chain[, pos] <- zeta_1
      zeta_2.chain[, pos] <- zeta_2
      Sigma_zeta_1.chain[pos] <- sigma_zeta_1
      Sigma_zeta_2.chain[pos] <- sigma_zeta_2
      Mu_1.chain[, pos] <- Mu_1
      Mu_2.chain[, pos] <- Mu_2
    }
  }

  close(pb)

  list(
    Beta_1.chain = Beta_1.chain,
    Beta_2.chain = Beta_2.chain,
    kappa.chain = kappa.chain,
    eta.chain = eta.chain,
    zeta_1.chain = zeta_1.chain,
    zeta_2.chain = zeta_2.chain,
    Sigma_zeta_1.chain = Sigma_zeta_1.chain,
    Sigma_zeta_2.chain = Sigma_zeta_2.chain,
    Mu_1.chain = Mu_1.chain,
    Mu_2.chain = Mu_2.chain,
    tau_1.chain = tau_1.chain,
    tau_2.chain = tau_2.chain,
    tau_3.chain = tau_3.chain
  )
}

SRE_sampler <- function(
  X, Y, D, S, nburn, nsim, nthin,
  zeta.init = NULL,
  a = 0.1, b = 0.1,
  a.zeta = 0.5, b.zeta = 0.5,
  beta_prec = 0.2
) {
  if (nsim %% nthin != 0) {
    stop("nsim must be divisible by nthin.")
  }

  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(S)

  if (length(Y) != n || length(D) != n || nrow(S) != n) {
    stop("Input dimensions do not match.")
  }

  if (is.null(zeta.init)) zeta.init <- rep(0, n)

  Beta <- rep(1, p)
  Sigma2_u <- 1
  sigma_zeta <- 1
  U <- rep(0, r)
  Mu <- rep(0, n)
  zeta <- zeta.init

  I_r <- Matrix::Diagonal(r)
  I_p <- Matrix::Diagonal(p)
  D_inv <- 1 / D
  M <- crossprod(S, D_inv * S)
  XD <- crossprod(X, D_inv * X)

  keep <- nsim %/% nthin
  Beta.chain <- matrix(0, nrow = p, ncol = keep)
  Sigma2_u.chain <- numeric(keep)
  U.chain <- matrix(0, nrow = r, ncol = keep)
  Mu.chain <- matrix(0, nrow = n, ncol = keep)
  zeta.chain <- matrix(0, nrow = n, ncol = keep)
  Sigma_zeta.chain <- numeric(keep)

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    resid_no_beta <- Y - S %*% U - zeta
    prec_Beta <- Matrix::forceSymmetric(XD + beta_prec * I_p)
    rhs_Beta <- crossprod(X, D_inv * resid_no_beta)
    Beta <- sample_from_prec(prec_Beta, rhs_Beta)

    resid_no_u <- Y - X %*% Beta - zeta
    prec_U <- Matrix::forceSymmetric(M + I_r / Sigma2_u)
    rhs_U <- crossprod(S, D_inv * resid_no_u)
    U <- sample_from_prec(prec_U, rhs_U)

    prec_diag_z <- D_inv + 1 / sigma_zeta
    var_vec_z <- 1 / prec_diag_z
    mean_vec_z <- var_vec_z * D_inv * as.numeric(Y - X %*% Beta - S %*% U)
    zeta <- stats::rnorm(n, mean_vec_z, sqrt(var_vec_z))

    a_prime <- a + r / 2
    b_prime <- b + 0.5 * as.numeric(crossprod(U))
    Sigma2_u <- 1 / stats::rgamma(1, shape = a_prime, rate = b_prime)

    a_prime_zeta <- a.zeta + n / 2
    b_prime_zeta <- b.zeta + 0.5 * as.numeric(crossprod(zeta))
    sigma_zeta <- 1 / stats::rgamma(1, shape = a_prime_zeta, rate = b_prime_zeta)

    Mu <- as.vector(X %*% Beta + S %*% U + zeta)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) %/% nthin
      Beta.chain[, pos] <- Beta
      U.chain[, pos] <- U
      Sigma2_u.chain[pos] <- Sigma2_u
      Mu.chain[, pos] <- Mu
      zeta.chain[, pos] <- zeta
      Sigma_zeta.chain[pos] <- sigma_zeta
    }
  }

  close(pb)

  list(
    Beta.chain = Beta.chain,
    U.chain = U.chain,
    Sigma2_u.chain = Sigma2_u.chain,
    Mu.chain = Mu.chain,
    zeta.chain = zeta.chain,
    Sigma_zeta.chain = Sigma_zeta.chain
  )
}

SRE_binomial <- function(
  X, Y, m, S, nburn, nsim, nthin,
  zeta.init = NULL,
  a = 0.1, b = 0.1,
  a.zeta = 0.1, b.zeta = 0.1,
  beta_prec = 0.2
) {
  if (nsim %% nthin != 0) {
    stop("nsim must be divisible by nthin.")
  }

  N <- length(Y)
  p <- ncol(X)
  r <- ncol(S)

  if (nrow(X) != N || nrow(S) != N || length(m) != N) {
    stop("Input dimensions do not match.")
  }

  if (is.null(zeta.init)) zeta.init <- rep(0, N)

  Mu <- rep(0, N)
  Beta <- rep(1, p)
  zeta <- zeta.init
  sigma_zeta <- 1
  Sigma2_u <- 1
  U <- rep(0, r)

  I_r <- Matrix::Diagonal(r)
  I_p <- Matrix::Diagonal(p)

  keep <- nsim %/% nthin
  Beta.chain <- matrix(0, nrow = p, ncol = keep)
  Sigma2_u.chain <- numeric(keep)
  U.chain <- matrix(0, nrow = r, ncol = keep)
  Mu.chain <- matrix(0, nrow = N, ncol = keep)
  zeta.chain <- matrix(0, nrow = N, ncol = keep)
  Sigma_zeta.chain <- numeric(keep)

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    Z_m <- Y - m / 2
    w <- BayesLogit::rpg(N, pmax(m, 1e-8), Mu)
    OM <- crossprod(S, w * S)

    a_prime <- a + r / 2
    b_prime <- b + 0.5 * as.numeric(crossprod(U))
    Sigma2_u <- 1 / stats::rgamma(1, shape = a_prime, rate = b_prime)

    resid_no_beta <- Z_m - w * as.numeric(S %*% U + zeta)
    prec_Beta <- Matrix::forceSymmetric(crossprod(X, w * X) + beta_prec * I_p)
    rhs_Beta <- crossprod(X, resid_no_beta)
    Beta <- sample_from_prec(prec_Beta, rhs_Beta)

    resid_no_u <- Z_m - w * as.numeric(X %*% Beta + zeta)
    prec_U <- Matrix::forceSymmetric(OM + I_r / Sigma2_u)
    rhs_U <- crossprod(S, resid_no_u)
    U <- sample_from_prec(prec_U, rhs_U)

    prec_diag_z <- w + 1 / sigma_zeta
    var_vec_z <- 1 / prec_diag_z
    mean_vec_z <- var_vec_z * (
      Z_m - w * as.numeric(X %*% Beta + S %*% U)
    )
    zeta <- stats::rnorm(N, mean_vec_z, sqrt(var_vec_z))

    a_prime_zeta <- a.zeta + N / 2
    b_prime_zeta <- b.zeta + 0.5 * as.numeric(crossprod(zeta))
    sigma_zeta <- 1 / stats::rgamma(1, shape = a_prime_zeta, rate = b_prime_zeta)

    Mu <- as.vector(X %*% Beta + S %*% U + zeta)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) %/% nthin
      Beta.chain[, pos] <- Beta
      U.chain[, pos] <- U
      Sigma2_u.chain[pos] <- Sigma2_u
      Mu.chain[, pos] <- Mu
      zeta.chain[, pos] <- zeta
      Sigma_zeta.chain[pos] <- sigma_zeta
    }
  }

  close(pb)

  list(
    Beta.chain = Beta.chain,
    U.chain = U.chain,
    Sigma2_u.chain = Sigma2_u.chain,
    Mu.chain = Mu.chain,
    zeta.chain = zeta.chain,
    Sigma_zeta.chain = Sigma_zeta.chain
  )
}

draw_direct_independent <- function(Mu1_true, Mu2_true, var_D, v_logit) {
  n <- length(Mu1_true)
  s1 <- sqrt(pmax(var_D, 1e-12))
  s2 <- sqrt(pmax(v_logit, 1e-12))
  e1 <- s1 * stats::rnorm(n)
  e2 <- s2 * stats::rnorm(n)
  mu1 <- Mu1_true + e1
  mu2 <- pmin(pmax(Mu2_true + e2, -12), 12)
  list(mu1 = mu1, mu2 = mu2)
}

draw_direct_correlated <- function(Mu1_true, Mu2_true, var_D, v_logit, rho = 0.8) {
  n <- length(Mu1_true)
  mu1 <- numeric(n)
  mu2 <- numeric(n)

  for (i in seq_len(n)) {
    sd1 <- sqrt(pmax(var_D[i], 1e-12))
    sd2 <- sqrt(pmax(v_logit[i], 1e-12))
    cov12 <- rho * sd1 * sd2

    Sigma <- matrix(
      c(sd1^2, cov12,
        cov12, sd2^2),
      nrow = 2, byrow = TRUE
    )

    e <- MASS::mvrnorm(1, mu = c(0, 0), Sigma = Sigma)
    mu1[i] <- Mu1_true[i] + e[1]
    mu2[i] <- pmin(pmax(Mu2_true[i] + e[2], -12), 12)
  }

  list(mu1 = mu1, mu2 = mu2)
}

prep_design <- function(mu1, mu2, var_D, v) {
  eps_p <- 1e-6
  eps_v <- 1e-12
  eps_m <- 1e-8

  pi_hat <- pmin(pmax(plogis(mu2), eps_p), 1 - eps_p)
  v_pi_hat <- (pi_hat * (1 - pi_hat))^2 * v
  v_pi_hat[!is.finite(v_pi_hat) | v_pi_hat <= 0] <- eps_v

  m <- pi_hat * (1 - pi_hat) / v_pi_hat
  m[!is.finite(m) | m <= 0] <- eps_m

  Z2 <- m * pi_hat

  mu1_mean <- mean(mu1)
  mu1_sd <- stats::sd(mu1)
  if (!is.finite(mu1_sd) || mu1_sd <= 0) mu1_sd <- 1

  Z1 <- (mu1 - mu1_mean) / mu1_sd
  D <- var_D / (mu1_sd^2)
  D[!is.finite(D) | D <= 0] <- eps_v

  list(
    Z1 = Z1,
    Z2 = Z2,
    D = D,
    m = m,
    pi_hat = pi_hat,
    mu1_mean = mu1_mean,
    mu1_sd = mu1_sd
  )
}

run_models <- function(
  X_1, X_2, Z1, Z2, D, m, S, nburn, nsim, nthin,
  tau1 = 1, tau2 = 1, tau3 = 0,
  a_zeta_1 = 1, b_zeta_1 = 1,
  a_zeta_2 = 1, b_zeta_2 = 1,
  beta_prec = 0.2
) {
  cat("  -> [1/3] Fitting Multi-type Shared Model...\n")
  fit_shared <- block_sampler_share_RE(
    X_1, X_2, Z1, Z2, D, m, S,
    nburn = nburn, nsim = nsim, nthin = nthin,
    tau_1_init = tau1, tau_2_init = tau2, tau_3_init = tau3,
    a1 = a_zeta_1, b1 = b_zeta_1,
    a2 = a_zeta_2, b2 = b_zeta_2,
    beta_prec = beta_prec
  )

  cat("\n  -> [2/3] Fitting Univariate Gaussian Model...\n")
  fit_gauss_uni <- SRE_sampler(
    X_1, Z1, D, S, nburn, nsim, nthin,
    a.zeta = a_zeta_1, b.zeta = b_zeta_1,
    beta_prec = beta_prec
  )

  cat("\n  -> [3/3] Fitting Univariate Binomial Model...\n")
  fit_binom_uni <- SRE_binomial(
    X_2, Z2, m, S, nburn, nsim, nthin,
    a.zeta = a_zeta_2, b.zeta = b_zeta_2,
    beta_prec = beta_prec
  )

  list(
    shared = fit_shared,
    gauss_uni = fit_gauss_uni,
    binom_uni = fit_binom_uni
  )
}

extract_estimates <- function(
  fit_shared,
  fit_gauss_uni,
  fit_binom_uni,
  mu1_est = NULL,
  mu1_mean = NULL,
  mu1_sd = NULL
) {
  q_gibbs_g <- apply(fit_shared$Mu_1.chain, 1, stats::quantile, c(0.025, 0.975))
  q_ind_g <- apply(fit_gauss_uni$Mu.chain, 1, stats::quantile, c(0.025, 0.975))
  q_gibbs_b <- apply(plogis(fit_shared$Mu_2.chain), 1, stats::quantile, c(0.025, 0.975))
  q_ind_b <- apply(plogis(fit_binom_uni$Mu.chain), 1, stats::quantile, c(0.025, 0.975))

  if (is.null(mu1_mean) || is.null(mu1_sd)) {
    if (is.null(mu1_est)) {
      stop("Provide either mu1_est or both mu1_mean and mu1_sd.")
    }
    mu1_mean <- mean(mu1_est)
    mu1_sd <- stats::sd(mu1_est)
    if (!is.finite(mu1_sd) || mu1_sd <= 0) mu1_sd <- 1
  }

  list(
    gibbs_g_l = q_gibbs_g[1, ] * mu1_sd + mu1_mean,
    gibbs_g_u = q_gibbs_g[2, ] * mu1_sd + mu1_mean,
    ind_g_l = q_ind_g[1, ] * mu1_sd + mu1_mean,
    ind_g_u = q_ind_g[2, ] * mu1_sd + mu1_mean,
    gibbs_b_l = q_gibbs_b[1, ],
    gibbs_b_u = q_gibbs_b[2, ],
    ind_b_l = q_ind_b[1, ],
    ind_b_u = q_ind_b[2, ],
    mu1_gibbs = rowMeans(fit_shared$Mu_1.chain) * mu1_sd + mu1_mean,
    mu1_ind = rowMeans(fit_gauss_uni$Mu.chain) * mu1_sd + mu1_mean,
    p_gibbs = rowMeans(plogis(fit_shared$Mu_2.chain)),
    p_ind = rowMeans(plogis(fit_binom_uni$Mu.chain))
  )
}
fetch_acs_region <- function(variables, summary_var = NULL, geometry = FALSE, survey = "acs5") {
  map_dfr(
    states_use,
    ~ get_acs(
      geography = "tract",
      variables = variables,
      summary_var = summary_var,
      state = .x,
      geometry = geometry,
      year = year_use,
      survey = survey
    )
  )
}
make_facet_sf <- function(sf_obj, value_list, value_name) {
  labs <- names(value_list)

  build_one <- function(values, lab) {
    out <- sf_obj
    out[[value_name]] <- as.numeric(values)
    out$model_type <- lab
    out
  }

  dplyr::bind_rows(Map(build_one, value_list, labs))
}

plot_facets <- function(
  sf_obj,
  value_list,
  value_name,
  title,
  legend_lab = value_name,
  prob = 0.95,
  choose_pal = RColorBrewer::brewer.pal(9, "RdBu", rev = TRUE)
) {
  plot_df <- make_facet_sf(sf_obj, value_list, value_name)
  vals <- plot_df[[value_name]]
  vals <- vals[is.finite(vals)]
  vmin <- min(vals, na.rm = TRUE)
  vmax <- stats::quantile(vals, prob, na.rm = TRUE)

  ggplot2::ggplot(plot_df) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[value_name]]), colour = NA) +
    ggplot2::facet_wrap(~model_type) +
    ggplot2::scale_fill_gradientn(
      colours = choose_pal,
      limits = c(vmin, vmax),
      oob = scales::squish,
      name = legend_lab
    ) +
    ggplot2::labs(
      title = title,
      caption = "Data source: 2017–2021 ACS 5-year estimates, U.S. Census Bureau"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}