# utils_clean.R

large_abs <- function(v, q){
  av  <- abs(v)
  thr <- as.numeric(stats::quantile(av, q, names = FALSE, type = 7))
  which((v > 0) | (av > thr))
}

MSE <- function(real, prediction){
  mean((real - prediction)^2)
}

rmvn <- function(n, mu, Sigma) {
  p <- nrow(Sigma)
  L <- Matrix::chol(Matrix::forceSymmetric(Sigma))
  Z <- matrix(stats::rnorm(p * n), p, n)
  t(mu + L %*% Z)
}

block_sampler_share_RE <- function(
  X_1, X_2, Z_1, Z_2, D, m, S,
  nburn, nsim, nthin,
  tau_1_init = 1, tau_2_init = 1, tau_3_init = 1,
  Beta_1_init = rep(1, ncol(X_1)),
  Beta_2_init = rep(1, ncol(X_2)),
  kappa_init  = rep(0, ncol(S)),
  eta_init    = rep(0, ncol(S)),
  zeta_1_init = rep(0, nrow(X_1)),
  zeta_2_init = rep(0, nrow(X_2)),
  a1 = 0.5, b1 = 0.5, a2 = 0.1, b2 = 0.1
){
  N   <- length(Z_1)
  p_1 <- ncol(X_1); p_2 <- ncol(X_2); r <- ncol(S)

  inv_sigma_tau    <- 0.25
  inv_sigma_random <- Matrix::Diagonal(r)
  inv_sigma_zeta   <- Matrix::Diagonal(N)
  d  <- Matrix::Diagonal(x = 1 / D)
  tS <- t(S)
  SD <- tS %*% d
  M  <- Matrix::forceSymmetric(tS %*% d %*% S)
  tX1 <- t(X_1); tX2 <- t(X_2)
  XD  <- Matrix::forceSymmetric(tX1 %*% d %*% X_1)

  Z_m <- Z_2 - m/2
  SZ2 <- tS %*% Z_m

  tau_1 <- tau_1_init; tau_2 <- tau_2_init; tau_3 <- tau_3_init
  Beta_1 <- Beta_1_init; Beta_2 <- Beta_2_init
  kappa  <- kappa_init;  eta    <- eta_init
  zeta_1 <- zeta_1_init; zeta_2 <- zeta_2_init
  Mu_1 <- rep(1, N); Mu_2 <- rep(1, N)

  keep <- nsim / nthin
  Beta_1_chain <- array(0, c(p_1, keep))
  Beta_2_chain <- array(0, c(p_2, keep))
  kappa_chain  <- array(0, c(r, keep))
  eta_chain    <- array(0, c(r, keep))
  zeta_1_chain <- array(0, c(N, keep))
  zeta_2_chain <- array(0, c(N, keep))
  Mu_1_chain   <- array(0, c(N, keep))
  Mu_2_chain   <- array(0, c(N, keep))
  tau_1_chain  <- array(0, c(1, keep))
  tau_2_chain  <- array(0, c(1, keep))
  tau_3_chain  <- array(0, c(1, keep))

  cat("Starting ", keep, " iterations.\n", sep = "")
  pb <- utils::txtProgressBar(min = 0, max = (nsim + nburn) / nthin, style = 3)

  for (index in 1:(nsim + nburn)) {
    w <- BayesLogit::rpg(N, m, Mu_2)
    Omega <- Matrix::Diagonal(x = w)
    SO <- tS %*% Omega
    OM <- Matrix::forceSymmetric(SO %*% S)

    A_eta <- tau_1^2 * M + inv_sigma_random + tau_2^2 * OM
    b_eta <- tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1 - zeta_1) +
             tau_2 * SZ2 - tau_2 * SO %*% (X_2 %*% Beta_2 + tau_3 * S %*% kappa + zeta_2)
    V_eta  <- Matrix::solve(A_eta)
    m_eta  <- Matrix::solve(A_eta, b_eta)
    eta <- as.vector(rmvn(1, as.vector(m_eta), V_eta))

    Vb1 <- Matrix::solve(XD + 0.2 * Matrix::Diagonal(p_1))
    mb1 <- Matrix::solve(XD + 0.2 * Matrix::Diagonal(p_1),
                         tX1 %*% d %*% (Z_1 - tau_1 * S %*% eta - zeta_1))
    Beta_1 <- as.vector(MASS::mvrnorm(1, mu = mb1, Sigma = Vb1))

    v_t1 <- 1 / as.numeric(t(eta) %*% M %*% eta + inv_sigma_tau)
    m_t1 <- v_t1 * as.numeric(t(eta) %*% SD %*% (Z_1 - X_1 %*% Beta_1 - zeta_1))
    tau_1 <- stats::rnorm(1, m_t1, sqrt(v_t1))

    a_gaus <- a1 + N / 2
    b_gaus <- b1 + 0.5 * as.numeric(crossprod(zeta_1))
    sigma_zeta_1 <- 1 / stats::rgamma(1, shape = a_gaus, rate = b_gaus)

    a_bio <- a2 + N / 2
    b_bio <- b2 + 0.5 * as.numeric(crossprod(zeta_2))
    sigma_zeta_2 <- 1 / stats::rgamma(1, shape = a_bio, rate = b_bio)

    A_z1 <- d + inv_sigma_zeta / sigma_zeta_1
    v_z1 <- 1 / Matrix::diag(A_z1)
    m_z1 <- v_z1 * Matrix::diag(d) * (Z_1 - X_1 %*% Beta_1 - tau_1 * S %*% eta)
    zeta_1 <- stats::rnorm(N, m_z1, sqrt(v_z1))

    A_z2 <- Omega + inv_sigma_zeta / sigma_zeta_2
    v_z2 <- 1 / Matrix::diag(A_z2)
    m_z2 <- v_z2 * w * ((Z_m / w) - X_2 %*% Beta_2 - tau_2 * S %*% eta - tau_3 * S %*% kappa)
    zeta_2 <- stats::rnorm(N, m_z2, sqrt(v_z2))

    v_t2 <- 1 / as.numeric(t(eta) %*% OM %*% eta + inv_sigma_tau)
    m_t2 <- v_t2 * as.numeric(t(eta) %*% SZ2 - t(eta) %*% SO %*% (X_2 %*% Beta_2 + tau_3 * S %*% kappa + zeta_2))
    tau_2 <- stats::rnorm(1, m_t2, sqrt(v_t2))

    v_t3 <- 1 / as.numeric(t(kappa) %*% OM %*% kappa + inv_sigma_tau)
    m_t3 <- v_t3 * as.numeric(t(kappa) %*% SO %*% ((Z_m / w) - X_2 %*% Beta_2 - tau_2 * S %*% eta - zeta_2))
    tau_3 <- stats::rnorm(1, m_t3, sqrt(v_t3))

    A_k <- tau_3^2 * OM + inv_sigma_random
    b_k <- tau_3 * SO %*% ((Z_m / w) - X_2 %*% Beta_2 - tau_2 * S %*% eta - zeta_2)
    V_k <- Matrix::solve(A_k)
    m_k <- Matrix::solve(A_k, b_k)
    kappa <- as.vector(rmvn(1, as.vector(m_k), V_k))

    XtOX <- Matrix::forceSymmetric(tX2 %*% Omega %*% X_2)
    Vb2  <- Matrix::solve(XtOX + 0.2 * Matrix::Diagonal(p_2))
    mb2  <- Matrix::solve(XtOX + 0.2 * Matrix::Diagonal(p_2),
                          tX2 %*% Omega %*% ((Z_m / w) -
                            (tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2)))
    Beta_2 <- as.vector(MASS::mvrnorm(1, mu = mb2, Sigma = Vb2))

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta + zeta_1)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + tau_3 * S %*% kappa + zeta_2)

    utils::setTxtProgressBar(pb, index)
    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta_1_chain[, pos] <- Beta_1
      Beta_2_chain[, pos] <- Beta_2
      kappa_chain[, pos]  <- kappa
      eta_chain[, pos]    <- eta
      zeta_1_chain[, pos] <- zeta_1
      zeta_2_chain[, pos] <- zeta_2
      Mu_1_chain[, pos]   <- Mu_1
      Mu_2_chain[, pos]   <- Mu_2
      tau_1_chain[, pos]  <- tau_1
      tau_2_chain[, pos]  <- tau_2
      tau_3_chain[, pos]  <- tau_3
    }
  }

  list(
    Beta_1.chain = Beta_1_chain, Beta_2.chain = Beta_2_chain,
    kappa.chain = kappa_chain, eta.chain = eta_chain,
    zeta_1.chain = zeta_1_chain, zeta_2.chain = zeta_2_chain,
    Mu_1.chain = Mu_1_chain, Mu_2.chain = Mu_2_chain,
    tau_1.chain = tau_1_chain, tau_2.chain = tau_2_chain, tau_3.chain = tau_3_chain
  )
}

SRE_sampler <- function(
  X, Y, D, S, nburn, nsim, nthin,
  zeta.init = rep(0, nrow(X)),
  a = 0.1, b = 0.1, a.zeta = 0.5, b.zeta = 0.5
){
  m <- nrow(X); p <- ncol(X); r <- ncol(S)

  Beta <- rep(1, p)
  Sigma2_u <- 1
  Sigma2_zeta <- 1
  U <- rep(0, r)
  Mu <- rep(0, m)
  zeta <- zeta.init

  inv_sigma_zeta <- Matrix::Diagonal(m)
  inv_sigma_u    <- Matrix::Diagonal(r)
  d  <- Matrix::Diagonal(x = 1 / D)
  SD <- t(S) %*% d
  M  <- Matrix::forceSymmetric(SD %*% S)
  tX <- t(X)
  XD <- Matrix::forceSymmetric(tX %*% d %*% X)

  keep <- nsim / nthin
  Beta_chain <- array(0, c(p, keep))
  Sigma2_u_chain <- numeric(keep)
  U_chain <- array(0, c(r, keep))
  Mu_chain <- array(0, c(m, keep))
  zeta_chain <- array(0, c(m, keep))
  Sigma2_zeta_chain <- numeric(keep)

  cat("Starting ", keep, " iterations.\n", sep = "")
  pb <- utils::txtProgressBar(min = 0, max = (nsim + nburn) / nthin, style = 3)

  for (index in 1:(nsim + nburn)) {
    Vb <- Matrix::solve(XD)
    mb <- Matrix::solve(XD, tX %*% d %*% (Y - S %*% U - zeta))
    Beta <- as.vector(MASS::mvrnorm(1, mu = mb, Sigma = Vb))

    Au <- M + inv_sigma_u / Sigma2_u
    bu <- SD %*% (Y - X %*% Beta - zeta)
    Vu <- Matrix::solve(Au)
    mu_u <- Matrix::solve(Au, bu)
    U <- as.vector(rmvn(1, as.vector(mu_u), Vu))

    a_prime <- a + r / 2
    b_prime <- b + 0.5 * as.numeric(crossprod(U))
    Sigma2_u <- 1 / stats::rgamma(1, shape = a_prime, rate = b_prime)

    Az <- d + inv_sigma_zeta / Sigma2_zeta
    vz <- 1 / Matrix::diag(Az)
    mz <- vz * Matrix::diag(d) * (Y - X %*% Beta - S %*% U)
    zeta <- stats::rnorm(m, mz, sqrt(vz))

    a_prime_zeta <- a.zeta + m / 2
    b_prime_zeta <- b.zeta + 0.5 * as.numeric(crossprod(zeta))
    Sigma2_zeta <- 1 / stats::rgamma(1, shape = a_prime_zeta, rate = b_prime_zeta)

    Mu <- X %*% Beta + S %*% U + zeta

    utils::setTxtProgressBar(pb, index)
    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta_chain[, pos] <- Beta
      U_chain[, pos]    <- U
      Sigma2_u_chain[pos] <- Sigma2_u
      Mu_chain[, pos]   <- Mu
      zeta_chain[, pos] <- zeta
      Sigma2_zeta_chain[pos] <- Sigma2_zeta
    }
  }

  list(Beta.chain = Beta_chain, U.chain = U_chain, Sigma2_u.chain = Sigma2_u_chain,
       Mu.chain = Mu_chain, zeta.chain = zeta_chain, Sigma2_zeta.chain = Sigma2_zeta_chain)
}

SRE_binomial <- function(
  X, Y, m, S, nburn, nsim, nthin,
  zeta.init = rep(0, nrow(X)),
  a = 0.1, b = 0.1, a.zeta = 0.1, b.zeta = 0.1
){
  N <- length(Y); p <- ncol(X); r <- ncol(S)

  Mu <- rep(0, N)
  Beta <- rep(1, p)
  zeta <- zeta.init
  Sigma2_zeta <- 1
  Sigma2_u <- 1
  U <- rep(1, r)

  inv_sigma_random <- Matrix::Diagonal(r)
  inv_sigma_zeta   <- Matrix::Diagonal(N)
  tS <- t(S); tX <- t(X)

  keep <- nsim / nthin
  Beta_chain <- array(0, c(p, keep))
  Sigma2_u_chain <- numeric(keep)
  U_chain <- array(0, c(r, keep))
  Mu_chain <- array(0, c(N, keep))
  zeta_chain <- array(0, c(N, keep))
  Sigma2_zeta_chain <- numeric(keep)

  cat("Starting ", keep, " iterations.\n", sep = "")
  pb <- utils::txtProgressBar(min = 0, max = (nsim + nburn) / nthin, style = 3)

  for (index in 1:(nsim + nburn)) {
    w <- BayesLogit::rpg(N, m, Mu)
    Omega <- Matrix::Diagonal(x = w)
    SO <- tS %*% Omega
    OM <- Matrix::forceSymmetric(SO %*% S)

    a_prime <- a + r / 2
    b_prime <- b + 0.5 * as.numeric(crossprod(U))
    Sigma2_u <- 1 / stats::rgamma(1, shape = a_prime, rate = b_prime)

    XtOX <- Matrix::forceSymmetric(tX %*% Omega %*% X)
    Vb <- Matrix::solve(XtOX + 0.01 * Matrix::Diagonal(p))
    mb <- Matrix::solve(XtOX + 0.01 * Matrix::Diagonal(p),
                        tX %*% Omega %*% ((Y - m/2) / w - S %*% U - zeta))
    Beta <- as.vector(MASS::mvrnorm(1, mu = mb, Sigma = Vb))

    Au <- OM + inv_sigma_random / Sigma2_u
    bu <- SO %*% (((Y - m/2) / w) - X %*% Beta - zeta)
    Vu <- Matrix::solve(Au)
    mu_u <- Matrix::solve(Au, bu)
    U <- as.vector(rmvn(1, as.vector(mu_u), Vu))

    a_prime_zeta <- a.zeta + N / 2
    b_prime_zeta <- b.zeta + 0.5 * as.numeric(crossprod(zeta))
    Sigma2_zeta <- 1 / stats::rgamma(1, shape = a_prime_zeta, rate = b_prime_zeta)

    Az <- Omega + inv_sigma_zeta / Sigma2_zeta
    vz <- 1 / Matrix::diag(Az)
    mz <- vz * w * (((Y - m/2) / w) - X %*% Beta - S %*% U)
    zeta <- stats::rnorm(N, mz, sqrt(vz))

    Mu <- as.vector(X %*% Beta + S %*% U + zeta)

    utils::setTxtProgressBar(pb, index)
    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta_chain[, pos] <- Beta
      U_chain[, pos]    <- U
      Sigma2_u_chain[pos] <- Sigma2_u
      Mu_chain[, pos]   <- Mu
      zeta_chain[, pos] <- zeta
      Sigma2_zeta_chain[pos] <- Sigma2_zeta
    }
  }

  list(Beta.chain = Beta_chain, U.chain = U_chain, Sigma2_u.chain = Sigma2_u_chain,
       Mu.chain = Mu_chain, zeta.chain = zeta_chain, Sigma2_zeta.chain = Sigma2_zeta_chain)
}

draw_direct <- function(N, Mu1_true, Mu2_true, var_D, v) {
  mu1 <- stats::rnorm(N, Mu1_true, sqrt(var_D))
  mu2 <- stats::rnorm(N, Mu2_true, sqrt(v))
  list(mu1 = mu1, mu2 = mu2)
}

prep_design <- function(mu1, mu2, var_D, v) {
  pi_hat   <- plogis(mu2)
  v_pi_hat <- (pi_hat*(1 - pi_hat))^2 * v
  m        <- pi_hat*(1 - pi_hat) / pmax(v_pi_hat, 1e-12)
  Z2       <- m * pi_hat
  Z1       <- as.numeric(scale(mu1))
  D        <- var_D / stats::var(mu1)
  list(Z1 = Z1, Z2 = Z2, D = D, m = m, pi_hat = pi_hat)
}

run_models <- function(X_1, X_2, Z1, Z2, D, m, S, nburn, nsim, nthin,
                       tau1 = 2, tau2 = -2, tau3 = 0) {
  fit_shared     <- block_sampler_share_RE(X_1, X_2, Z1, Z2, D, m, S,
                                           nburn = nburn, nsim = nsim, nthin = nthin,
                                           tau_1_init = tau1, tau_2_init = tau2, tau_3_init = tau3)
  fit_gauss_uni  <- SRE_sampler(X_1, Z1, D, S, nburn, nsim, nthin)
  fit_binom_uni  <- SRE_binomial(X_2, Z2, m, S, nburn, nsim, nthin)
  list(shared = fit_shared, gauss_uni = fit_gauss_uni, binom_uni = fit_binom_uni)
}

extract_estimates <- function(fit_shared, fit_gauss_uni, fit_binom_uni, mu1_est) {
  q_gibbs_g <- apply(fit_shared$Mu_1.chain, 1, stats::quantile, c(0.025, 0.975))
  q_ind_g   <- apply(fit_gauss_uni$Mu.chain, 1, stats::quantile, c(0.025, 0.975))
  q_gibbs_b <- apply(plogis(fit_shared$Mu_2.chain), 1, stats::quantile, c(0.025, 0.975))
  q_ind_b   <- apply(plogis(fit_binom_uni$Mu.chain), 1, stats::quantile, c(0.025, 0.975))

  sd1  <- stats::sd(mu1_est); mn1 <- mean(mu1_est)
  gibbs_g_l <- q_gibbs_g[1,]*sd1 + mn1; gibbs_g_u <- q_gibbs_g[2,]*sd1 + mn1
  ind_g_l   <- q_ind_g[1,]*sd1   + mn1; ind_g_u   <- q_ind_g[2,]*sd1   + mn1

  mu1_gibbs <- rowMeans(fit_shared$Mu_1.chain)*sd1 + mn1
  mu1_ind   <- rowMeans(fit_gauss_uni$Mu.chain)   *sd1 + mn1
  p_gibbs   <- rowMeans(plogis(fit_shared$Mu_2.chain))
  p_ind     <- rowMeans(plogis(fit_binom_uni$Mu.chain))

  list(
    gibbs_g_l = gibbs_g_l, gibbs_g_u = gibbs_g_u,
    ind_g_l = ind_g_l,     ind_g_u   = ind_g_u,
    gibbs_b_l = q_gibbs_b[1,], gibbs_b_u = q_gibbs_b[2,],
    ind_b_l   = q_ind_b[1,],   ind_b_u   = q_ind_b[2,],
    mu1_gibbs = mu1_gibbs, mu1_ind = mu1_ind,
    p_gibbs = p_gibbs, p_ind = p_ind
  )
}

interval_score <- function(lower, upper, truth, alpha = 0.05) {
  w <- 2/(alpha/2)
  (upper - lower) + w*pmax(lower - truth, 0) + w*pmax(truth - upper, 0)
}

coverage_count <- function(lower, upper, truth) {
  rowSums((truth >= lower) & (truth <= upper))
}
