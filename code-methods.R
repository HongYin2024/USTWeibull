# UE 
NewU_Fast <- function(dat){
  n <- length(dat)
  dat <- sort(dat) 
  w <- (n - 1:n)
  N_pairs <- n * (n - 1) / 2
  sum_min <- sum(w * dat)
  sum_min_sq <- sum(w * (dat^2))
  mu_m <- sum_min / N_pairs
  mu_m2 <- sum_min_sq / N_pairs
  var_min_est <- mu_m2 - (mu_m)^2
  var_dat <- var(dat)
  ratio <- var_dat / var_min_est
  if(is.na(ratio) || ratio <= 1.000001) return(c(NA, NA, NA))
  shape_est <- (2 * log(2)) / log(ratio)
  if(is.na(shape_est) || shape_est <= 0) return(c(NA, NA, NA))
  g1 <- gamma(1 + 1/shape_est)
  g2 <- gamma(1 + 2/shape_est)
  denom <- g2 - g1^2
  if(is.na(denom) || denom <= 0) return(c(shape_est, NA, NA))
  scale_est <- sqrt(var_dat / denom)
  loc_est <- dat[1] - (n^(-1/shape_est)) * scale_est * g1
  return(c(shape_est, scale_est, loc_est))
}

#BLE(Hall & Wang 2005)
BL_HallWang <- function(dat, par0=NULL){
  n <- length(dat); x1 <- dat[1]; x2 <- dat[2]
  nll_BL <- function(par){
    shape <- par[1]; scale <- par[2]; theta <- par[3]
    if(shape <= 1e-5 || scale <= 1e-5 || theta >= x1 - 1e-9)
      return(1e10)
    xt <- dat - theta
    if(any(xt <= 0)) return(1e10)
    logL <- n*log(shape) - n*shape*log(scale) +
      (shape-1)*sum(log(xt)) - sum((xt/scale)^shape)
    penalty <- log(x1 - theta) - log(x2 - theta)
    return(-(logL + penalty))
  }
  
  if(is.null(par0)) par0 <- c(1, sd(dat), x1 - 0.1*sd(dat))
  res <- try(optim(par=par0, fn=nll_BL, method="Nelder-Mead",
                   control=list(maxit=2000)), silent=TRUE)
  if(inherits(res, "try-error")) return(c(NA, NA, NA))
  return(res$par)
}

#LSPF(Nagatsuka et al. 2013 CSDA) ---

LSPF_Nagatsuka2013 <- function(dat){
  n <- length(dat)
  dat <- sort(dat)
  x1 <- dat[1]
  xn <- dat[n]
  w <- (dat - x1) / (xn - x1)
  neg_log_lik_lspf <- function(beta_val) {
    if(beta_val <= 0.001) return(1e10)
    integrand <- function(y_vec) {
      val <- numeric(length(y_vec))
      for(j in seq_along(y_vec)) {
        y <- y_vec[j]
        if(y == 0) y <- 1e-12 
        M_y <- y + 1 
        term1 <- (beta_val - 1) * sum(log(y + w))
        sum_pow <- sum( ((y + w) / M_y)^beta_val )
        term2 <- n * beta_val * log(M_y) + n * log(sum_pow)
        val[j] <- exp(term1 - term2)
      }
      return(val)
    res <- try(integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-6, stop.on.error = FALSE)$value, silent = TRUE)
    if(inherits(res, "try-error") || is.na(res) || res <= 0) return(1e10)
    return(-( (n-1)*log(beta_val) + log(res) ))
  }
  opt <- try(optimize(neg_log_lik_lspf, interval = c(0.1, 50), tol = 1e-10), silent = TRUE)
  if(inherits(opt, "try-error")) return(c(NA, NA, NA))
  beta_hat <- opt$minimum
  alpha_init <- ( mean( (dat - x1)^beta_hat ) ) ^ (1 / beta_hat)
  gamma_hat <- x1 - (n^(-1/beta_hat)) * alpha_init * gamma(1 + 1/beta_hat)
  alpha_hat <- ( mean( (dat - gamma_hat)^beta_hat ) ) ^ (1 / beta_hat)
  return(c(beta = beta_hat, alpha = alpha_hat, gamma = gamma_hat))
  }
  
#BUE
New_PBootstrap <- function(dat, B=100){
  n <- length(dat)
  est_orig <- NewU_Fast(dat)
  if(any(is.na(est_orig)) || est_orig[1] <= 0) return(c(NA, NA, NA))
  
  boot_ests <- matrix(NA, nrow=B, ncol=3)
  for(b in 1:B){
    sim_dat <- sort(rweibull3(n, est_orig[1], est_orig[2], est_orig[3]))
    res_b <- NewU_Fast(sim_dat)
    if(!any(is.na(res_b))) boot_ests[b, ] <- res_b
  }
  valid_b <- complete.cases(boot_ests)
  if(sum(valid_b) < 10) return(est_orig)
  
  est_mean <- colMeans(boot_ests[valid_b, , drop=FALSE])
  est_bc <- 2 * est_orig - est_mean
  if(est_bc[1] <= 1e-3) est_bc[1] <- est_orig[1]
  if(est_bc[2] <= 1e-3) est_bc[2] <- est_orig[2]
  if(est_bc[3] >= dat[1]) est_bc[3] <- dat[1] - 1e-5
  return(est_bc)
}
#HUE
NewHybrid_BL <- function(dat){
  u_res <- NewU_Fast(dat)
  shape_fixed <- u_res[1]
  if(is.na(shape_fixed) || shape_fixed <= 0) return(c(NA, NA, NA))
  
  n <- length(dat); x1 <- dat[1]; x2 <- dat[2]
  nll_hybrid_bl <- function(par){
    scale <- par[1]; theta <- par[2]
    if(scale <= 1e-5 || theta >= x1 - 1e-9) return(1e10)
    xt <- dat - theta
    if(any(xt <= 0)) return(1e10)
    penalty <- log(x1 - theta) - log(x2 - theta)
    logL <- n*log(shape_fixed) - n*shape_fixed*log(scale) + 
      (shape_fixed-1)*sum(log(xt)) - sum((xt/scale)^shape_fixed)
    return(-(logL + penalty))
  }
  init_par <- c(u_res[2], u_res[3])
  if(is.na(init_par[1])) init_par <- c(sd(dat), x1 - 0.1*sd(dat))
  
  res <- try(optim(par=init_par, fn=nll_hybrid_bl, method="Nelder-Mead"), silent=TRUE)
  if(inherits(res, "try-error")) return(c(shape_fixed, NA, NA))
  return(c(shape_fixed, res$par[1], res$par[2]))
}

#BHUE

NewHybrid_PBootstrap <- function(dat, B=100){
  n <- length(dat)
  est_orig <- NewHybrid_BL(dat)
  if(any(is.na(est_orig)) || est_orig[1] <= 0) return(c(NA, NA, NA))
  boot_ests <- matrix(NA, nrow=B, ncol=3)
  for(b in 1:B){
    sim_dat <- sort(rweibull3(n, est_orig[1], est_orig[2], est_orig[3]))
    res_b <- NewHybrid_BL(sim_dat)
    if(!any(is.na(res_b))) boot_ests[b, ] <- res_b
  }
  valid_b <- complete.cases(boot_ests)
  if(sum(valid_b) < 10) return(est_orig)
  est_mean <- colMeans(boot_ests[valid_b, , drop=FALSE])
  est_bc <- 2 * est_orig - est_mean
  
  if(est_bc[1] <= 1e-3) est_bc[1] <- est_orig[1]
  if(est_bc[2] <= 1e-3) est_bc[2] <- est_orig[2]
  if(est_bc[3] >= dat[1]) est_bc[3] <- dat[1] - 1e-5
  return(est_bc)
  }
