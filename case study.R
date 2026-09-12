
rweibull3 <- function(n, shape, scale, loc) {
  return(loc + rweibull(n, shape, scale))
}
evaluate_fit <- function(dat, pars, method_name) {
  if(any(is.na(pars))) {
    return(data.frame(Method = method_name, Shape = NA, Scale = NA, Loc = NA, 
                      KS = NA, AD = NA, MTTF = NA, B10 = NA))
  }
  
  shape <- pars[1]; scale <- pars[2]; loc <- pars[3]
  dat <- sort(dat)
  n <- length(dat)
  z <- (dat - loc) / scale
  z[z < 0] <- 0 
  p <- 1 - exp(-(z^shape))
  ecdf_steps <- (1:n) / n
  ecdf_steps_prev <- (0:(n-1)) / n
  ks_stat <- max(pmax(abs(p - ecdf_steps), abs(p - ecdf_steps_prev)))
  p_ad <- pmax(pmin(p, 1 - 1e-15), 1e-15) 
  p_ad <- sort(p_ad) 
  i <- 1:n
  ad_stat <- -n - (1/n) * sum((2*i - 1) * (log(p_ad) + log(1 - rev(p_ad))))
  mttf <- loc + scale * gamma(1 + 1/shape)
  b10 <- loc + scale * (-log(0.9))^(1/shape)
  return(data.frame(Method = method_name, 
                    Shape = round(shape, 4), 
                    Scale = round(scale, 4), 
                    Loc = round(loc, 4), 
                    KS = round(ks_stat, 4), 
                    AD = round(ad_stat, 4), 
                    MTTF = round(mttf, 4), 
                    B10 = round(b10, 4)))
}
#
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
#BL
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
#LSPF
library(pracma)
LSPF_Nagatsuka2013 <- function(dat) {
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
    }
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

run_case_study <- function(dat, dataset_name) {
  cat(sprintf("\n=======================================================\n"))
  cat(sprintf(" start %s (size n = %d)\n", dataset_name, length(dat)))
  cat(sprintf("=======================================================\n"))
  results <- list()

  # 1. BLE (Hall & Wang)
  res_ble <- tryCatch(BL_HallWang(dat), error=function(e) c(NA, NA, NA))
  results[[1]] <- evaluate_fit(dat, res_ble, "BLE")
  
  # 2. LSPF
  res_lspf <- tryCatch(LSPF_Nagatsuka2013(dat), error=function(e) c(NA, NA, NA))
  results[[2]] <- evaluate_fit(dat, as.numeric(res_lspf), "LSPF")
  
  # 3. UE
  res_ue <- tryCatch(NewU_Fast(dat), error=function(e) c(NA, NA, NA))
  results[[3]] <- evaluate_fit(dat, res_ue, "UE")
  
  # 4. BUE (Bootstrap)
  res_bue <- tryCatch(New_PBootstrap(dat, B=100), error=function(e) c(NA, NA, NA))
  results[[4]] <- evaluate_fit(dat, res_bue, "BUE")
  
  # 5. HUE
  res_hue <- tryCatch(NewHybrid_BL(dat), error=function(e) c(NA, NA, NA))
  results[[5]] <- evaluate_fit(dat, res_hue, "HUE")
  
  # 6. BHUE
  res_bhue <- tryCatch(NewHybrid_PBootstrap(dat, B=100), error=function(e) c(NA, NA, NA))
  results[[6]] <- evaluate_fit(dat, res_bhue, "BHUE")
  
  # joint
  final_df <- do.call(rbind, results)
  print(final_df)
  return(final_df)
}
# case

# case2(n=19)
data_fluid <- c(0.19, 0.78, 0.96, 1.31, 2.78, 3.16, 4.15, 4.67, 4.85, 6.50, 
                7.35, 8.01, 8.27, 12.06, 31.75, 32.52, 33.91, 36.71, 72.89)

res_case2 <- run_case_study(data_fluid, "Dielectric Breakdown of Insulating Fluid")
NewU_Fast(data_fluid)

#case1: Rockette et al. (1974)
dat_Rockette <- c(3.1, 4.6, 5.6, 6.8)
res_caseD <- run_case_study(dat_Rockette, "LSPF2:Rockette数据")

#case3: (n=128) -----------------

data_coetzee <- c(0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02,
                  0.03,0.04,0.06,0.08,0.10,0.10,0.12,0.12,0.12,0.13,
                  0.14,0.15,0.15,0.15,0.16,0.16,0.17,0.18,0.18,0.19,
                  0.20,0.21,0.22,0.23,0.25,0.26,0.28,0.28,0.30,0.32,
                  0.34,0.36,0.38,0.39,0.41,0.41,0.42,0.43,0.44,0.44,
                  0.45,0.45,0.50,0.53,0.56,0.58,0.58,0.61,0.62,0.62,
                  0.62,0.64,0.66,0.70,0.70,0.70,0.72,0.77,0.78,0.78,
                  0.80,0.82,0.83,0.85,0.86,0.96,0.97,0.98,0.99,1.05,
                  1.06,1.07,1.18,1.35,1.36,1.42,1.55,1.59,1.65,1.73,
                  1.77,1.79,1.80,1.91,2.09,2.14,2.15,2.15,2.31,2.33,
                  2.36,2.36,2.43,2.45,2.50,2.51,2.58,2.64,2.68,3.08,
                  3.94,4.12,4.33,4.42,4.53,4.88,4.97,5.11,5.32,5.55,
                  6.63,6.89,7.62,11.41,11.76,11.85,12.36,13.22)
res_coetzee <- run_case_study(data_coetzee, "Coetzee Equipment Degradation (n=128)")

