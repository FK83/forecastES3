t_test_mat <- function(dat, lag_length = 0, c = 0, nm1 = FALSE, eps = 1e-6){
  
  # Dimensions
  n <- nrow(dat)
  k <- ncol(dat)
  
  # Column means
  m <- colMeans(dat)

  # De-mean data 
  dat_dm <- scale(dat, scale = FALSE)
  
  # Columnwise HAC 
  v <- vHACC_mat(x = dat_dm, k = lag_length)
  
  # Option to divide by n-1 rather than n (as done by sandwich package)
  if (nm1){
    v <- v * (n/(n-1))
  }
  
  # Replace very small variances by NAs
  v[v < eps] <- NA
  
  # t-test stat for all data columns
  t_test <- sqrt(n)*(m-c)/sqrt(v)
  
  return(t_test)
  
}


#' @importFrom RcppRoll roll_sum
# Goetze/Kuensch Bootstrap procedure for simple t-test
# (simultaneously across multiple columns of a data matrix)
gk_boot <- function(dat, l, n_boot = 5000){
  
  # Enforce matrix input
  if (!is.matrix(dat)){
    dat <- matrix(dat)
  }
  
  # Dimensions
  n <- nrow(dat)
  k <- ncol(dat)
  n_blocks <- n/l
  
  if (n_blocks != as.integer(n_blocks)){
    stop("code assumes that sample size is a multiple of block length")
  }
  
  # Compute mean and second moment for each of the n-l+ overlapping blocks
  # (values will be looked up later during bootstrap iterations)
  # Variance formula exploits that obs are independent across blocks
  # See G/K, p.1918 -> Subtract squared overall mean later
  # (corresponds to HAC with rectangular kernel)
  m <- m2 <- matrix(NA, n-l+1, k)
  for (dd in 1:k){
    tmp <- roll_sum(dat[, dd], n = l, align = "left")
    m[, dd] <- tmp/l
    m2[, dd] <- (tmp^2)/l
  }
  
  # Compute bootstrap expectations for centering
  # See G/K, p.1917
  cent <- rep(0, k)
  for (jj in 1:n){
    c <- min(jj/l, 1, (n + 1 - jj)/l)
    cent <- cent + c*dat[jj, ]
  }
  cent <- cent / (n - l + 1)
  
  # Bootstrap t stats will be saved here
  t_boot <- matrix(NA, n_boot, k)
  
  # Loop over bootstrap iterations
  for (bb in 1:n_boot){
    
    # Draw starting points of blocks 
    sp <- sample.int(n-l+1, n_blocks, replace = TRUE)
    
    # Mean of current bootstrap sample in each column
    mn <- colMeans(m[sp, , drop = FALSE])
    
    # Numerator of t stat for each column
    # (mean of bootstrap sample minus bootstrap expectation)
    t_num <- mn - cent
    
    # Denominator of t stat for each column
    # (root of variance of bootstrap sample)
    t_denom <- sqrt(colMeans(m2[sp, , drop = FALSE]) - mn^2)
    
    # Record t stat for each column
    # (need to multiply by root of sample size, see G/K, p.1918)
    t_boot[bb, ] <- (sqrt(n)*t_num)/t_denom
    
    # Print progress
    if (bb %% 1000 == 0) print(paste("Now at iteration", bb))
    
  }
  
  # Compute t-stats for original sample (using standard HAC estimator)
  t_orig <- t_test_mat(dat = dat, lag_length = l - 1)
  
  # Transform to matrix 
  t_orig_mat <- t(matrix(t_orig, ncol(dat), n_boot))
  
  # P-values for two-sided test
  pv_two_sided <- colMeans(abs(t_orig_mat) > abs(t_boot))
  
  # P-values for one-sided test (smaller or equal)
  pv_one_sided <- colMeans(t_orig_mat > t_boot)
  
  return(list(pv_two_sided = pv_two_sided, pv_one_sided = pv_one_sided,
              t_orig = t_orig))
  
}

# Naive procedure for t_test
naive_boot <- function(dat, l, lag_length = l-1, n_boot = 5000){
  
  # Enforce matrix input
  if (!is.matrix(dat)){
    dat <- matrix(dat)
  }
  
  # Dimensions
  n <- nrow(dat)
  k <- ncol(dat)
  n_blocks <- n/l
  
  if (n_blocks != as.integer(n_blocks)){
    stop("code assumes that sample size is a multiple of block length")
  }
  
  # Column means (needed for centering)
  m <- colMeans(dat)
  
  # Bootstrap t stats will be saved here
  t_boot <- matrix(NA, n_boot, k)
  
  # Loop over bootstrap iterations
  for (bb in 1:n_boot){
    
    # Draw starting points of blocks 
    sp <- sample.int(n-l+1, n_blocks, replace = TRUE)
    
    # Construct bootstrap indices and bootstrap sample
    ind <- unlist(lapply(sp, function(z) seq(z, length.out = l)))
    dat_boot <- dat[ind, ]
    
    # Record t stat for each column
    # (need to multiply by root of sample size, see G/K, p.1918)
    t_boot[bb, ] <- t_test_mat(dat_boot, lag_length = lag_length, 
                               c = m)
    
    # Print progress
    if (bb %% 1000 == 0) print(paste("Now at iteration", bb))
    
  }
  
  # Compute t-stats for original sample (using standard HAC estimator)
  t_orig <- t_test_mat(dat = dat, lag_length = lag_length)
  
  # Transform to matrix 
  t_orig_mat <- t(matrix(t_orig, ncol(dat), n_boot))
  
  # P-values for two-sided test
  pv_two_sided <- colMeans(abs(t_orig_mat) > abs(t_boot))
  
  # P-values for one-sided test (smaller or equal)
  pv_one_sided <- colMeans(t_orig_mat > t_boot)
  
  return(list(pv_two_sided = pv_two_sided, pv_one_sided = pv_one_sided, t_orig = t_orig))
  
}

#' Simulate many AR(1) processes with zero mean
#' 
#' @importFrom stats arima.sim
#' @param a AR parameter
#' @param n_sample Sample size
#' @param n_mc Nr of MC iterations
#' @return matrix, each col is a simulated AR process
ar_sim <- function(a = 0.5, n_sample = 500, n_mc = 1000){
  
  out <- matrix(0, n_sample, n_mc)
  
  for (jj in 1:n_mc){
    
    out[, jj] <- arima.sim(list(ar = a), n = n_sample)
    
  }
  
  out
  
}

#' T-tests using stationary bootstrap
#' 
#' @importFrom matrixStats colMeans2 rowMaxs
#' @importFrom boot tsboot
#' @param dat Matrix, each column is a time series
#' @param n_boot Nr of bootstrap iterations
#' @param block_length Expected block length of stationary bootstrap
#' @param eps Threshold for setting t-stat to NA
#' @return t-stats 
stationary_boot <- function(dat, n_boot, block_length, eps = 1e-6){
  
  # Sample size
  n <- nrow(dat)
  
  # Parameter for geometric dist
  p <- 1/block_length
  
  # Column means of data
  m <- colMeans2(dat)
  
  # Column-wise bootstrap standard deviation of x
  # (Politis and Romano, p.1305)
  s <- sqrt(vSB_mat(dat, p))
  s[s^2 < eps] <- NA
  
  # Draw bootstrap indices
  boot_ind <- t(tsboot(1:n, identity, R = n_boot, l = block_length, 
                     sim = "geom")$t)
  
  # Enter t-stats here
  t_boot <- matrix(0, n_boot, ncol(dat))
  
  # Bootstrap loop
  for (jj in 1:n_boot){
    
    # Column means for this iteration
    cM <- colMeans2(dat, rows = boot_ind[, jj]) 
    
    # T-stats for this iteration
    t_boot[jj, ] <- sqrt(n)*(cM - m)/s
    
  }
  
  # Compute t-stats for original sample
  # (using SB variance estimator, as in Hansen 2005)
  t_orig <- sqrt(n)*m/s
  
  # Transform to matrix 
  t_orig_mat <- t(matrix(t_orig, ncol(dat), n_boot))
  
  # P-values for two-sided test
  pv_two_sided <- colMeans(abs(t_orig_mat) > abs(t_boot))
  
  # P-values for one-sided test (smaller or equal)
  pv_one_sided <- colMeans(t_orig_mat > t_boot)
  
  return(list(pv_two_sided = pv_two_sided, pv_one_sided = pv_one_sided,
              t_orig = t_orig, t_boot = t_boot))
  
}
