#' @useDynLib forecastES3
#' @importFrom Rcpp sourceCpp
NULL
#' @importFrom sandwich NeweyWest
#' @importFrom VaRES varT esT
#' @importFrom magrittr %>%
#' @importFrom graphics abline legend lines matplot plot polygon par axis
#' @importFrom stats coefficients lm pnorm qnorm rnorm rt sd
#' @importFrom utils head tail

# Function for first extremal score
S1 <- function(x1, y, v1, alpha = 0.05, modified = FALSE){

  if (length(x1) != length(y) | length(alpha) > 1){
  
    stop("Invalid input..")
	
  }

  if (length(v1) > 1){
	n <- length(x1)
	x1 <- matrix(x1, nrow = n, ncol = length(v1))
	y <- matrix(y, nrow = n, ncol = length(v1))
	v1 <- t(matrix(v1, nrow = length(v1), ncol = n))
  }

  if (modified){
    out <- ((y <= x1) - alpha) * (v1 <= x1) - (y <= x1) * (v1 <= y) 
  } else {
    out <- ((y <= x1) - alpha) * ((v1 <= x1) - (v1 <= y))
  }
  out
}

# Function for second extremal score
S2 <- function(x1, x2, y, v2, alpha = 0.05, modified = FALSE){

  if (length(x1) != length(y) | length(x2) != length(y) | length(alpha) > 1){
  
    stop("Invalid input..")
	
  }

  if (length(v2) > 1){
	n <- length(x1)
	x1 <- matrix(x1, nrow = n, ncol = length(v2))
	x2 <- matrix(x2, nrow = n, ncol = length(v2))
	y <- matrix(y, nrow = n, ncol = length(v2))
	v2 <- t(matrix(v2, nrow = length(v2), ncol = n))
  }

  if (modified){
    out1 <- (v2 <= x2) * ( (y <= x1) * (x1 - y)/alpha - (x1 - v2) )  
    out2 <- 0
  } else {
    out1 <- (v2 <= x2) * ( (y <= x1) * (x1 - y)/alpha - (x1 - v2) )  
    out2 <- (v2 <= y)*(y - v2)
  }
  (out1 + out2)
}

get_grid <- function(TSfc, plot_nr){
  # Determine grid
  if (plot_nr == 1){
    xax <- TSfc[, c(1, 3, 5)] %>% as.numeric %>% 
      unique %>% sort
  } else {
    xax <- TSfc[, c(2, 4)] %>% as.numeric %>% 
      unique %>% sort
  }  
  xax
}



#' Murphy Diagrams for VaR and Expected Shortfall
#' @param TSfc Matrix containing forecasts and realizations. First two columns contain VaR/ES forecasts from method A; cols 3-4 contain VaR/ES forecasts from method B; col 5 contains realizations.
#' @param alpha scalar, level of VaR/ES
#' @param plot_nr Either 1 (first elementary score) or 2 (second elementary score)
#' @param labels Vector of labels (set to NULL to omit labels)
#' @param grid Either NULL (to use grid at all data points) or numeric vector of customized grid points.
#' @param legend_loc Location of legend
#' @param murphy_cols Colors to be used for first plot
#' @param cex_gen Font size
#' @param ... other plotting parameters
#' @rdname murphy
#' @export
murphy_VaRES <- function(TSfc, alpha = 0.025, plot_nr = 2, labels = NULL, 
                         grid = NULL, legend_loc = NULL, murphy_cols = NULL,
                         cex_gen = 1.6, ...){
  
  # Get grid
  if (is.null(grid)){
    xax <- get_grid(TSfc, plot_nr)  
  } else {
    xax <- grid
  }
  
  # Enter data here
  dat <- matrix(0, length(xax), 2)
    
  # Compute scores, depending on which plot is produced
  if (plot_nr == 1){
    dat[, 1] <- colMeans(S1(TSfc[, 1], TSfc[, 5], xax, alpha = alpha, modified = FALSE))
    dat[, 2] <- colMeans(S1(TSfc[, 3], TSfc[, 5], xax, alpha = alpha, modified = FALSE))
  } else if (plot_nr == 2){
    dat[, 1] <- colMeans(S2(TSfc[, 1], TSfc[, 2], TSfc[, 5], xax, alpha = alpha, modified = FALSE))
    dat[, 2] <- colMeans(S2(TSfc[, 3], TSfc[, 4], TSfc[, 5], xax, alpha = alpha, modified = FALSE))
  }

  # Cosmetics
  if (plot_nr == 1){
    xlab <- expression(v[1])
  } else {
    xlab <- expression(v[2])
  }
  if (is.null(murphy_cols)) murphy_cols <- c("tomato1", "steelblue4")
  if (is.null(labels)) labels <- c("Method 1", "Method 2")
  
  # Make plot
  matplot(x = xax, y = dat, type = "l", bty = "n", col = murphy_cols, 
          xlab = xlab, lwd = 3, ylab = "Score", cex.lab = cex_gen,
          cex.axis = cex_gen, ...)
  if (is.null(legend_loc)){
    legend("bottomleft", labels, col = murphy_cols, lwd = 3, bty = "n", 
           cex = cex_gen)
  }
  
  dat
  
}

#' @param HAC_lags Truncation lag for Newey-West variance estimator (needed for confidence bands)
#' @param conf_level Level of confidence bands
#' @rdname murphy
#' @export
murphy_VaRES_diff <- function(TSfc, alpha = 0.025, plot_nr = 2, grid = NULL,
                              labels = NULL, legend_loc = NULL, HAC_lags = 9, 
                              conf_level = 0.95, cex_gen = 1.6, ...){

  # Get grid
  if (is.null(grid)){
    xax <- get_grid(TSfc, plot_nr)  
  } else {
    xax <- grid
  }
  
  # Sample size
  nobs <- length(TSfc[, 5])
  
  # Quantile of limiting normal dist
  scl <- abs(qnorm(0.5 * (1 - conf_level)))
  
  # Enter data here
  df <- matrix(0, length(xax), 5)
  df[, 1] <- xax
    
  # Compute scores, depending on which plot is produced
  if (plot_nr == 1){
    sc_a <- S1(TSfc[, 1], TSfc[, 5], xax, alpha = alpha)
    sc_b <- S1(TSfc[, 3], TSfc[, 5], xax, alpha = alpha)
  } else if (plot_nr == 2){
    sc_a <- S2(TSfc[, 1], TSfc[, 2], TSfc[, 5], xax, alpha = alpha)
    sc_b <- S2(TSfc[, 3], TSfc[, 4], TSfc[, 5], xax, alpha = alpha)
  }
  
  # Mean scores
  df[, 2] <- colMeans(sc_a)
  df[, 3] <- colMeans(sc_b)
    
  # DM test for all threshold values	
  sc_s <- (sc_a - sc_b) %>% (function(z) scale(z, scale = FALSE)) %>%
    (function(z) vSB_mat(z, p = 0.1)/nrow(z)) %>% sqrt %>%
    (function(z) ifelse(!is.na(z), z, 0))
  df[, 4] <- (df[, 2] - df[, 3]) - scl*sc_s
  df[, 5] <- (df[, 2] - df[, 3]) + scl*sc_s

  # Make plot
  matplot(x = df[, 1], y = df[, 4:5], type = "n", ylab = "", 
          bty = "n", cex.axis = cex_gen, cex.lab = cex_gen, 
          col = 1, xaxt = "n", ...)
  aux <- seq(min(df[, 1]), max(df[, 1]), length.out = 5)
  axis(1, at = aux, labels = round(aux), cex.axis = cex_gen)
  polygon(c(df[, 1], rev(df[, 1])), c(df[, 5], rev(df[, 4])), col = "grey", 
          border = NA)
  lines(x = df[, 1], y = (df[, 2] - df[, 3]), type = "l", col = 1, 
        lwd = 2.5)
  abline(h = 0, lty = 2)
  
  # Return data used in plot
  df
  
}

# Function for DM test
t_test <- function(x, eps = 1e-10, k = 2, out = "s"){
  if (sd(x) < eps){
    t <- nw_se <- p <- NA
  } else {
    fit <- lm(x~1)
    nw_se <- sqrt(unname(NeweyWest(fit, lag = k, prewhite = FALSE, adjust = TRUE)))
    t <- unname(coefficients(fit))/nw_se
    p <- pnorm(t)
  }
  if (out == "s"){
    nw_se
  } else if (out == "t"){
    t
  } else if (out == "p"){
    p 
  }
}

# DM stats for VAR-ES forecasts
# Compute DM stats (pointwise for many grid values, both extremal scores)
# ts_b is the bootstrap or permutation sample of the time series
# content of cols: VaR_1, ES_1, VaR_2, ES_2, Realization
# S1_diff are the extremal score differences for the first extremal score as computed by S1, same for S2_diff
# v1, v2 are the grids
# The function returns a one-sided p-value for each grid point of c(v1,v2) for the null hypothesis that the second model (HEAVY) is at most as good as the first one (HS).
DM_stat <- function(ts_b, S1_diff, S2_diff, v1, v2, HAC_lags, alpha = alpha, 
                    score_type = "12"){
  # Extremal scores from bootstrap (or permutation) sample
  S1_Bdiff <- S1(ts_b[,1], ts_b[,5], v1 = v1, alpha = alpha) - 
    S1(ts_b[,3], ts_b[,5], v1 = v1, alpha = alpha)
  S2_Bdiff <- S2(ts_b[,1], ts_b[,2], ts_b[,5], v2 = v2, alpha = alpha) - 
    S2(ts_b[,3], ts_b[,4], ts_b[,5], v2 = v2, alpha = alpha)
  
  # Subtract extremal scores from original sample
  dat1 <- S1_Bdiff - S1_diff
  dat2 <- S2_Bdiff - S2_diff
  
  # Diebold-Mariano tests
  if (score_type %in% c("12", "1")){
    Tstat1 <- t_matC(dat1, k = HAC_lags, eps = 1e-5)
  } else {
    Tstat1 <- c()
  }
  if (score_type %in% c("2", "12")){
    Tstat2 <- t_matC(dat2, k = HAC_lags, eps = 1e-5)
  } else {
    Tstat2 <- c()
  }
  return(c(Tstat1,Tstat2))
}

# Helper function to get subgrid
subgrid <- function(x, k){
  x[seq(from = k, by = k, length.out = floor(length(x)/k))]
}

#' Westfall-young procedure for VaR-ES forecasts
#' @details Tests the hypothesis that the first method weakly dominates the second method. \code{ptest_VaRES} employs a finite grid approximation, whereas \code{ptest_VaRES_full_2_boot} uses analytical calculations. See Section 3.2 in the paper for details.
#' @param TSfc Matrix containing forecasts and realizations. First two columns contain VaR/ES forecasts from method A; cols 3-4 contain VaR/ES forecasts from method B; col 5 contains realizations.
#' @param alpha Scalar, level of VaR/ES to be predicted
#' @param np Nr of iterations (for bootstrap or permutation)
#' @param grid_thin Thinning of grid (1 means no thinning, i.e. all forecasts in sample used as grid points)
#' @param block_length Mean block length for bootstrap or permutation
#' @param HAC_lags Nr of lags for HAC estimator
#' @param score_type Either "12" (both scores) or "2" (second score only)
#' @param grid_bound Upper bound for jumps (only needed if grid_thin > 1)
#' @param sampling_method Either "bootstrap" or "permutation"
#' @param small_output Dummy, set to TRUE to reduce list of outputs (e.g. for use in simulations)
#' @return \code{pval_WY} is the p-value of the full Westfall-Young procedure; \code{pval_OS} is the p-value of the simplified one-step procedure (which is equivalent to Hansen's test)
#' @export
#' @author Alexander Jordan, Fabian Krueger, Johanna Ziegel
ptest_VaRES <- function(TSfc, alpha = 0.025, np = 500, grid_thin = 1,
                        block_length = 0.7343665*(nrow(TSfc)^(1/3)), 
                        HAC_lags = 1, score_type = "2", 
                        grid_bound = NA, sampling_method = "bootstrap",
                        small_output = FALSE){
  
  # Starting time
  t0 <- Sys.time()
  
  # Determine endpoints for tests
  w01 <- sort(unique(unlist(TSfc[, c(1, 3, 5)])))
  w02 <- sort(unique(unlist(TSfc[, c(2, 4)])))
  
  if (identical(grid_thin, "equi")) {
    w02 <- range(w02)
    w02 <- seq(w02[1], w02[2], length.out = dim(TSfc)[1] / 5)
    grid_thin <- 1
  }
  grid_thin <- as.numeric(grid_thin)
  
  # Grid for v1, v2
  v1_grid <- subgrid(w01, grid_thin)
  v2_grid <- subgrid(w02, grid_thin)
  
  # Sample size
  T <- nrow(TSfc)
  
  # Compute score differences
  S1_diff <- S1(TSfc[, 1], TSfc[, 5], v1 = v1_grid, 
                alpha = alpha) - 
    S1(TSfc[, 3], TSfc[, 5], v1 = v1_grid,
       alpha = alpha)
  S2_diff <- S2(TSfc[, 1], TSfc[, 2], TSfc[, 5], v2 = v2_grid,
                alpha = alpha) - 
    S2(TSfc[, 3], TSfc[, 4], TSfc[, 5], v2 = v2_grid, 
       alpha = alpha) 
  
  if (sampling_method == "permutation"){
    rd <- permute_diffs(score_type, S1_diff, S2_diff, block_length, 
                        np, HAC_lags)
  } else if (sampling_method == "bootstrap") {
    rd <- bootstrap_diffs(score_type, S1_diff, S2_diff, block_length,
                          np)
  }
  
  # T-stats for original data
  ty <- rd$ty
  
  # T-stats for Bootstrap data
  Bty <- rd$Bty
  
  # Set NA values to minus infinity (-> no rejection of H0)
  ty[is.na(ty)] <- -Inf
  Bty[is.na(Bty)] <- -Inf

  # Order according to original t-stats
  spi <- sort(ty, index.return = TRUE)$ix  
  # equivalently: spi <- order(ty)
  
  # Nr of evaluation points 
  nc <- dim(Bty)[2]
  
  # Westfall-Young correction
  ry <- rep(NA, nc)
    
  # Two options: Discreteness correction (no/yes)
  if (is.na(grid_bound)){
      
    # Compute corrected p-values (no discreteness correction)
    ry[spi[1]] <- mean(Bty[,spi[1]] >= ty[spi[1]])
    for (k in 2:nc){
      ry[spi[k]] <- (rowMaxs(Bty, cols = spi[1:k]) >= ty[spi[k]]) %>% 
        mean
        
      # Old, less efficient version:
      # TMP <- rowSums(Bty[, spi[1:k]] >= ty[spi[k]])
      # ry[spi[k]] <- mean(TMP >= 1)
        
    }
      
  } else {
      
    # Compute corrected p-values (with discreteness correction)
    for (k in 1:nc){
        
      sel <- which(ty <= (ty[k] + grid_bound))
      ry[k] <- (rowMaxs(Bty, cols = sel) >= (ty[k] - grid_bound)) %>%
        mean
        
      # Old, less efficient version:
      # TMP <- (Bty[, sel] >= (ty[k] - grid_bound))
      # if (length(sel) == 1){
      #   TMP <- as.matrix(TMP)
      # }
      # ry[k] <- mean(rowSums(TMP) >= 1)
        
    }
  
  }
  
  # Time
  t1 <- Sys.time()
  
  ################
  # Returns
  ################
  
  if(small_output){
    list(pval_WY = min(ry), pval_OS = ry[which.max(ty)])
  } else {
    list(run_time = (t1 - t0), pval_WY = ry, pval_OS = ry[which.max(ty)], 
         v1_grid = v1_grid, v2_grid = v2_grid, np = np, 
         block_length = rd$block_length, 
         HAC_lags = HAC_lags, alpha = alpha, simulated_t_stats = Bty, 
         S1_diff = S1_diff, S2_diff = S2_diff)  
  }
}


permute_diffs <- function(score_type, S1_diff, S2_diff, block_length, 
                          np, HAC_lags, eps = 1e-6){
  
  # Sample size
  T <- nrow(S1_diff)
  
  # T-stats for original sample
  if (score_type %in% c("1", "12")){
    ty1 <- t_matC(S1_diff, k = HAC_lags, eps)
  } else {
    ty1 <- c()
  }
  if (score_type %in% c("2", "12")){
    ty2 <- t_matC(S2_diff, k = HAC_lags, eps)
  } else {
    ty2 <- c()
  }
  ty <- c(ty1, ty2)
  
  # Draw signs (permutation)
  T_aux <- ifelse(T %% block_length == 0, T, T - T %% block_length + block_length)
  signs_aux <- matrix(sample(c(-1, 1), size = T_aux*np/block_length, 
                             replace = TRUE), nrow = T_aux/block_length, 
                      ncol = np)
  signs <- apply(signs_aux, 2, function(z) rep(z, each = block_length))[1:T, ]
  
  # Loop
  if (score_type == "12"){
    Bty <- matrix(-1, np, ncol(S1_diff) + ncol(S2_diff))  
  } else {
    Bty <- matrix(-1, np, ncol(S2_diff))  
  }
  
  for (jj in 1:np){  
    if (score_type %in% c("1", "12")){
      S1_tmp <- t_matC(signs[, jj] * S1_diff, k = HAC_lags, eps)
    } else {
      S1_tmp <- c()
    }
    if (score_type %in% c("2", "12")){
      S2_tmp <- t_matC(signs[, jj] * S2_diff, k = HAC_lags, eps)
    } else {
      S2_tmp <- c()
    }
    Bty[jj, ] <- c(S1_tmp, S2_tmp)
  }

  return(list(ty = ty, Bty = Bty))
}

bootstrap_diffs <- function(score_type, S1_diff, S2_diff, block_length, 
                            np, eps = 1e-6){
  
  if (score_type == "2"){
    dat <- S2_diff
  } else if (score_type == "12"){
    dat <- cbind(S1_diff, S2_diff)
  }
  
  #if (block_length == "auto"){
    #block_length <- b.star(dat)[, 1] %>% (function(z) max(c(z, 1)))
  #}
  
  sb <- stationary_boot(dat, np, block_length, eps)
  
  return(list(ty = sb$t_orig, Bty = sb$t_boot, block_length = block_length))
  
}


#' @rdname ptest_VaRES
#' @export
ptest_VaRES_full_2_boot <- function(TSfc, alpha = 0.025, np = 500,
                                    block_length = 10 * (nrow(TSfc) / 2525)^(1/3)) {
  
  # Get forecasts for VaR (v) and Expected Shortfall (e), as well as 
  # realizations (y)
  v <- TSfc[, c(1, 3)]
  e <- TSfc[, c(2, 4)]
  y <- TSfc[, 5]
  
  # Determine endpoints for tests
  pre_grid <- sort(unique(as.vector(unlist(e))))
  
  # Population
  pre_a <- (y < v) * (v - y) / alpha - v
  pre_b <- outer(e, pre_grid, ">=")
  
  a_t <- pre_b[, 1, ] * pre_a[, 1] - pre_b[, 2, ] * pre_a[, 2]
  b_t <- pre_b[, 1, ]              - pre_b[, 2, ]
  
  # Parameters for original data
  ## T statistic numerator sqrt(n) * (a + b * x)
  a_orig <- colMeans2(a_t)
  b_orig <- colMeans2(b_t)
  
  ## T statistic denominator sqrt(c + 2 * d * x + e * x^2) (Politis & Romano)
  n <- length(y)
  i <- outer(1:n, 1:n, function(a, b) abs(a - b))  # lag matrix
  q <- 1 / block_length
  k <- (n - i / n) * (1 - q)^i + i / n * (1 - q)^(n - i)  # kernel weight matrix
  a_t_demean <- a_t - a_orig
  b_t_demean <- b_t - b_orig
  
  cov_param_PR <- function(a, b = NULL) {
    if (is.null(b)) b <- a
    sum(k * (a %o% b)) / n
  }
  
  c_orig <- apply(a_t_demean, 2, cov_param_PR)
  d_orig <-
    sapply(seq_along(pre_grid),
           function(i) cov_param_PR(a_t_demean[, i], b_t_demean[, i]))
  e_orig <- apply(b_t_demean, 2, cov_param_PR)
  
  ##
  params_orig <- data.frame(a = a_orig, b = b_orig, c = c_orig, d = d_orig, e = e_orig)
  
  # Parameters 
  
  # Parameters for bootstrap samples (T statistic numerator)
  ind <- tsboot(1:n, identity, R = np, l = block_length, sim = "geom")$t
  params_boot <- apply(ind, 1, function(i) {
    data.frame(
      a = colMeans2(a_t, rows = i) - a_orig,
      b = colMeans2(b_t, rows = i) - b_orig
    )})
  
  # Find additional grid points
  expr <- expression((a * d - b * c) / (b * d - a * e))
  candidates <-
    with(params_orig,
         cbind(
           eval(expr),
           sapply(seq_along(params_boot),
                  function(i) with(params_boot[[i]], eval(expr)))
         ))
  bins <- seq_along(pre_grid)
  which_bin_cand <- apply(candidates, 2, findInterval, pre_grid, left.open = TRUE) + 1
  accepted_cand <- candidates[which(bins == which_bin_cand)]
  
  # Compute T statistics
  grid <- list(
    upperlimits = pre_grid,
    lowerlimits = head(pre_grid, -1),
    inbetweeners = accepted_cand)
  
  helpTnumerator <- function(grid_type, params) {
    x <- grid[[grid_type]]
    expr <- expression(sqrt(n) * (a + b * x))
    if (grid_type == 1L) {
      with(params, eval(expr))
    } else if (grid_type == 2L) {
      with(tail(params, -1), eval(expr))
    } else if (grid_type == 3L) {
      with(params[findInterval(x, pre_grid) + 1, ], eval(expr))
    }
  }
  calcTnum <- function(params) {
    do.call(c, lapply(seq_along(grid), helpTnumerator, params))
  }
  helpTdenominator <- function(grid_type, params) {
    x <- grid[[grid_type]]
    expr <- expression(sqrt(c + 2 * d * x + e * x^2))
    if (grid_type == 1L) {
      with(params, eval(expr))
    } else if (grid_type == 2L) {
      with(tail(params, -1), eval(expr))
    } else if (grid_type == 3L) {
      with(params[findInterval(x, pre_grid) + 1, ], eval(expr))
    }
  }
  Tdenominator <- do.call(c, lapply(seq_along(grid), helpTdenominator, params_orig))
  
  ## original data
  ty <- calcTnum(params_orig) / Tdenominator
  
  ## bootstrap data
  Bty <- do.call(cbind, lapply(params_boot, calcTnum)) / Tdenominator
  
  # Set NA values to minus infinity (-> no rejection of H0)
  ty[is.na(ty)] <- -Inf
  Bty[is.na(Bty)] <- -Inf
  
  # Order according to original t-stats
  spi <- sort(ty, index.return = TRUE)$ix  
  # equivalently: spi <- order(ty)
  
  # Nr of evaluation points 
  nc <- nrow(Bty)
  
  # Westfall-Young correction
  ry <- rep(NA, nc)
  
  # Compute corrected p-values (no discreteness correction)
  currentBMax <- Bty[spi[1], ]
  ry[spi[1]] <- mean(currentBMax >= ty[spi[1]])
  for (k in 2:nc) {
    currentBMax <- pmax(currentBMax, Bty[spi[k], ])
    ry[spi[k]] <- mean(currentBMax >= ty[spi[k]])
  }
  
  ################
  # Returns
  ################
  
  list(pval_WY = min(ry), pval_OS = ry[which.max(ty)])
}


#' Simulate data from the design of Section 4 in the paper.
#' @param T Length of series to be simulated. 
#' @param s2_0,a,b,t_df Parameters of time series model which generates 'true' forecasts
#' @param zeta1,zeta2 Variances of disturbances
#' @param alpha Level of VaR/ES forecasts
#' @details Simulate time series of VaR/ES forecasts and corresponding realizations, as described in the paper by Ziegel et al (2018). 
#' @export
sim_VaRES <- function(T, s2_0 = 0.35, a = 0.5, b = 0.7, 
                      zeta1 = 1, zeta2 = 1, t_df = 6, alpha = 0.025){
  
  # Simulate 'realized measure' (coefficients from AR1 fitted to log of SP500 RK)
  rm <- (-0.6223316 + arima.sim(list(ar = 0.8316009), n = T,
                                rand.gen = function(n) rnorm(n, sd = 0.6185419),
                                n.start = 500)) %>%
    exp

  # Simulate variance sequence
  s2 <- rep(0, T)
  s2[1] <- s2_0
  for (jj in 2:T){
    s2[jj] <- a*rm[jj-1] + b*s2[jj-1]
  }
  # Scale factor to make t variance equal to one
  scale_fac <- sqrt((t_df-2)/t_df)
  # VaR and ES of t distributed variable
  t_quant <- varT(p = alpha, n = t_df)
  t_es <- esT(p = alpha, n = t_df)

  # Simulate returns
  y <- sqrt(s2) * scale_fac * rt(T, df = t_df)

  # True VaR and ES values
  var_true <- sqrt(s2) * scale_fac * t_quant
  es_true <- sqrt(s2) * scale_fac * t_es

  # Errors (epsilon_t) of the two forecasters
  e1 <- rnorm(T, sd = sqrt(zeta1))
  e2 <- rnorm(T, sd = sqrt(zeta2))

  # Forecasts
  var1 <- var_true + e1
  es1 <- es_true + e1
  var2 <- var_true + e2
  es2 <- es_true + e2
  TSfc <- cbind(var1, es1, var2, es2, y)

  # Return
  return(TSfc)

}

# Helper function to get forecast data into 
# format required by some functions
make_TSfc <- function(m1, m2, dat){
  sel <- c(paste0(c("var_", "es_"), m1),
           paste0(c("var_", "es_"), m2),
           "rlz")
  as.matrix(dat[, sel])
}
