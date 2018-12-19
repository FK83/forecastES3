# Objective function for fitting HEAVY model
ff <- function(theta, rm, ret, df_fix = 7){

  # Construct parameters
  pars <- heavy_helper(theta, rm, ret, df_fix)
  
  # Log likelihood
  ll <- dts_log(ret[-1], pars$s_2t, pars$df)
  
  # Return negative log-likelihood
  return(-sum(ll))
}

# Function to compute VaR and ES of HEAVY model
heavy_fc <- function(theta, rm, ret, df_fix = 7, alpha = 0.025){
  
  # Construct parameters
  pars <- heavy_helper(theta, rm, ret, df_fix)
  
  # Compute and return (vola, VaR, ES)
  c(pars$s_tp1, vets(pars$s_tp1, pars$df, alpha))
  
}

vets <- function(s, df, alpha){
  
  # Get (VaR, ES) of standard t dist
  tv <- varT(p = alpha, n = df)
  te <- esT(p = alpha, n = df)
  
  # Return (VaR, ES)
  scale_fac <- sqrt((df-2)/df)
  s*scale_fac*c(tv, te)
  
}

# Helper function to construct HEAVY stuff
heavy_helper <- function(theta, rm, ret, df_fix){
  
  # Length of time series
  t <- length(rm)
  
  # Truncation lag
  trunc <- floor(sqrt(t))
  
  # Initialization
  init <- (1/sqrt(t))*sum(ret[1:trunc]^2)  
  
  # Parameters
  w <- exp(theta[1])
  a <- theta[2]
  b <- theta[3]
  
  if (length(theta) == 4){
    # df for t-distribution 
    df <- theta[4]
  } else {
    df <- df_fix
  }
  
  # Construct variances
  all <- rep(0, t+1)
  all[1] <- init
  for (jj in 2:(t+1)){
    all[jj] <- w + a*rm[jj-1] + b*all[jj-1]
  }
  
  # Standard deviation for periods 2 to t (relevant for likelihood)
  s_2t <- sqrt(all[-c(1, t+1)])
  
  # Standard deviation for period t+1 (relevant for forecast)
  s_tp1 <- sqrt(all[t+1])
  
  # Return outputs
  list(s_2t = s_2t, s_tp1 = s_tp1, df = df)
  
}

# Density of scaled t-distribution 
# (with sd s and df degrees of freedom)
#' @importFrom stats dt integrate quantile
dts <- function(x, s, df){
  a <- sqrt((df-2)/df)
  as <- a*s
  dt(x/as, df = df)/as
}

# Log density of scaled t-dist
dts_log <- function(x, s, df){
  a <- sqrt((df-2)/df)
  as <- a*s
  (dt(x/as, df = df, log = TRUE) - log(as))
}

# Helper function
m <- function(x, order = 1, s = 1, df = 7){
  dts(x, s = s, df = df)*(x^order)
}

# Compute variance of scaled t-dist
# (Numerically, to check parametrization)
check_s <- function(s, df){
  x2 <- integrate(m, -Inf, Inf, order = 2, s = s, df = df)$value
  x1 <- integrate(m, -Inf, Inf, order = 1, s = s, df = df)$value
  sqrt(x2 - x1^2)
}

# Empirical ES
es_emp <- function(data, p){
  q <- quantile(data, p)
  mean(data[data < q])
}