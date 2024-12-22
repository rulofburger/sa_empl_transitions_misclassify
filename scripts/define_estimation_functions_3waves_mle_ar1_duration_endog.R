# DEFINE FUNCTIONS ====

#> General functions that restrict parameter values to unit interval

minval <- 1e-100

logit_transform <- function(param0) {
  log(param0/(1 - param0))
}

logit_inverse <- function(param_input0) {
  1/(1 + exp(-param_input0))
}

# partial_exGauss_lambda <- function(x, sigma, lambda) {
#   # Compute intermediate values
#   u <- x / sigma - sigma / lambda  # u = x/sigma - sigma/lambda
#   phi_u <- dnorm(u)                # PDF of the standard normal at u
#   Phi_u <- pnorm(u)                # CDF of the standard normal at u
#   
#   # Compute the exponential term
#   exp_term <- exp(-x / lambda + (sigma^2) / (2 * lambda^2))
#   
#   # Term 1: Derivative of the exponential part
#   dh_dlambda <- x / (lambda^2) - sigma^2 / (lambda^3)
#   term1 <- exp_term * (-1 / (lambda^2) + 1 / lambda * dh_dlambda)
#   
#   # Term 2: Derivative of the CDF part
#   du_dlambda <- sigma / (lambda^2)  # Derivative of u with respect to lambda
#   term2 <- (1 / lambda) * exp_term * phi_u * du_dlambda
#   
#   # Total derivative
#   term1 + term2
# }

# Function to compute the partial derivative of the exGaussian PDF with respect to lambda
partial_exGauss_lambda <- function(x, sigma, lambda) {
  # Exponential term
  exp_term <- exp(-x / lambda + sigma^2 / (2 * lambda^2))
  
  # CDF term
  cdf_term <- pnorm(x / sigma - sigma / lambda)
  
  # PDF term
  pdf_term <- dnorm(x / sigma - sigma / lambda)
  
  # Term 1: Contribution from the exponential term
  term1 <- exp_term * cdf_term * (-1 / lambda^2 + x / lambda^3 - sigma^2 / lambda^4)
  
  # Term 2: Contribution from the CDF term
  term2 <- (1 / lambda) * exp_term * pdf_term * (sigma / lambda^2)
  
  # Total derivative
  return(term1 + term2)
}

partial_exGauss_sigma <- function(x, sigma, lambda) {
  # Compute intermediate values
  u <- x / sigma - sigma / lambda  # u = x/sigma - sigma/lambda
  phi_u <- dnorm(u)                # PDF of the standard normal at u
  Phi_u <- pnorm(u)                # CDF of the standard normal at u
  
  # Compute Term 1: Derivative of the exponential part
  term1 <- (1 / lambda) * exp(-x / lambda + (sigma^2) / (2 * lambda^2)) * (sigma / lambda^2) * Phi_u
  
  # Compute Term 2: Derivative of the CDF part
  du_dsigma <- -x / (sigma^2) - 1 / lambda  # Derivative of u with respect to sigma
  term2 <- (1 / lambda) * exp(-x / lambda + (sigma^2) / (2 * lambda^2)) * phi_u * du_dsigma
  
  # Total derivative
  term1 + term2
  
}

partial_phi_sigma <- function(x, sigma) {
  dnorm(x, mean = 0, sd = sigma)*((x^2)/(sigma^3) - 1/sigma)
}

partial_diff_phi_sigma <- function(x, sigma) {
  dnorm(x, mean = 0, sd = sqrt(2*(sigma)^2))*((x^2)/(2*(sigma^3)) - 1/sigma)
}


# Load the stats package for normal distribution functions
library(stats)





# DURATION DEPENDENCE (timegap + tenure) ====


calc_lli_3waves_ar1_duration_endog <- function(param_transformed, pi_fixed = NULL, only_pi = NULL) {
  
  param <- param_transformed
  if(!is.null(pi_fixed)) {
    param$pi <- 0.5*logit_inverse(pi_fixed$pi)
  } else {
    param$pi <- 0.5*logit_inverse(param_transformed$pi)
  }
  if(!is.null(only_pi)) {
    param$theta_0 <- logit_inverse(only_pi$theta_0)
    param$theta_1 <- logit_inverse(only_pi$theta_1)
    param$lambda_g <- exp(only_pi$lambda_g) + 0.1
    param$lambda_h <- exp(only_pi$lambda_h) + 0.1
    param$sigma_g <- exp(only_pi$sigma_g) + 0.1
    param$sigma_h <- exp(only_pi$sigma_h) + 0.1
  } else {
    param$theta_0 <- logit_inverse(param_transformed$theta_0)
    param$theta_1 <- logit_inverse(param_transformed$theta_1)
    param$lambda_g <- exp(param_transformed$lambda_g) + 0.1
    param$lambda_h <- exp(param_transformed$lambda_h) + 0.1
    param$sigma_g <- exp(param_transformed$sigma_g) + 0.1
    param$sigma_h <- exp(param_transformed$sigma_h) + 0.1
  }
  
  mu <- param$theta_0/(param$theta_1 + param$theta_0)

  df_probs_temp <- df_template_duration_endog %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu, 1 - mu),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
        y1_star == 0 & y2_star == 1 ~ param$theta_0,
        y1_star == 1 & y2_star == 0 ~ param$theta_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta_0,
        y2_star == 0 & y3_star == 1 ~ param$theta_0,
        y2_star == 1 & y3_star == 0 ~ param$theta_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta_1
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      p_g1 = case_when(
        y1 == 0 ~ 1,
        y1 == 1 ~ gamlss.dist::dexGAUS(tenure1, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)
      ),
      p_g2 = case_when(
        y2 == 0 ~ 1,
        y2 == 1 & y2_star == 0 ~ gamlss.dist::dexGAUS(tenure2, mu = 0, sigma = param$sigma_g, nu = param$lambda_g),
        y2 == 1 & y2_star == 1 & y1_star == 0 ~ dnorm(tenure2 - 0.125, mean = 0, sd = param$sigma_g),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 1 ~ dnorm(tenure2 - tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure2 - 0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)
      ),
      p_g3 = case_when(
        y3 == 0 ~ 1,
        y3 == 1 & y3_star == 0 ~ gamlss.dist::dexGAUS(tenure3, mu = 0, sigma = param$sigma_g, nu = param$lambda_g),
        y3 == 1 & y3_star == 1 & y2_star == 0 ~ dnorm(tenure3 - 0.125, mean = 0, sd = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 1 ~ dnorm(tenure3 - tenure2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 0 ~ dnorm(tenure3 - 0.375, mean = 0, sd = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 1 ~ dnorm(tenure3 - tenure1 - 0.5, mean = 0, sd = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure3 - 0.5, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)
      ),
      p_h1 = case_when(
        y1 == 1 ~ 1,
        y1 == 0 ~ gamlss.dist::dexGAUS(timegap1, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)
      ),
      p_h2 = case_when(
        y2 == 1 ~ 1,
        y2 == 0 & y2_star == 1 ~ gamlss.dist::dexGAUS(timegap2, mu = 0, sigma = param$sigma_h, nu = param$lambda_h),
        y2 == 0 & y2_star == 0 & y1_star == 1 ~ dnorm(timegap2 - 0.125, mean = 0, sd = param$sigma_h),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 0 ~ dnorm(timegap2 - timegap1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap2 - 0.25, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)
      ),
      p_h3 = case_when(
        y3 == 1 ~ 1,
        y3 == 0 & y3_star == 1 ~ gamlss.dist::dexGAUS(timegap3, mu = 0, sigma = param$sigma_h, nu = param$lambda_h),
        y3 == 0 & y3_star == 0 & y2_star == 1 ~ dnorm(timegap3 - 0.125, mean = 0, sd = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 0 ~ dnorm(timegap3 - timegap2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 1 ~ dnorm(timegap3 - 0.375, mean = 0, sd = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 0 ~ dnorm(timegap3 - timegap1 - 0.5, mean = 0, sd = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap3 - 0.5, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)
      ),
      p1_star = if_else(p1_star < minval, minval, p1_star),
      p2_star = if_else(p2_star < minval, minval, p2_star),
      p3_star = if_else(p3_star < minval, minval, p3_star),
      p1 = if_else(p1 < minval, minval, p1),
      p2 = if_else(p2 < minval, minval, p2),
      p3 = if_else(p3 < minval, minval, p3),
      p_g1 = if_else(p_g1 < minval, minval, p_g1),
      p_g2 = if_else(p_g2 < minval, minval, p_g2),
      p_g3 = if_else(p_g3 < minval, minval, p_g3),
      p_h1 = if_else(p_h1 < minval, minval, p_h1),
      p_h2 = if_else(p_h2 < minval, minval, p_h2),
      p_h3 = if_else(p_h3 < minval, minval, p_h3),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3*p_g1*p_g2*p_g3*p_h1*p_h2*p_h3
    ) 
  
  df_probs <- df_probs_temp %>% 
    group_by(y1, y2, y3, timegap1, timegap2, timegap3, tenure1, tenure2, tenure3) %>% 
    summarise(joint_p = sum(joint_p), .groups = "drop") 
  
  df_lli <- df_estimate %>% 
    left_join(df_probs, by = c('y1', 'y2', 'y3', 'timegap1', 'timegap2', 'timegap3', 'tenure1', 'tenure2', 'tenure3')) %>% 
    mutate(lli = weight*log(joint_p)) %>%
    pull(lli)

}

calc_lli_derivatives_3waves_ar1_duration_endog <- function(param_transformed, pi_fixed = NULL, only_pi = NULL) {
 
  param <- param_transformed
  if(!is.null(pi_fixed)) {
    param$pi <- 0.5*logit_inverse(pi_fixed$pi)
  } else {
    param$pi <- 0.5*logit_inverse(param_transformed$pi)
  }
  if(!is.null(only_pi)) {
    param$theta_0 <- logit_inverse(only_pi$theta_0)
    param$theta_1 <- logit_inverse(only_pi$theta_1)
    param$lambda_g <- exp(only_pi$lambda_g) + minval
    param$lambda_h <- exp(only_pi$lambda_h) + minval
    param$sigma_g <- exp(only_pi$sigma_g) + minval
    param$sigma_h <- exp(only_pi$sigma_h) + minval
  } else {
    param$theta_0 <- logit_inverse(param_transformed$theta_0)
    param$theta_1 <- logit_inverse(param_transformed$theta_1)
    param$lambda_g <- exp(param_transformed$lambda_g) + minval
    param$lambda_h <- exp(param_transformed$lambda_h) + minval
    param$sigma_g <- exp(param_transformed$sigma_g) + minval
    param$sigma_h <- exp(param_transformed$sigma_h) + minval
  }
  
  mu <- param$theta_0/(param$theta_1 + param$theta_0)
  
  df_probs_temp <- df_template_duration_endog %>% 
    mutate(
      p1_star = if_else(y1_star == 1, mu, 1 - mu),
      p2_star = case_when(
        y1_star == 0 & y2_star == 0 ~ 1 - param$theta_0,
        y1_star == 0 & y2_star == 1 ~ param$theta_0,
        y1_star == 1 & y2_star == 0 ~ param$theta_1,
        y1_star == 1 & y2_star == 1 ~ 1 - param$theta_1
      ),
      p3_star = case_when(
        y2_star == 0 & y3_star == 0 ~ 1 - param$theta_0,
        y2_star == 0 & y3_star == 1 ~ param$theta_0,
        y2_star == 1 & y3_star == 0 ~ param$theta_1,
        y2_star == 1 & y3_star == 1 ~ 1 - param$theta_1
      ),
      p1 = if_else(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = if_else(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = if_else(y3 == y3_star, 1 - param$pi, param$pi),
      p_g1 = case_when(
        y1 == 0 ~ 1,
        y1 == 1 ~ gamlss.dist::dexGAUS(tenure1, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)
      ),
      p_g2 = case_when(
        y2 == 0 ~ 1,
        y2 == 1 & y2_star == 0 ~ gamlss.dist::dexGAUS(tenure2, mu = 0, sigma = param$sigma_g, nu = param$lambda_g),
        y2 == 1 & y2_star == 1 & y1_star == 0 ~ dnorm(tenure2 - 0.125, mean = 0, sd = param$sigma_g),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 1 ~ dnorm(tenure2 - tenure1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure2 - 0.25, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)
      ),
      p_g3 = case_when(
        y3 == 0 ~ 1,
        y3 == 1 & y3_star == 0 ~ gamlss.dist::dexGAUS(tenure3, mu = 0, sigma = param$sigma_g, nu = param$lambda_g),
        y3 == 1 & y3_star == 1 & y2_star == 0 ~ dnorm(tenure3 - 0.125, mean = 0, sd = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 1 ~ dnorm(tenure3 - tenure2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_g)^2)),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 0 ~ dnorm(tenure3 - 0.375, mean = 0, sd = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 1 ~ dnorm(tenure3 - tenure1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_g)^2)),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 0 ~ gamlss.dist::dexGAUS(tenure3 - 0.5, mu = 0, sigma = param$sigma_g, nu = param$lambda_g)
      ),
      p_h1 = case_when(
        y1 == 1 ~ 1,
        y1 == 0 ~ gamlss.dist::dexGAUS(timegap1, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)
      ),
      p_h2 = case_when(
        y2 == 1 ~ 1,
        y2 == 0 & y2_star == 1 ~ gamlss.dist::dexGAUS(timegap2, mu = 0, sigma = param$sigma_h, nu = param$lambda_h),
        y2 == 0 & y2_star == 0 & y1_star == 1 ~ dnorm(timegap2 - 0.125, mean = 0, sd = param$sigma_h),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 0 ~ dnorm(timegap2 - timegap1 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap2 - 0.25, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)
      ),
      p_h3 = case_when(
        y3 == 1 ~ 1,
        y3 == 0 & y3_star == 1 ~ gamlss.dist::dexGAUS(timegap3, mu = 0, sigma = param$sigma_h, nu = param$lambda_h),
        y3 == 0 & y3_star == 0 & y2_star == 1 ~ dnorm(timegap3 - 0.125, mean = 0, sd = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 0 ~ dnorm(timegap3 - timegap2 - 0.25, mean = 0, sd = sqrt(2*(param$sigma_h)^2)),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 1 ~ dnorm(timegap3 - 0.375, mean = 0, sd = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 0 ~ dnorm(timegap3 - timegap1 - 0.5, mean = 0, sd = sqrt(2*(param$sigma_h)^2)),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 1 ~ gamlss.dist::dexGAUS(timegap3 - 0.5, mu = 0, sigma = param$sigma_h, nu = param$lambda_h)
      ),
      p1_star = if_else(p1_star < minval, minval, p1_star),
      p2_star = if_else(p2_star < minval, minval, p2_star),
      p3_star = if_else(p3_star < minval, minval, p3_star),
      p1 = if_else(p1 < minval, minval, p1),
      p2 = if_else(p2 < minval, minval, p2),
      p3 = if_else(p3 < minval, minval, p3),
      p_g1 = if_else(p_g1 < minval, minval, p_g1),
      p_g2 = if_else(p_g2 < minval, minval, p_g2),
      p_g3 = if_else(p_g3 < minval, minval, p_g3),
      p_h1 = if_else(p_h1 < minval, minval, p_h1),
      p_h2 = if_else(p_h2 < minval, minval, p_h2),
      p_h3 = if_else(p_h3 < minval, minval, p_h3),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3*p_g1*p_g2*p_g3*p_h1*p_h2*p_h3,
      d1_star_theta_0 = case_when(
        y1_star == 0 ~ -param$theta_1/((param$theta_1 + param$theta_0)^2),
        y1_star == 1 ~ param$theta_1/((param$theta_1 + param$theta_0)^2)
      ),
      d1_star_theta_1 = case_when(
        y1_star == 0 ~ param$theta_0/((param$theta_1 + param$theta_0)^2),
        y1_star == 1 ~ -param$theta_0/((param$theta_1 + param$theta_0)^2)
      ),
      d2_star_theta_0 = case_when(
        y1_star == 0 & y2_star == 0 ~ -1,
        y1_star == 0 & y2_star == 1 ~ 1,
        y1_star == 1 & y2_star == 0 ~ 0,
        y1_star == 1 & y2_star == 1 ~ 0
      ),
      d2_star_theta_1 = case_when(
        y1_star == 0 & y2_star == 0 ~ 0,
        y1_star == 0 & y2_star == 1 ~ 0,
        y1_star == 1 & y2_star == 0 ~ 1,
        y1_star == 1 & y2_star == 1 ~ -1
      ),
      d3_star_theta_0 = case_when(
        y2_star == 0 & y3_star == 0 ~ -1,
        y2_star == 0 & y3_star == 1 ~ 1,
        y2_star == 1 & y3_star == 0 ~ 0,
        y2_star == 1 & y3_star == 1 ~ 0
      ),
      d3_star_theta_1 = case_when(
        y2_star == 0 & y3_star == 0 ~ 0,
        y2_star == 0 & y3_star == 1 ~ 0,
        y2_star == 1 & y3_star == 0 ~ 1,
        y2_star == 1 & y3_star == 1 ~ -1
      ),
      d1_pi = if_else(y1 == y1_star, -1, 1),
      d2_pi = if_else(y2 == y2_star, -1, 1),
      d3_pi = if_else(y3 == y3_star, -1, 1),
      d1_lambda_g = case_when(
        y1 == 0 ~ 0,
        y1 == 1 ~ partial_exGauss_lambda(x = tenure1, sigma = param$sigma_g, lambda = param$lambda_g)
      ),
      d1_sigma_g = case_when(
        y1 == 0 ~ 0,
        y1 == 1 ~ partial_exGauss_sigma(x = tenure1, sigma = param$sigma_g, lambda = param$lambda_g)
      ),
      d1_lambda_h = case_when(
        y1 == 1 ~ 0,
        y1 == 0 ~ partial_exGauss_lambda(x = timegap1, sigma = param$sigma_h, lambda = param$lambda_h)
      ),
      d1_sigma_h = case_when(
        y1 == 1 ~ 0,
        y1 == 0 ~ partial_exGauss_sigma(x = timegap1, sigma = param$sigma_h, lambda = param$lambda_h)
      ),
      d2_lambda_g = case_when(
        y2 == 0 ~ 0,
        y2 == 1 & y2_star == 0 ~ partial_exGauss_lambda(x = tenure2, sigma = param$sigma_g, lambda = param$lambda_g),
        y2 == 1 & y2_star == 1 & y1_star == 0 ~ 0,
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 1 ~ 0,
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 0 ~ partial_exGauss_lambda(x = tenure2 - 0.25, sigma = param$sigma_g, lambda = param$lambda_g)
      ),
      d2_sigma_g = case_when(
        y2 == 0 ~ 0,
        y2 == 1 & y2_star == 0 ~ partial_exGauss_sigma(x = tenure2, sigma = param$sigma_g, lambda = param$lambda_g),
        y2 == 1 & y2_star == 1 & y1_star == 0 ~ partial_phi_sigma(x = tenure2 - 0.125, sigma = param$sigma_g),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 1 ~ partial_diff_phi_sigma(x = tenure2 - tenure1 - 0.25, sigma = param$sigma_g),
        y2 == 1 & y2_star == 1 & y1_star == 1 & y1 == 0 ~ partial_exGauss_sigma(x = tenure2 - 0.25, sigma = param$sigma_g, lambda = param$lambda_g)
      ),
      d2_lambda_h = case_when(
        y2 == 1 ~ 0,
        y2 == 0 & y2_star == 1 ~ partial_exGauss_lambda(x = timegap2, sigma = param$sigma_h, lambda = param$lambda_h),
        y2 == 0 & y2_star == 0 & y1_star == 1 ~ 0,
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 0 ~ 0,
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 1 ~ partial_exGauss_lambda(x = timegap2 - 0.25, sigma = param$sigma_h, lambda = param$lambda_h)
      ),
      d2_sigma_h = case_when(
        y2 == 1 ~ 0,
        y2 == 0 & y2_star == 1 ~ partial_exGauss_sigma(x = timegap2, sigma = param$sigma_h, lambda = param$lambda_h),
        y2 == 0 & y2_star == 0 & y1_star == 1 ~ partial_phi_sigma(x = timegap2 - 0.125, sigma = param$sigma_h),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 0 ~ partial_diff_phi_sigma(x = timegap2 - timegap1 - 0.25, sigma = param$sigma_h),
        y2 == 0 & y2_star == 0 & y1_star == 0 & y1 == 1 ~ partial_exGauss_sigma(x = timegap2 - 0.25, sigma = param$sigma_h, lambda = param$lambda_h)
      ),
      d3_lambda_g = case_when(
        y3 == 0 ~ 0,
        y3 == 1 & y3_star == 0 ~ partial_exGauss_lambda(x = tenure3, sigma = param$sigma_g, lambda = param$lambda_g),
        y3 == 1 & y3_star == 1 & y2_star == 0 ~ 0,
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 1 ~ 0,
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 0 ~ 0,
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 1 ~ 0,
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 0 ~ partial_exGauss_lambda(x = tenure3 - 0.5, sigma = param$sigma_g, lambda = param$lambda_g)
      ),
      d3_sigma_g = case_when(
        y3 == 0 ~ 0,
        y3 == 1 & y3_star == 0 ~ partial_exGauss_sigma(x = tenure3, sigma = param$sigma_g, lambda = param$lambda_g),
        y3 == 1 & y3_star == 1 & y2_star == 0 ~ partial_phi_sigma(x = tenure3 - 0.125, sigma = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 1 ~ partial_diff_phi_sigma(x = tenure3 - tenure2 - 0.25, sigma = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 0 ~ partial_phi_sigma(x = tenure3 - 0.375, sigma = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 1 ~ partial_diff_phi_sigma(x = tenure3 - tenure1 - 0.5, sigma = param$sigma_g),
        y3 == 1 & y3_star == 1 & y2_star == 1 & y2 == 0 & y1_star == 1 & y1 == 0 ~ partial_exGauss_sigma(x = tenure3 - 0.5, sigma = param$sigma_g, lambda = param$lambda_g)
      ),
      d3_lambda_h = case_when(
        y3 == 1 ~ 0,
        y3 == 0 & y3_star == 1 ~ partial_exGauss_lambda(x = timegap3, sigma = param$sigma_h, lambda = param$lambda_h),
        y3 == 0 & y3_star == 0 & y2_star == 1 ~ 0,
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 0 ~ 0,
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 1 ~ 0,
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 0 ~ 0,
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 1 ~ partial_exGauss_lambda(x = timegap3 - 0.5, sigma = param$sigma_h, lambda = param$lambda_h)
      ),
      d3_sigma_h = case_when(
        y3 == 1 ~ 0,
        y3 == 0 & y3_star == 1 ~ partial_exGauss_sigma(x = timegap3, sigma = param$sigma_h, lambda = param$lambda_h),
        y3 == 0 & y3_star == 0 & y2_star == 1 ~ partial_phi_sigma(x = timegap3 - 0.125, sigma = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 0 ~ partial_diff_phi_sigma(x = timegap3 - timegap2 - 0.25, sigma = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 1 ~ partial_phi_sigma(x = timegap3 - 0.375, sigma = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 0 ~ partial_diff_phi_sigma(x = timegap3 - timegap1 - 0.5, sigma = param$sigma_h),
        y3 == 0 & y3_star == 0 & y2_star == 0 & y2 == 1 & y1_star == 0 & y1 == 1 ~ partial_exGauss_sigma(x = timegap3 - 0.5, sigma = param$sigma_h, lambda = param$lambda_h)
      ),      
      joint_d_theta_0 = 
        d1_star_theta_0*joint_p/p1_star + 
        d2_star_theta_0*joint_p/p2_star + 
        d3_star_theta_0*joint_p/p3_star,
      joint_d_theta_1 = 
        d1_star_theta_1*joint_p/p1_star + 
        d2_star_theta_1*joint_p/p2_star + 
        d3_star_theta_1*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3,
      joint_d_lambda_g = 
        d1_lambda_g*joint_p/p_g1 + 
        d2_lambda_g*joint_p/p_g2 + 
        d3_lambda_g*joint_p/p_g3,
      joint_d_sigma_g = 
        d1_sigma_g*joint_p/p_g1 + 
        d2_sigma_g*joint_p/p_g2 + 
        d3_sigma_g*joint_p/p_g3,
      joint_d_lambda_h = 
        d1_lambda_h*joint_p/p_h1 + 
        d2_lambda_h*joint_p/p_h2 + 
        d3_lambda_h*joint_p/p_h3,
      joint_d_sigma_h = 
        d1_sigma_h*joint_p/p_h1 + 
        d2_sigma_h*joint_p/p_h2 + 
        d3_sigma_h*joint_p/p_h3
    ) 
  
  df_grad <- df_probs_temp %>% 
    group_by(y1, y2, y3, timegap1, timegap2, timegap3, tenure1, tenure2, tenure3) %>% 
    summarise(
      joint_d_theta_0 = sum(joint_d_theta_0), 
      joint_d_theta_1 = sum(joint_d_theta_1),
      joint_d_pi = sum(joint_d_pi),
      joint_d_lambda_g = sum(joint_d_lambda_g),
      joint_d_sigma_g = sum(joint_d_sigma_g),
      joint_d_lambda_h = sum(joint_d_lambda_h),
      joint_d_sigma_h = sum(joint_d_sigma_h),
      joint_p = sum(joint_p),
      .groups = "drop") %>% 
    mutate(
      joint_d_theta_0 = joint_d_theta_0*param$theta_0*(1 - param$theta_0), 
      joint_d_theta_1 = joint_d_theta_1*param$theta_1*(1 - param$theta_1),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)*0.5,
      joint_d_lambda_g = joint_d_lambda_g*param$lambda_g,
      joint_d_sigma_g = joint_d_sigma_g*param$sigma_g,
      joint_d_lambda_h = joint_d_lambda_h*param$lambda_h,
      joint_d_sigma_h = joint_d_sigma_h*param$sigma_h
    )
  
  df_gi <- df_estimate %>% 
    # left_join(df_grad, by = c('y1', 'y2', 'y3', 'timegap1', 'timegap2', 'tenure1', 'tenure2')) %>% 
    left_join(df_grad, by = c('y1', 'y2', 'y3', 'timegap1', 'timegap2', 'timegap3', 'tenure1', 'tenure2', 'tenure3')) %>% 
    mutate(
      lgi_theta_0 = weight*joint_d_theta_0/joint_p,
      lgi_theta_1 = weight*joint_d_theta_1/joint_p,
      lgi_pi = weight*joint_d_pi/joint_p,
      lgi_lambda_g = weight*joint_d_lambda_g/joint_p,
      lgi_sigma_g = weight*joint_d_sigma_g/joint_p,
      lgi_lambda_h = weight*joint_d_lambda_h/joint_p,
      lgi_sigma_h = weight*joint_d_sigma_h/joint_p,
    ) %>%
    select(lgi_theta_0, lgi_theta_1, lgi_lambda_g, lgi_lambda_h, lgi_sigma_g, lgi_sigma_h, lgi_pi)
  
  if(!is.null(pi_fixed)) df_gi <- df_gi %>% select(lgi_theta_0, lgi_theta_1, lgi_lambda_g, lgi_lambda_h, lgi_sigma_g, lgi_sigma_h)
  if(!is.null(only_pi)) df_gi <- df_gi %>% select(lgi_pi)
  return(df_gi)

}

calc_mle_3waves_ar1_duration_endog <- function(param_transformed) {
  print(param_transformed)  
  ll <- sum(calc_lli_3waves_ar1_duration_endog(param_transformed))
  print(ll)
  return(ll)
}

calc_mle_derivatives_3waves_ar1_duration_endog <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_duration_endog(param_transformed))
}

calc_mle_3waves_ar1_duration_endog_pi_fixed <- function(param_transformed, pi_fixed) {
  print(param_transformed)
  print(pi_fixed)
  ll <- sum(calc_lli_3waves_ar1_duration_endog(param_transformed = param_transformed, pi_fixed = pi_fixed))
  print(ll)
  return(ll)
}

calc_mle_derivatives_3waves_ar1_duration_endog_pi_fixed <- function(param_transformed, pi_fixed) {
  colSums(calc_lli_derivatives_3waves_ar1_duration_endog(param_transformed, pi_fixed = pi_fixed))
}

calc_mle_3waves_ar1_duration_endog_only_pi <- function(param_transformed, only_pi) {
  print(param_transformed)
  print(only_pi)
  ll <- sum(calc_lli_3waves_ar1_duration_endog(param_transformed = param_transformed, only_pi = only_pi))
  print(ll)
  return(ll)
}

calc_mle_derivatives_3waves_ar1_duration_endog_only_pi <- function(param_transformed, only_pi) {
  colSums(calc_lli_derivatives_3waves_ar1_duration_endog(param_transformed, only_pi = only_pi))
}



