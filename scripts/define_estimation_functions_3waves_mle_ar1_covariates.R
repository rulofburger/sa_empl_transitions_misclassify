# DEFINE FUNCTIONS ====

#> General functions that restrict parameter values to unit interval

# Transform parameter from (0,1) to (-Inf,Inf) using logit transformation
# This allows optimization over unconstrained parameters

logit_transform <- function(param0) {
  log(param0/(1 - param0))
}

# Transform parameter from (-Inf,Inf) to (0,1) using inverse logit (sigmoid)
# This maps optimized unconstrained parameters back to constrained space

logit_inverse <- function(param_input0) {
  1/(1 + exp(-param_input0))
}

# COVARIATE DEPENDENT TRANSITION RATES: EDUC + AGE ====


# Vectorized log-likelihood
calc_lli_3waves_ar1_covariates_age_educ <- function(param_transformed, pi0 = FALSE) {
  param <- param_transformed
  if (pi0) param$pi <- 0
  else     param$pi <- logit_inverse(param$pi)
  
  df_probs_temp <- fmutate(
    df_template_covariates_age_educ,
    
    # Time-specific transition rates
    theta0_1 = logit_inverse(param$intercept_0 + 
                               param$age_0*age1 +
                               param$age2_0*age1^2 + 
                               param$educ_0*educ1),
    theta0_2 = logit_inverse(param$intercept_0 +
                               param$age_0*age2 +
                               param$age2_0*age2^2 + 
                               param$educ_0*educ2),
    theta1_1 = logit_inverse(param$intercept_1 + 
                               param$age_1*age1 + 
                               param$age2_1*age1^2 + 
                               param$educ_1*educ1),
    theta1_2 = logit_inverse(param$intercept_1 +
                               param$age_1*age2 + 
                               param$age2_1*age2^2 + 
                               param$educ_1*educ2),
    
    mu_1 = theta0_1 / (theta0_1 + theta1_1),
    
    # Probabilities
    p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
    p2_star = fcase(
      y1_star == 0 & y2_star == 0, 1 - theta0_1,
      y1_star == 0 & y2_star == 1, theta0_1,
      y1_star == 1 & y2_star == 0, theta1_1,
      y1_star == 1 & y2_star == 1, 1 - theta1_1),
    p3_star = fcase(
      y2_star == 0 & y3_star == 0, 1 - theta0_2,
      y2_star == 0 & y3_star == 1, theta0_2,
      y2_star == 1 & y3_star == 0, theta1_2,
      y2_star == 1 & y3_star == 1, 1 - theta1_2),
    
    p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
    p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
    p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
    
    joint_p = p1_star * p1 * p2_star * p2 * p3_star * p3)
  
  df_probs <- df_probs_temp |>
    fgroup_by(y1, 
              y2, 
              y3, 
              age1, 
              age2, 
              educ1, 
              educ2) |>
    fsummarise(joint_p = fsum(joint_p)) |>
    fungroup()
  
  df_estimate |>
    join(df_probs, 
         on = c("y1", 
                "y2", 
                "y3", 
                "age1", 
                "age2", 
                "educ1", 
                "educ2"), 
         verbose = FALSE) |>
    fmutate(lli = weight * log(joint_p)) |>
    pull(lli) |> 
    as.numeric()
}

# Vectorized gradient
calc_lli_derivatives_3waves_ar1_covariates_age_educ <- function(param_transformed, pi0 = FALSE) {
  param <- param_transformed
  if (pi0) param$pi <- 0
  else     param$pi <- logit_inverse(param$pi)
  
  df_probs_temp <- df_template_covariates_age_educ |> 
    fmutate(
      theta0_1 = logit_inverse(param$intercept_0 +
                                 param$age_0*age1 + 
                                 param$age2_0*age1^2 +
                                 param$educ_0*educ1),
      theta0_2 = logit_inverse(param$intercept_0 + 
                                 param$age_0*age2 +
                                 param$age2_0*age2^2 + 
                                 param$educ_0*educ2),
      theta1_1 = logit_inverse(param$intercept_1 + 
                                 param$age_1*age1 +
                                 param$age2_1*age1^2 + 
                                 param$educ_1*educ1),
      theta1_2 = logit_inverse(param$intercept_1 + 
                                 param$age_1*age2 + 
                                 param$age2_1*age2^2 + 
                                 param$educ_1*educ2),
      
      mu_1 = theta0_1 / (theta0_1 + theta1_1),
      
      p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - theta0_1,
        y1_star == 0 & y2_star == 1, theta0_1,
        y1_star == 1 & y2_star == 0, theta1_1,
        y1_star == 1 & y2_star == 1, 1 - theta1_1),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - theta0_2,
        y2_star == 0 & y3_star == 1, theta0_2,
        y2_star == 1 & y3_star == 0, theta1_2,
        y2_star == 1 & y3_star == 1, 1 - theta1_2),
      
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      
      joint_p = p1_star * p1 * p2_star * p2 * p3_star * p3,
      
      # Derivatives w.r.t. transition rates
      d1_star_theta0_1 = fifelse(y1_star == 0, -theta1_1 / (theta1_1 + theta0_1)^2,
                                 theta1_1 / (theta1_1 + theta0_1)^2),
      d1_star_theta1_1 = fifelse(y1_star == 0,  theta0_1 / (theta1_1 + theta0_1)^2,
                                 -theta0_1 / (theta1_1 + theta0_1)^2),
      
      d2_star_theta0_1 = as.numeric(y1_star == 0 & y2_star == 0) * -1 +
        as.numeric(y1_star == 0 & y2_star == 1) * 1,
      d2_star_theta1_1 = as.numeric(y1_star == 1 & y2_star == 0) * 1 +
        as.numeric(y1_star == 1 & y2_star == 1) * -1,
      
      d3_star_theta0_2 = as.numeric(y2_star == 0 & y3_star == 0) * -1 +
        as.numeric(y2_star == 0 & y3_star == 1) * 1,
      d3_star_theta1_2 = as.numeric(y2_star == 1 & y3_star == 0) * 1 +
        as.numeric(y2_star == 1 & y3_star == 1) * -1,
      
      d1_pi = fifelse(y1 == y1_star, -1, 1),
      d2_pi = fifelse(y2 == y2_star, -1, 1),
      d3_pi = fifelse(y3 == y3_star, -1, 1),
      
      # Gradient components
      joint_d_theta0_1 = (d1_star_theta0_1 / p1_star + d2_star_theta0_1 / p2_star) * joint_p,
      joint_d_theta0_2 = d3_star_theta0_2 / p3_star * joint_p,
      joint_d_theta1_1 = (d1_star_theta1_1 / p1_star + d2_star_theta1_1 / p2_star) * joint_p,
      joint_d_theta1_2 = d3_star_theta1_2 / p3_star * joint_p,
      joint_d_pi       = (d1_pi / p1 + d2_pi / p2 + d3_pi / p3) * joint_p)
  
  df_grad <- df_probs_temp |>
    fgroup_by(y1, y2, y3, age1, age2, educ1, educ2) |>
    fsummarise(
      joint_d_theta0_1 = fsum(joint_d_theta0_1),
      joint_d_theta0_2 = fsum(joint_d_theta0_2),
      joint_d_theta1_1 = fsum(joint_d_theta1_1),
      joint_d_theta1_2 = fsum(joint_d_theta1_2),
      joint_d_pi       = fsum(joint_d_pi),
      joint_p          = fsum(joint_p),
      theta0_1         = fmean(theta0_1),
      theta0_2         = fmean(theta0_2),
      theta1_1         = fmean(theta1_1),
      theta1_2         = fmean(theta1_2)) |>
    fmutate(
      joint_d_theta0_1    = joint_d_theta0_1 * theta0_1 * (1 - theta0_1),
      joint_d_theta0_2    = joint_d_theta0_2 * theta0_2 * (1 - theta0_2),
      joint_d_theta1_1    = joint_d_theta1_1 * theta1_1 * (1 - theta1_1),
      joint_d_theta1_2    = joint_d_theta1_2 * theta1_2 * (1 - theta1_2),
      joint_d_pi          = joint_d_pi * param$pi * (1 - param$pi),
      
      # Chain rule for covariate-linked coefficients
      joint_d_intercept_0 = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_age_0       = joint_d_theta0_1 * age1 + joint_d_theta0_2 * age2,
      joint_d_age2_0      = joint_d_theta0_1 * age1^2 + joint_d_theta0_2 * age2^2,
      joint_d_educ_0      = joint_d_theta0_1 * educ1 + joint_d_theta0_2 * educ2,
      
      joint_d_intercept_1 = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_age_1       = joint_d_theta1_1 * age1 + joint_d_theta1_2 * age2,
      joint_d_age2_1      = joint_d_theta1_1 * age1^2 + joint_d_theta1_2 * age2^2,
      joint_d_educ_1      = joint_d_theta1_1 * educ1 + joint_d_theta1_2 * educ2)
  
  df_gi <- df_estimate |>
    join(df_grad, on = c("y1", "y2", "y3", "age1", "age2", "educ1", "educ2"), verbose = FALSE) |>
    fmutate(
      lgi_intercept_0 = weight * joint_d_intercept_0 / joint_p,
      lgi_age_0       = weight * joint_d_age_0       / joint_p,
      lgi_age2_0      = weight * joint_d_age2_0      / joint_p,
      lgi_educ_0      = weight * joint_d_educ_0      / joint_p,
      lgi_intercept_1 = weight * joint_d_intercept_1 / joint_p,
      lgi_age_1       = weight * joint_d_age_1       / joint_p,
      lgi_age2_1      = weight * joint_d_age2_1      / joint_p,
      lgi_educ_1      = weight * joint_d_educ_1      / joint_p,
      lgi_pi          = weight * joint_d_pi          / joint_p)
  
  if (pi0) {
    df_gi <- df_gi |> 
      fselect(lgi_intercept_0, 
              lgi_age_0, 
              lgi_age2_0, 
              lgi_educ_0,
              lgi_intercept_1, 
              lgi_age_1, 
              lgi_age2_1, 
              lgi_educ_1)
  } else {
    df_gi <- df_gi |> 
      fselect(lgi_intercept_0, 
              lgi_age_0, 
              lgi_age2_0, 
              lgi_educ_0,
              lgi_intercept_1, 
              lgi_age_1, 
              lgi_age2_1, 
              lgi_educ_1, 
              lgi_pi)
  }
  
  df_gi
  
}

# Wrappers
calc_mle_3waves_ar1_covariates_age_educ <- function(param_transformed) {
  fsum(calc_lli_3waves_ar1_covariates_age_educ(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates_age_educ <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_age_educ(param_transformed))
}

calc_mle_3waves_ar1_covariates_age_educ_pi0 <- function(param_transformed) {
  fsum(calc_lli_3waves_ar1_covariates_age_educ(param_transformed, pi0 = TRUE))
}

calc_mle_derivatives_3waves_ar1_covariates_age_educ_pi0 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_age_educ(param_transformed, pi0 = TRUE))
}





# COVARIATE DEPENDENT TRANSITION RATES: EDUC + AGE + RACE + FEMALE ====
calc_lli_3waves_ar1_covariates_age_educ_female_race <- function(param_transformed, pi0 = FALSE) {
  param <- param_transformed
  if (pi0) param$pi <- 0
  if (!pi0) param$pi <- 1 / (1 + exp(-param$pi))  # inverse logit
  
  # Copy base data
  df <- df_template_covariates_age_educ_female_race
  
  # Compute logits and probabilities
  lp00 <- with(df, 
               param$intercept_0 + param$age_0 * age1 + param$age2_0 * age1^2 +
                 param$educ_0 * educ1 + param$female_0 * female1 +
                 param$race2_0 * (race1 == 2) + param$race3_0 * (race1 == 3) + param$race4_0 * (race1 == 4))
  lp01 <- with(df, 
               param$intercept_0 + param$age_0 * age2 + param$age2_0 * age2^2 +
                 param$educ_0 * educ2 + param$female_0 * female2 +
                 param$race2_0 * (race2 == 2) + param$race3_0 * (race2 == 3) + param$race4_0 * (race2 == 4))
  lp10 <- with(df, 
               param$intercept_1 + param$age_1 * age1 + param$age2_1 * age1^2 +
                 param$educ_1 * educ1 + param$female_1 * female1 +
                 param$race2_1 * (race1 == 2) + param$race3_1 * (race1 == 3) + param$race4_1 * (race1 == 4))
  lp11 <- with(df, 
               param$intercept_1 + param$age_1 * age2 + param$age2_1 * age2^2 +
                 param$educ_1 * educ2 + param$female_1 * female2 +
                 param$race2_1 * (race2 == 2) + param$race3_1 * (race2 == 3) + param$race4_1 * (race2 == 4))
  
  # Inverse logit
  th00 <- 1 / (1 + exp(-lp00))
  th01 <- 1 / (1 + exp(-lp01))
  th10 <- 1 / (1 + exp(-lp10))
  th11 <- 1 / (1 + exp(-lp11))
  
  mu1 <- th00 / (th00 + th10)
  
  # Index-based conditional probabilities
  idx12 <- df$y1_star * 2 + df$y2_star + 1L
  idx23 <- df$y2_star * 2 + df$y3_star + 1L
  
  p2_star_mat <- cbind(1 - th00, th00, th10, 1 - th10)
  p3_star_mat <- cbind(1 - th01, th01, th11, 1 - th11)
  
  n <- nrow(df)
  
  # Compute full probabilities
  p1_star <- fifelse(df$y1_star == 1, mu1, 1 - mu1)
  p2_star <- p2_star_mat[cbind(seq_len(n), idx12)]
  p3_star <- p3_star_mat[cbind(seq_len(n), idx23)]
  
  p1 <- fifelse(df$y1 == df$y1_star, 1 - param$pi, param$pi)
  p2 <- fifelse(df$y2 == df$y2_star, 1 - param$pi, param$pi)
  p3 <- fifelse(df$y3 == df$y3_star, 1 - param$pi, param$pi)
  
  joint_p <- p1_star * p1 * p2_star * p2 * p3_star * p3
  
  # Collapse group summary
  df_probs <- df |> 
    fgroup_by(y1, 
              y2,
              y3, 
              age1, 
              age2, 
              educ1, 
              educ2, 
              female1, 
              female2, 
              race1, 
              race2) |>
    fsummarise(joint_p = fsum(joint_p)) |>
    fungroup()
  
  # Final log-likelihood contributions
  lli <- 
    join(df_estimate, 
         df_probs,
         on       = c("y1", 
                      "y2", 
                      "y3", 
                      "age1",
                      "age2", 
                      "educ1",
                      "educ2", 
                      "female1", 
                      "female2", 
                      "race1", 
                      "race2"),
         verbose  = FALSE, 
         validate = F) |> 
  fmutate(lli = weight * log(joint_p)) |> 
    pull(lli)
  
  lli
}



calc_lli_derivatives_3waves_ar1_covariates_age_educ_female_race <- function(param_transformed, pi0 = FALSE) {
  
  
  
  param <- param_transformed
  if (pi0) param$pi <- 0
  if (!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates_age_educ_female_race %>% 
    fmutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$female_0*female1 + param$race2_0*(race1 == 2) + param$race3_0*(race1 == 3) + param$race4_0*(race1 == 4)),
      theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$female_0*female2 + param$race2_0*(race2 == 2) + param$race3_0*(race2 == 3) + param$race4_0*(race2 == 4)),
      theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$female_1*female1 + param$race2_1*(race1 == 2) + param$race3_1*(race1 == 3) + param$race4_1*(race1 == 4)),
      theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$female_1*female2 + param$race2_1*(race2 == 2) + param$race3_1*(race2 == 3) + param$race4_1*(race2 == 4)),
      mu_1     = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - theta0_1,
        y1_star == 0 & y2_star == 1, theta0_1,
        y1_star == 1 & y2_star == 0, theta1_1,
        y1_star == 1 & y2_star == 1, 1 - theta1_1
      ),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - theta0_2,
        y2_star == 0 & y3_star == 1, theta0_2,
        y2_star == 1 & y3_star == 0, theta1_2,
        y2_star == 1 & y3_star == 1, 1 - theta1_2
      ),
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) %>% 
    fmutate(
      d1_star_theta0_1 = fcase(
        y1_star == 0,  -theta1_1/((theta1_1 + theta0_1)^2),
        y1_star == 1,  theta1_1/((theta1_1 + theta0_1)^2)
      ),
      d1_star_theta1_1 = fcase(
        y1_star == 0, theta0_1/((theta1_1 + theta0_1)^2),
        y1_star == 1, -theta0_1/((theta1_1 + theta0_1)^2)
      ),
      d2_star_theta0_1 = fcase(
        y1_star == 0 & y2_star == 0, -1,
        y1_star == 0 & y2_star == 1, 1,
        y1_star == 1 & y2_star == 0, 0,
        y1_star == 1 & y2_star == 1, 0
      ),
      d2_star_theta1_1 = fcase(
        y1_star == 0 & y2_star == 0, 0,
        y1_star == 0 & y2_star == 1, 0,
        y1_star == 1 & y2_star == 0, 1,
        y1_star == 1 & y2_star == 1, -1
      ),
      d3_star_theta0_2 = fcase(
        y2_star == 0 & y3_star == 0, -1,
        y2_star == 0 & y3_star == 1, 1,
        y2_star == 1 & y3_star == 0, 0,
        y2_star == 1 & y3_star == 1, 0
      ),
      d3_star_theta1_2 = fcase(
        y2_star == 0 & y3_star == 0, 0,
        y2_star == 0 & y3_star == 1, 0,
        y2_star == 1 & y3_star == 0, 1,
        y2_star == 1 & y3_star == 1, -1
      ),
      d1_pi = fifelse(y1 == y1_star, -1, 1),
      d2_pi = fifelse(y2 == y2_star, -1, 1),
      d3_pi = fifelse(y3 == y3_star, -1, 1),
      joint_d_theta0_1 = 
        d1_star_theta0_1*joint_p/p1_star + 
        d2_star_theta0_1*joint_p/p2_star,
      joint_d_theta0_2 = 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star,
      joint_d_theta1_2 = 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    fgroup_by(y1, y2, y3, age1, age2, educ1, educ2, female1, female2, race1, race2) %>% 
    fsummarise(
      joint_d_theta0_1 = fsum(joint_d_theta0_1),
      joint_d_theta0_2 = fsum(joint_d_theta0_2),
      joint_d_theta1_1 = fsum(joint_d_theta1_1),
      joint_d_theta1_2 = fsum(joint_d_theta1_2), 
      joint_d_pi       = fsum(joint_d_pi),
      joint_p          = fsum(joint_p),
      theta0_1         = fmean(theta0_1),
      theta0_2         = fmean(theta0_2),
      theta1_1         = fmean(theta1_1),
      theta1_2         = fmean(theta1_2)) %>% 
    fungroup() |> 
    fmutate(
      joint_d_theta0_1    = joint_d_theta0_1*theta0_1*(1 - theta0_1), 
      joint_d_theta0_2    = joint_d_theta0_2*theta0_2*(1 - theta0_2),
      joint_d_theta1_1    = joint_d_theta1_1*theta1_1*(1 - theta1_1), 
      joint_d_theta1_2    = joint_d_theta1_2*theta1_2*(1 - theta1_2),
      joint_d_pi          = joint_d_pi*param$pi*(1 - param$pi), 
      joint_d_intercept_0 = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_age_0       = joint_d_theta0_1*age1 + joint_d_theta0_2*age2,
      joint_d_age2_0      = joint_d_theta0_1*(age1^2) + joint_d_theta0_2*(age2^2),
      joint_d_educ_0      = joint_d_theta0_1*educ1 + joint_d_theta0_2*educ2,
      joint_d_female_0    = joint_d_theta0_1*female1 + joint_d_theta0_2*female2,
      joint_d_race2_0     = joint_d_theta0_1*(race1 == 2) + joint_d_theta0_2*(race2 == 2),
      joint_d_race3_0     = joint_d_theta0_1*(race1 == 3) + joint_d_theta0_2*(race2 == 3),
      joint_d_race4_0     = joint_d_theta0_1*(race1 == 4) + joint_d_theta0_2*(race2 == 4),
      joint_d_intercept_1 = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_age_1       = joint_d_theta1_1*age1 + joint_d_theta1_2*age2,
      joint_d_age2_1      = joint_d_theta1_1*(age1^2) + joint_d_theta1_2*(age2^2),
      joint_d_educ_1      = joint_d_theta1_1*educ1 + joint_d_theta1_2*educ2,
      joint_d_female_1    = joint_d_theta1_1*female1 + joint_d_theta1_2*female2,
      joint_d_race2_1     = joint_d_theta1_1*(race1 == 2) + joint_d_theta1_2*(race2 == 2),
      joint_d_race3_1     = joint_d_theta1_1*(race1 == 3) + joint_d_theta1_2*(race2 == 3),
      joint_d_race4_1     = joint_d_theta1_1*(race1 == 4) + joint_d_theta1_2*(race2 == 4))
  
  df_gi <- df_estimate %>% 
    join(df_grad, 
         how = "left",
                    on           = c('y1', 
                                     'y2', 
                                     'y3', 
                                     'age1', 
                                     'age2', 
                                     'educ1', 
                                     'educ2', 
                                     'female1', 
                                     'female2', 
                                     'race1', 
                                     'race2'), 
                     verbose      = F) %>% 
    fmutate(
      lgi_intercept_0 = weight*joint_d_intercept_0/joint_p,
      lgi_age_0       = weight*joint_d_age_0/joint_p,
      lgi_age2_0      = weight*joint_d_age2_0/joint_p,
      lgi_educ_0      = weight*joint_d_educ_0/joint_p,
      lgi_female_0    = weight*joint_d_female_0/joint_p,
      lgi_race2_0     = weight*joint_d_race2_0/joint_p,
      lgi_race3_0     = weight*joint_d_race3_0/joint_p,
      lgi_race4_0     = weight*joint_d_race4_0/joint_p,
      lgi_intercept_1 = weight*joint_d_intercept_1/joint_p,
      lgi_age_1       = weight*joint_d_age_1/joint_p,
      lgi_age2_1      = weight*joint_d_age2_1/joint_p,
      lgi_educ_1      = weight*joint_d_educ_1/joint_p,
      lgi_female_1    = weight*joint_d_female_1/joint_p,
      lgi_race2_1     = weight*joint_d_race2_1/joint_p,
      lgi_race3_1     = weight*joint_d_race3_1/joint_p,
      lgi_race4_1     = weight*joint_d_race4_1/joint_p,
      lgi_pi          = weight*joint_d_pi/joint_p) %>%
    fselect(lgi_intercept_0, 
           lgi_age_0, 
           lgi_age2_0, 
           lgi_educ_0, 
           lgi_female_0, 
           lgi_race2_0, 
           lgi_race3_0, 
           lgi_race4_0, 
           lgi_intercept_1, 
           lgi_age_1, 
           lgi_age2_1,
           lgi_educ_1, 
           lgi_female_1,
           lgi_race2_1, 
           lgi_race3_1, 
           lgi_race4_1, 
           lgi_pi)
  
  if (pi0) df_gi <- df_gi %>% fselect(-lgi_pi)
  
  df_gi
 
}




calc_mle_3waves_ar1_covariates_age_educ_female_race <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates_age_educ_female_race(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_age_educ_female_race(param_transformed))
}

calc_mle_3waves_ar1_covariates_age_educ_female_race_pi0 <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates_age_educ_female_race(param_transformed, pi0 = TRUE))
}

calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_pi0 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_age_educ_female_race(param_transformed, pi0 = TRUE))
}





# COVARIATE DEPENDENT TRANSITION RATES: EDUC + AGE + RACE + FEMALE + CONTRACT ====
calc_lli_3waves_ar1_covariates_age_educ_female_race_contract <- function(param_transformed, pi0 = FALSE) {
  param    <- param_transformed
  param$pi <- if (pi0) 0 else logit_inverse(param$pi)
  
  df <- df_template_covariates_age_educ_female_race_contract
  
  # Precompute linear predictors
  linpred <- ftransform(
    df,
    lp00 = param$intercept_0 + 
      param$age_0 * age1 + 
      param$age2_0 * age1^2 + 
      param$educ_0 * educ1 +
      param$female_0 * female1 +
      param$race2_0 * (race1 == 2) + 
      param$race3_0 * (race1 == 3) + 
      param$race4_0 * (race1 == 4),
    lp01 = param$intercept_0 + 
      param$age_0 * age2 + 
      param$age2_0 * age2^2 + 
      param$educ_0 * educ2 +
      param$female_0 * female2 + 
      param$race2_0 * (race2 == 2) +
      param$race3_0 * (race2 == 3) + 
      param$race4_0 * (race2 == 4),
    lp10 = param$intercept_1 + 
      param$age_1 * age1 + 
      param$age2_1 * age1^2 + 
      param$educ_1 * educ1 +
      param$female_1 * female1 + 
      param$race2_1 * (race1 == 2) + 
      param$race3_1 * (race1 == 3) + 
      param$race4_1 * (race1 == 4) +
      param$contract * contracttype1 + 
      param$contract * param$contract_missing * contracttype1_missing,
    lp11 = param$intercept_1 + 
      param$age_1 * age2 + 
      param$age2_1 * age2^2 + 
      param$educ_1 * educ2 +
      param$female_1 * female2 + 
      param$race2_1 * (race2 == 2) + 
      param$race3_1 * (race2 == 3) + 
      param$race4_1 * (race2 == 4) +
      param$contract * contracttype2 + 
      param$contract * param$contract_missing * contracttype2_missing)
  
  # Compute transition probabilities
  linpred <- linpred |> 
    fmutate(
      th00 = logit_inverse(lp00),
      th01 = logit_inverse(lp01),
      th10 = logit_inverse(lp10),
      th11 = logit_inverse(lp11),
      mu1  = th00 / (th00 + th10))
  
  # Indexed transitions
  idx12 <- linpred$y1_star * 2 + linpred$y2_star + 1L
  idx23 <- linpred$y2_star * 2 + linpred$y3_star + 1L
  
  mat_p2 <- cbind(1 - linpred$th00, 
                  linpred$th00, 
                  linpred$th10, 
                  1 - linpred$th10)
  mat_p3 <- cbind(1 - linpred$th01, 
                  linpred$th01, 
                  linpred$th11, 
                  1 - linpred$th11)
  
  df_probs_temp <- linpred |> 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu1, 1 - mu1),
      p2_star = mat_p2[cbind(seq_len(nrow(df)), idx12)],
      p3_star = mat_p3[cbind(seq_len(nrow(df)), idx23)],
      p1      = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2      = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3      = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star * p1 * p2_star * p2 * p3_star * p3)
  
  # Collapse-style group summarise
  df_probs <- df_probs_temp |> 
    fgroup_by(y1, 
              y2, 
              y3, 
              age1, 
              age2, 
              educ1, 
              educ2, 
              female1, 
              female2, 
              race1, 
              race2,
              contracttype1, 
              contracttype2, 
              contracttype1_missing, 
              contracttype2_missing) |> 
    fsummarise(joint_p = fsum(joint_p)) |> 
    fungroup()
  
  # Final log-likelihood vector
  lli <- 
    join(df_estimate, 
       df_probs,
       on = c("y1", 
              "y2", 
              "y3", 
              "age1",
              "age2", 
              "educ1", 
              "educ2", 
              "female1", 
              "female2",
              "race1", 
              "race2", 
              "contracttype1", 
              "contracttype2", 
              "contracttype1_missing", 
              "contracttype2_missing"),
       verbose  = FALSE, 
       validate = F) |> 
  ftransform(
    lli = weight * log(joint_p)) |> 
    pull(lli)
  
  lli
}



calc_lli_derivatives_3waves_ar1_covariates_age_educ_female_race_contract <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  param$pi <- if (pi0) 0 else logit_inverse(param$pi)
  
  df <- df_template_covariates_age_educ_female_race_contract
  
  # Precompute linear predictors
  linpred <- df |> 
    ftransform(
    lp00 = param$intercept_0 + 
      param$age_0 * age1 + 
      param$age2_0 * age1^2 + 
      param$educ_0 * educ1 +
      param$female_0 * female1 + 
      param$race2_0 * (race1 == 2) + 
      param$race3_0 * (race1 == 3) + 
      param$race4_0 * (race1 == 4),
    lp01 = param$intercept_0 + 
      param$age_0 * age2 + 
      param$age2_0 * age2^2 + 
      param$educ_0 * educ2 +
      param$female_0 * female2 + 
      param$race2_0 * (race2 == 2) +
      param$race3_0 * (race2 == 3) +
      param$race4_0 * (race2 == 4),
    lp10 = param$intercept_1 + 
      param$age_1 * age1 + 
      param$age2_1 * age1^2 + 
      param$educ_1 * educ1 +
      param$female_1 * female1 + 
      param$race2_1 * (race1 == 2) + 
      param$race3_1 * (race1 == 3) + 
      param$race4_1 * (race1 == 4) +
      param$contract * contracttype1 + 
      param$contract * param$contract_missing * contracttype1_missing,
    lp11 = param$intercept_1 + 
      param$age_1 * age2 + 
      param$age2_1 * age2^2 + 
      param$educ_1 * educ2 +
      param$female_1 * female2 + 
      param$race2_1 * (race2 == 2) +
      param$race3_1 * (race2 == 3) + 
      param$race4_1 * (race2 == 4) +
      param$contract * contracttype2 + 
      param$contract * param$contract_missing * contracttype2_missing
  )
  
  # Transition probabilities
  linpred <- linpred |> 
    fmutate(
      th00 = logit_inverse(lp00),
      th01 = logit_inverse(lp01),
      th10 = logit_inverse(lp10),
      th11 = logit_inverse(lp11),
      mu1  = th00 / (th00 + th10))
  
  # Indices for lookup
  idx12 <- linpred$y1_star * 2 + linpred$y2_star + 1L
  idx23 <- linpred$y2_star * 2 + linpred$y3_star + 1L
  
  mat_p2 <- cbind(1 - linpred$th00, 
                  linpred$th00, 
                  linpred$th10, 
                  1 - linpred$th10)
  mat_p3 <- cbind(1 - linpred$th01,
                  linpred$th01, 
                  linpred$th11, 
                  1 - linpred$th11)
  
  # Likelihood components
  df_probs_temp <- linpred |> 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu1, 1 - mu1),
      p2_star = mat_p2[cbind(seq_len(nrow(df)), idx12)],
      p3_star = mat_p3[cbind(seq_len(nrow(df)), idx23)],
      p1      = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2      = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3      = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star * p1 * p2_star * p2 * p3_star * p3)
  
  # Derivatives of latent probs
  df_probs_temp <- df_probs_temp |> 
    fmutate(
      d1_star_theta0_1 = (2 * y1_star - 1) * th10 / (th00 + th10)^2,
      d1_star_theta1_1 = -(2 * y1_star - 1) * th00 / (th00 + th10)^2,
      d2_star_theta0_1 = (y1_star == 0) * (2 * y2_star - 1),
      d2_star_theta1_1 = (y1_star == 1) * (1 - 2 * y2_star),
      d3_star_theta0_2 = (y2_star == 0) * (2 * y3_star - 1),
      d3_star_theta1_2 = (y2_star == 1) * (1 - 2 * y3_star),
      d1_pi            = ifelse(y1 == y1_star, -1, 1),
      d2_pi            = ifelse(y2 == y2_star, -1, 1),
      d3_pi            = ifelse(y3 == y3_star, -1, 1),
      joint_d_theta0_1 = d1_star_theta0_1 * joint_p / p1_star + d2_star_theta0_1 * joint_p / p2_star,
      joint_d_theta0_2 = d3_star_theta0_2 * joint_p / p3_star,
      joint_d_theta1_1 = d1_star_theta1_1 * joint_p / p1_star + d2_star_theta1_1 * joint_p / p2_star,
      joint_d_theta1_2 = d3_star_theta1_2 * joint_p / p3_star,
      joint_d_pi       = d1_pi * joint_p / p1 + d2_pi * joint_p / p2 + d3_pi * joint_p / p3)
  
  # Aggregate derivatives
  df_grad <- df_probs_temp |> 
    fgroup_by(
      y1, 
      y2, 
      y3, 
      age1, 
      age2, 
      educ1, 
      educ2, 
      female1, 
      female2, 
      race1, 
      race2,
      contracttype1, 
      contracttype2, 
      contracttype1_missing, 
      contracttype2_missing) |> 
    fsummarise(
      joint_d_theta0_1 = fsum(joint_d_theta0_1),
      joint_d_theta0_2 = fsum(joint_d_theta0_2),
      joint_d_theta1_1 = fsum(joint_d_theta1_1),
      joint_d_theta1_2 = fsum(joint_d_theta1_2),
      joint_d_pi       = fsum(joint_d_pi),
      joint_p          = fsum(joint_p),
      th00             = fmean(th00), 
      th01             = fmean(th01), 
      th10             = fmean(th10), 
      th11             = fmean(th11)) |> 
    ftransform(
      joint_d_theta0_1 = joint_d_theta0_1 * th00 * (1 - th00),
      joint_d_theta0_2 = joint_d_theta0_2 * th01 * (1 - th01),
      joint_d_theta1_1 = joint_d_theta1_1 * th10 * (1 - th10),
      joint_d_theta1_2 = joint_d_theta1_2 * th11 * (1 - th11),
      joint_d_pi       = joint_d_pi * param$pi * (1 - param$pi))
  
  # Compute gradients per parameter
  df_grad <- df_grad |> 
    ftransform(
      joint_d_intercept_0      = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_age_0            = joint_d_theta0_1 * age1 + joint_d_theta0_2 * age2,
      joint_d_age2_0           = joint_d_theta0_1 * age1^2 + joint_d_theta0_2 * age2^2,
      joint_d_educ_0           = joint_d_theta0_1 * educ1 + joint_d_theta0_2 * educ2,
      joint_d_female_0         = joint_d_theta0_1 * female1 + joint_d_theta0_2 * female2,
      joint_d_race2_0          = joint_d_theta0_1 * (race1 == 2) + joint_d_theta0_2 * (race2 == 2),
      joint_d_race3_0          = joint_d_theta0_1 * (race1 == 3) + joint_d_theta0_2 * (race2 == 3),
      joint_d_race4_0          = joint_d_theta0_1 * (race1 == 4) + joint_d_theta0_2 * (race2 == 4),
      joint_d_intercept_1      = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_age_1            = joint_d_theta1_1 * age1 + joint_d_theta1_2 * age2,
      joint_d_age2_1           = joint_d_theta1_1 * age1^2 + joint_d_theta1_2 * age2^2,
      joint_d_educ_1           = joint_d_theta1_1 * educ1 + joint_d_theta1_2 * educ2,
      joint_d_female_1         = joint_d_theta1_1 * female1 + joint_d_theta1_2 * female2,
      joint_d_race2_1          = joint_d_theta1_1 * (race1 == 2) + joint_d_theta1_2 * (race2 == 2),
      joint_d_race3_1          = joint_d_theta1_1 * (race1 == 3) + joint_d_theta1_2 * (race2 == 3),
      joint_d_race4_1          = joint_d_theta1_1 * (race1 == 4) + joint_d_theta1_2 * (race2 == 4),
      joint_d_contract         = joint_d_theta1_1 * contracttype1 + joint_d_theta1_2 * contracttype2 +
        param$contract_missing * (joint_d_theta1_1 * contracttype1_missing + joint_d_theta1_2 * contracttype2_missing),
      joint_d_contract_missing = param$contract * (joint_d_theta1_1 * contracttype1_missing + joint_d_theta1_2 * contracttype2_missing))
  
  # Join with estimate df and compute log-gradient contributions
  df_gi <- join(df_estimate, 
                df_grad,
                on = c("y1", 
                       "y2", 
                       "y3", 
                       "age1", 
                       "age2", 
                       "educ1", 
                       "educ2", 
                       "female1", 
                       "female2",
                       "race1", 
                       "race2", 
                       "contracttype1",
                       "contracttype2",
                       "contracttype1_missing", 
                       "contracttype2_missing"),
                verbose  = FALSE, 
                validate = F) |>
    ftransform(
      lgi_intercept_0      = weight * joint_d_intercept_0 / joint_p,
      lgi_age_0            = weight * joint_d_age_0 / joint_p,
      lgi_age2_0           = weight * joint_d_age2_0 / joint_p,
      lgi_educ_0           = weight * joint_d_educ_0 / joint_p,
      lgi_female_0         = weight * joint_d_female_0 / joint_p,
      lgi_race2_0          = weight * joint_d_race2_0 / joint_p,
      lgi_race3_0          = weight * joint_d_race3_0 / joint_p,
      lgi_race4_0          = weight * joint_d_race4_0 / joint_p,
      lgi_intercept_1      = weight * joint_d_intercept_1 / joint_p,
      lgi_age_1            = weight * joint_d_age_1 / joint_p,
      lgi_age2_1           = weight * joint_d_age2_1 / joint_p,
      lgi_educ_1           = weight * joint_d_educ_1 / joint_p,
      lgi_female_1         = weight * joint_d_female_1 / joint_p,
      lgi_race2_1          = weight * joint_d_race2_1 / joint_p,
      lgi_race3_1          = weight * joint_d_race3_1 / joint_p,
      lgi_race4_1          = weight * joint_d_race4_1 / joint_p,
      lgi_contract         = weight * joint_d_contract / joint_p,
      lgi_contract_missing = weight * joint_d_contract_missing / joint_p,
      lgi_pi               = weight * joint_d_pi / joint_p)
  
  # Keep only the desired columns that start with 'lgi_'
  df_gi <- df_gi[, startsWith(names(df_gi), "lgi_"), drop = FALSE]
  
  # Remove pi derivative if pi0 = TRUE
  if (pi0) df_gi <- df_gi[, names(df_gi) != "lgi_pi", drop = FALSE]
  

  return(df_gi)
}

calc_mle_3waves_ar1_covariates_age_educ_female_race_contract <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates_age_educ_female_race_contract(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_contract <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_age_educ_female_race_contract(param_transformed))
}

calc_mle_3waves_ar1_covariates_age_educ_female_race_contract_pi0 <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates_age_educ_female_race_contract(param_transformed, pi0 = TRUE))
}

calc_mle_derivatives_3waves_ar1_covariates_age_educ_female_race_contract_pi0 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_age_educ_female_race_contract(param_transformed, pi0 = TRUE))
}


# DURATION DEPENDENCE (timegap + tenure) ====


calc_lli_3waves_ar1_covariates_duration <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    fmutate(
      theta0_1 = logit_inverse(param$intercept_0 + (y1 == 0)*param$timegap*log(timegap1 + 0.125) + (y1 == 0)*param$timegap*log(param$timegap_missing + 0.125)),
      theta0_2 = logit_inverse(param$intercept_0 + (y2 == 0)*param$timegap*log(timegap2 + 0.125) + (y2 == 0)*param$timegap*log(param$timegap_missing + 0.125)),
      theta1_1 = logit_inverse(param$intercept_1 + (y1 == 1)*param$tenure*log(tenure1 + 0.125) + (y1 == 1)*param$tenure*log(param$tenure_missing + 0.125)),
      theta1_2 = logit_inverse(param$intercept_1 + (y2 == 1)*param$tenure*log(tenure2 + 0.125) + (y2 == 1)*param$tenure*log(param$tenure_missing + 0.125)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - theta0_1,
        y1_star == 0 & y2_star == 1, theta0_1,
        y1_star == 1 & y2_star == 0, theta1_1,
        y1_star == 1 & y2_star == 1, 1 - theta1_1
      ),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - theta0_2,
        y2_star == 0 & y3_star == 1, theta0_2,
        y2_star == 1 & y3_star == 0, theta1_2,
        y2_star == 1 & y3_star == 1, 1 - theta1_2
      ),
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 

  df_probs <- df_probs_temp %>% 
    fgroup_by(y1, y2, y3, timegap1, timegap2, tenure1, tenure2) %>% 
    fsummarise(joint_p = sum(joint_p)) |> 
    fungroup()
 
  
  lli <- df_estimate %>% 
    join(df_probs, 
         on       = c('y1', 
                      'y2', 
                      'y3', 
                      'timegap1', 
                      'timegap2', 
                      'tenure1', 
                      'tenure2'), 
         how      = "left", 
         validate = F, 
         verbose = F) %>% 
    fmutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
  lli
  
}

calc_lli_derivatives_3waves_ar1_covariates_duration <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    fmutate(
      theta0_1 = logit_inverse(param$intercept_0 + (y1 == 0)*param$timegap*log(timegap1 + 0.125) + (y1 == 0)*param$timegap*log(param$timegap_missing + 0.125)),
      theta0_2 = logit_inverse(param$intercept_0 + (y2 == 0)*param$timegap*log(timegap2 + 0.125) + (y2 == 0)*param$timegap*log(param$timegap_missing + 0.125)),
      theta1_1 = logit_inverse(param$intercept_1 + (y1 == 1)*param$tenure*log(tenure1 + 0.125) + (y1 == 1)*param$tenure*log(param$tenure_missing + 0.125)),
      theta1_2 = logit_inverse(param$intercept_1 + (y2 == 1)*param$tenure*log(tenure2 + 0.125) + (y2 == 1)*param$tenure*log(param$tenure_missing + 0.125)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)) %>% 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - theta0_1,
        y1_star == 0 & y2_star == 1, theta0_1,
        y1_star == 1 & y2_star == 0, theta1_1,
        y1_star == 1 & y2_star == 1, 1 - theta1_1
      ),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - theta0_2,
        y2_star == 0 & y3_star == 1, theta0_2,
        y2_star == 1 & y3_star == 0, theta1_2,
        y2_star == 1 & y3_star == 1, 1 - theta1_2
      ),
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) %>% 
    fmutate(
      d1_star_theta0_1 = fcase(
        y1_star == 0, -theta1_1/((theta1_1 + theta0_1)^2),
        y1_star == 1, theta1_1/((theta1_1 + theta0_1)^2)
      ),
      d1_star_theta1_1 = fcase(
        y1_star == 0, theta0_1/((theta1_1 + theta0_1)^2),
        y1_star == 1, -theta0_1/((theta1_1 + theta0_1)^2)
      ),
      d2_star_theta0_1 = fcase(
        y1_star == 0 & y2_star == 0, -1,
        y1_star == 0 & y2_star == 1, 1,
        y1_star == 1 & y2_star == 0, 0,
        y1_star == 1 & y2_star == 1, 0
      ),
      d2_star_theta1_1 = fcase(
        y1_star == 0 & y2_star == 0, 0,
        y1_star == 0 & y2_star == 1, 0,
        y1_star == 1 & y2_star == 0, 1,
        y1_star == 1 & y2_star == 1, -1
      ),
      d3_star_theta0_2 = fcase(
        y2_star == 0 & y3_star == 0, -1,
        y2_star == 0 & y3_star == 1, 1,
        y2_star == 1 & y3_star == 0, 0,
        y2_star == 1 & y3_star == 1, 0
      ),
      d3_star_theta1_2 = fcase(
        y2_star == 0 & y3_star == 0, 0,
        y2_star == 0 & y3_star == 1, 0,
        y2_star == 1 & y3_star == 0, 1,
        y2_star == 1 & y3_star == 1, -1
      ),
      d1_pi = fifelse(y1 == y1_star, -1, 1),
      d2_pi = fifelse(y2 == y2_star, -1, 1),
      d3_pi = fifelse(y3 == y3_star, -1, 1),
      joint_d_theta0_1 = 
        d1_star_theta0_1*joint_p/p1_star + 
        d2_star_theta0_1*joint_p/p2_star,
      joint_d_theta0_2 = 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star,
      joint_d_theta1_2 = 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3) 
  
  df_grad <- df_probs_temp %>% 
    fmutate(
      joint_d_timegap0_1 = joint_d_theta0_1*log(timegap1 + 0.125),
      joint_d_timegap0_2 = joint_d_theta0_2*log(timegap2 + 0.125),
      joint_d_tenure1_1  = joint_d_theta1_1*log(tenure1 + 0.125),
      joint_d_tenure1_2  = joint_d_theta1_2*log(tenure2 + 0.125)) %>% 
    fgroup_by(y1, y2, y3, timegap1, timegap2, tenure1, tenure2) %>% 
    fsummarise(
      joint_d_theta0_1   = sum(joint_d_theta0_1),
      joint_d_theta0_2   = sum(joint_d_theta0_2),
      joint_d_theta1_1   = sum(joint_d_theta1_1),
      joint_d_theta1_2   = sum(joint_d_theta1_2),
      joint_d_timegap0_1 = sum(joint_d_timegap0_1),
      joint_d_timegap0_2 = sum(joint_d_timegap0_2),
      joint_d_tenure1_1  = sum(joint_d_tenure1_1),
      joint_d_tenure1_2  = sum(joint_d_tenure1_2),
      joint_d_pi         = sum(joint_d_pi),
      joint_p            = sum(joint_p),
      theta0_1           = mean(theta0_1),
      theta0_2           = mean(theta0_2),
      theta1_1           = mean(theta1_1),
      theta1_2           = mean(theta1_2)) %>%
    fungroup() |> 
    fmutate(
      joint_d_theta0_1   = joint_d_theta0_1*theta0_1*(1 - theta0_1), 
      joint_d_theta0_2   = joint_d_theta0_2*theta0_2*(1 - theta0_2),
      joint_d_theta1_1   = joint_d_theta1_1*theta1_1*(1 - theta1_1), 
      joint_d_theta1_2   = joint_d_theta1_2*theta1_2*(1 - theta1_2),
      joint_d_timegap0_1 = joint_d_timegap0_1*theta0_1*(1 - theta0_1),
      joint_d_timegap0_2 = joint_d_timegap0_2*theta0_2*(1 - theta0_2),
      joint_d_tenure1_1  = joint_d_tenure1_1*theta1_1*(1 - theta1_1),
      joint_d_tenure1_2  = joint_d_tenure1_2*theta1_2*(1 - theta1_2),
      joint_d_pi         = joint_d_pi*param$pi*(1 - param$pi)) %>%
    fmutate(
      joint_d_intercept_0 = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_timegap_0   = joint_d_timegap0_1 + joint_d_timegap0_2,
      joint_d_intercept_1 = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_tenure_1    = joint_d_tenure1_1 + joint_d_tenure1_2) %>% 
    frename(
      timegap1 = timegap1_old,
      timegap2 = timegap2_old,
      tenure1  = tenure1_old,
      tenure2  = tenure2_old)
  
  df_gi <- df_estimate %>% 
    join(df_grad,
         on = c('y1', 
                'y2', 
                'y3', 
                'timegap1', 
                'timegap2', 
                'tenure1', 
                'tenure2'), 
         validate = F, 
         verbose = F) %>% 
    fmutate(
      lgi_intercept_0 = weight*joint_d_intercept_0/joint_p,
      lgi_timegap_0   = weight*joint_d_timegap_0/joint_p,
      lgi_intercept_1 = weight*joint_d_intercept_1/joint_p,
      lgi_tenure_1    = weight*joint_d_tenure_1/joint_p,
      lgi_pi          = weight*joint_d_pi/joint_p) %>%
    fselect(lgi_intercept_0, 
            lgi_timegap_0, 
            lgi_intercept_1, 
            lgi_tenure_1, 
            lgi_pi)
  
  df_gi
}

calc_mle_3waves_ar1_covariates_duration <- function(param_transformed) {
  
  sum(calc_lli_3waves_ar1_covariates_duration(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates_duration <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_duration(param_transformed))
}

calc_mle_3waves_ar1_covariates2_pi0 <- function(param_transformed) {
  sum(calc_lli_3waves_ar1_covariates_duration(param_transformed, pi0 = TRUE))
}

calc_mle_derivatives_3waves_ar1_covariates_duration_pi0 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates_duration(param_transformed, pi0 = TRUE))
}

# FUNCTIONS FOR AR(1) ML ESTIMATOR OVER 3 WAVES WITH MISCLASSIFICATION ERROR AND COVARIATE DEPENDENT TRANSITION RATES (educ + age + timegap + tenure) ====


calc_lli_3waves_ar1_covariates3 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    fmutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$timegap*log(timegap1 + 1.5)),
      theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$timegap*log(timegap2 + 1.5)),
      theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$tenure*log(tenure1 + 1.5)),
      theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$tenure*log(tenure2 + 1.5)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - theta0_1,
        y1_star == 0 & y2_star == 1, theta0_1,
        y1_star == 1 & y2_star == 0, theta1_1,
        y1_star == 1 & y2_star == 1, 1 - theta1_1
      ),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - theta0_2,
        y2_star == 0 & y3_star == 1, theta0_2,
        y2_star == 1 & y3_star == 0, theta1_2,
        y2_star == 1 & y3_star == 1, 1 - theta1_2
      ),
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) 
  
  df_probs <- df_probs_temp %>% 
    fgroup_by(y1, y2, y3, age1, age2, educ1, educ2, timegap1, timegap2, tenure1, tenure2) %>% 
    fsummarise(joint_p = sum(joint_p))
  
  df_lli <- df_estimate %>% 
    join(df_probs, 
         on = c('y1', 
                'y2', 
                'y3', 
                'age1', 
                'age2', 
                'educ1', 
                'educ2', 
                'timegap1', 
                'timegap2', 
                'tenure1',
                'tenure2'), 
         validate = F, 
         verbose = F) %>% 
    fmutate(lli = weight*log(joint_p)) %>%
    pull(lli)
  
  df_lli
  
}

calc_lli_derivatives_3waves_ar1_covariates3 <- function(param_transformed, pi0 = FALSE) {
  
  param <- param_transformed
  if(pi0) param$pi <- 0
  if(!pi0) param$pi <- logit_inverse(param_transformed$pi)
  
  df_probs_temp <- df_template_covariates %>% 
    fmutate(
      theta0_1 = logit_inverse(param$intercept_0 + param$age_0*age1 + param$age2_0*age1^2 + param$educ_0*educ1 + param$timegap*log(timegap1 + 1.5)),
      theta0_2 = logit_inverse(param$intercept_0 + param$age_0*age2 + param$age2_0*age2^2 + param$educ_0*educ2 + param$timegap*log(timegap2 + 1.5)),
      theta1_1 = logit_inverse(param$intercept_1 + param$age_1*age1 + param$age2_1*age1^2 + param$educ_1*educ1 + param$tenure*log(tenure1 + 1.5)),
      theta1_2 = logit_inverse(param$intercept_1 + param$age_1*age2 + param$age2_1*age2^2 + param$educ_1*educ2 + param$tenure*log(tenure2 + 1.5)),
      mu_1 = theta0_1/(theta1_1 + theta0_1)
    ) %>% 
    fmutate(
      p1_star = fifelse(y1_star == 1, mu_1, 1 - mu_1),
      p2_star = fcase(
        y1_star == 0 & y2_star == 0, 1 - theta0_1,
        y1_star == 0 & y2_star == 1, theta0_1,
        y1_star == 1 & y2_star == 0, theta1_1,
        y1_star == 1 & y2_star == 1, 1 - theta1_1
      ),
      p3_star = fcase(
        y2_star == 0 & y3_star == 0, 1 - theta0_2,
        y2_star == 0 & y3_star == 1, theta0_2,
        y2_star == 1 & y3_star == 0, theta1_2,
        y2_star == 1 & y3_star == 1, 1 - theta1_2
      ),
      p1 = fifelse(y1 == y1_star, 1 - param$pi, param$pi),
      p2 = fifelse(y2 == y2_star, 1 - param$pi, param$pi),
      p3 = fifelse(y3 == y3_star, 1 - param$pi, param$pi),
      joint_p = p1_star*p1*p2_star*p2*p3_star*p3
    ) %>% 
    fmutate(
      d1_star_theta0_1 = fcase(
        y1_star == 0, -theta1_1/((theta1_1 + theta0_1)^2),
        y1_star == 1, theta1_1/((theta1_1 + theta0_1)^2)
      ),
      d1_star_theta1_1 = fcase(
        y1_star == 0, theta0_1/((theta1_1 + theta0_1)^2),
        y1_star == 1, -theta0_1/((theta1_1 + theta0_1)^2)
      ),
      d2_star_theta0_1 = fcase(
        y1_star == 0 & y2_star == 0, -1,
        y1_star == 0 & y2_star == 1, 1,
        y1_star == 1 & y2_star == 0, 0,
        y1_star == 1 & y2_star == 1, 0
      ),
      d2_star_theta1_1 = fcase(
        y1_star == 0 & y2_star == 0, 0,
        y1_star == 0 & y2_star == 1, 0,
        y1_star == 1 & y2_star == 0, 1,
        y1_star == 1 & y2_star == 1, -1
      ),
      d3_star_theta0_2 = fcase(
        y2_star == 0 & y3_star == 0, -1,
        y2_star == 0 & y3_star == 1, 1,
        y2_star == 1 & y3_star == 0, 0,
        y2_star == 1 & y3_star == 1, 0
      ),
      d3_star_theta1_2 = fcase(
        y2_star == 0 & y3_star == 0, 0,
        y2_star == 0 & y3_star == 1, 0,
        y2_star == 1 & y3_star == 0, 1,
        y2_star == 1 & y3_star == 1, -1
      ),
      d1_pi = fifelse(y1 == y1_star, -1, 1),
      d2_pi = fifelse(y2 == y2_star, -1, 1),
      d3_pi = fifelse(y3 == y3_star, -1, 1),
      joint_d_theta0_1 = 
        d1_star_theta0_1*joint_p/p1_star + 
        d2_star_theta0_1*joint_p/p2_star,
      joint_d_theta0_2 = 
        d3_star_theta0_2*joint_p/p3_star,
      joint_d_theta1_1 = 
        d1_star_theta1_1*joint_p/p1_star + 
        d2_star_theta1_1*joint_p/p2_star,
      joint_d_theta1_2 = 
        d3_star_theta1_2*joint_p/p3_star,
      joint_d_pi = 
        d1_pi*joint_p/p1 + 
        d2_pi*joint_p/p2 + 
        d3_pi*joint_p/p3
    ) 
  
  df_grad <- df_probs_temp %>% 
    fgroup_by(y1, y2, y3, age1, age2, educ1, educ2, timegap1, timegap2, tenure1, tenure2) %>% 
    fsummarise(
      joint_d_theta0_1 = sum(joint_d_theta0_1),
      joint_d_theta0_2 = sum(joint_d_theta0_2),
      joint_d_theta1_1 = sum(joint_d_theta1_1),
      joint_d_theta1_2 = sum(joint_d_theta1_2), 
      joint_d_pi = sum(joint_d_pi),
      joint_p = sum(joint_p),
      theta0_1 = mean(theta0_1),
      theta0_2 = mean(theta0_2),
      theta1_1 = mean(theta1_1),
      theta1_2 = mean(theta1_2)) %>%
    fungroup() |> 
    fmutate(
      joint_d_theta0_1 = joint_d_theta0_1*theta0_1*(1 - theta0_1), 
      joint_d_theta0_2 = joint_d_theta0_2*theta0_2*(1 - theta0_2),
      joint_d_theta1_1 = joint_d_theta1_1*theta1_1*(1 - theta1_1), 
      joint_d_theta1_2 = joint_d_theta1_2*theta1_2*(1 - theta1_2),
      joint_d_pi = joint_d_pi*param$pi*(1 - param$pi)
    ) %>%
    fmutate(
      joint_d_intercept_0 = joint_d_theta0_1 + joint_d_theta0_2,
      joint_d_age_0 = joint_d_theta0_1*age1 + joint_d_theta0_2*age2,
      joint_d_age2_0 = joint_d_theta0_1*(age1^2) + joint_d_theta0_2*(age2^2),
      joint_d_educ_0 = joint_d_theta0_1*educ1 + joint_d_theta0_2*educ2,
      joint_d_timegap_0 = joint_d_theta0_1*log(timegap1 + 1.5) + joint_d_theta0_2*log(timegap2 + 1.5),
      joint_d_intercept_1 = joint_d_theta1_1 + joint_d_theta1_2,
      joint_d_age_1 = joint_d_theta1_1*age1 + joint_d_theta1_2*age2,
      joint_d_age2_1 = joint_d_theta1_1*(age1^2) + joint_d_theta1_2*(age2^2),
      joint_d_educ_1 = joint_d_theta1_1*educ1 + joint_d_theta1_2*educ2,
      joint_d_tenure_1 = joint_d_theta1_1*log(tenure1 + 1.5) + joint_d_theta1_2*log(tenure2 + 1.5)
    )
  
  df_gi <- df_estimate %>% 
    join(df_grad, 
         on       = c('y1', 
                      'y2', 
                      'y3', 
                      'age1', 
                      'age2', 
                      'educ1', 
                      'educ2', 
                      'timegap1', 
                      'timegap2', 
                      'tenure1', 
                      'tenure2'), 
         validate = F, 
         verbose  = F) %>% 
    fmutate(
      lgi_intercept_0 = weight*joint_d_intercept_0/joint_p,
      lgi_age_0       = weight*joint_d_age_0/joint_p,
      lgi_age2_0      = weight*joint_d_age2_0/joint_p,
      lgi_educ_0      = weight*joint_d_educ_0/joint_p,
      lgi_timegap_0   = weight*joint_d_timegap_0/joint_p,
      lgi_intercept_1 = weight*joint_d_intercept_1/joint_p,
      lgi_age_1       = weight*joint_d_age_1/joint_p,
      lgi_age2_1      = weight*joint_d_age2_1/joint_p,
      lgi_educ_1      = weight*joint_d_educ_1/joint_p,
      lgi_tenure_1    = weight*joint_d_tenure_1/joint_p,
      lgi_pi          = weight*joint_d_pi/joint_p) %>%
    fselect(lgi_intercept_0, 
            lgi_age_0, 
            lgi_age2_0, 
            lgi_educ_0, 
            lgi_timegap_0, 
            lgi_intercept_1, 
            lgi_age_1, 
            lgi_age2_1, 
            lgi_educ_1, 
            lgi_tenure_1, 
            lgi_pi)
  df_gi
  
}

calc_mle_3waves_ar1_covariates3 <- function(param_transformed) {
  sum(calc_lli_3waves_ar3_covariates3(param_transformed))
}

calc_mle_derivatives_3waves_ar1_covariates3 <- function(param_transformed) {
  colSums(calc_lli_derivatives_3waves_ar1_covariates3(param_transformed))
}
