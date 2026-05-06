# Extracted from test-estep.R:385

# prequel ----------------------------------------------------------------------
.make_df <- function(n = 50, theta0 = 0.1, theta1 = 0.9, pi = 0.05, seed = 42) {
  set.seed(seed)
  alpha <- theta0 / (theta0 + 1 - theta1)

  h1 <- rbinom(n, 1, alpha)
  h2 <- ifelse(h1 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))
  h3 <- ifelse(h2 == 1, rbinom(n, 1, theta1), rbinom(n, 1, theta0))

  flip <- function(h) ifelse(rbinom(length(h), 1, pi) == 1, 1 - h, h)
  data.frame(y1 = flip(h1), y2 = flip(h2), y3 = flip(h3), weight = rep(1, n))
}
.make_df_perfect <- function(n = 60, theta0 = 0.1, theta1 = 0.9) {
  alpha <- theta0 / (theta0 + 1 - theta1)
  # Alternate employed/nonemployed for deterministic mix
  n_emp  <- round(n * alpha)
  n_nemp <- n - n_emp
  # All employment spells stay employed (theta1=1 approx); nonemployment stays
  h1 <- c(rep(1L, n_emp), rep(0L, n_nemp))
  h2 <- h1  # deterministic: no transitions
  h3 <- h1
  data.frame(y1 = h1, y2 = h2, y3 = h3, weight = rep(1, n))
}

# test -------------------------------------------------------------------------
df     <- .make_df(n = 70)
params <- list(alpha = 0.5, theta0 = 0.1, theta1 = 0.9, pi = 0.05)
result <- e_step(df, params, model_type = "symmetric")
W <- sum(df$weight)
expect_equal(result$suff$H0 + result$suff$H1, 3 * W, tolerance = 1e-8)
