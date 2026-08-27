test_that("constant reliability design reproduces constant symmetric error", {
  set.seed(91)
  n <- 80L
  df <- data.frame(y1 = rbinom(n, 1, .5), y2 = rbinom(n, 1, .5),
                   y3 = rbinom(n, 1, .5), weight = runif(n, .5, 2))
  X <- matrix(1, n, 1, dimnames = list(NULL, "intercept"))
  attr(X, "entry_active") <- attr(X, "persistence_active") <- TRUE
  Z <- replicate(3L, matrix(1, n, 1,
    dimnames = list(NULL, "error_intercept")), simplify = FALSE)
  params <- list(beta0 = c(intercept = qnorm(.08)),
                 beta1 = c(intercept = qnorm(.93)), alpha = .48,
                 delta = c(error_intercept = qlogis(2 * .03)))
  eta <- pack_covariate_reliability(params, X)
  data <- list(df = df, X = X, Z = Z)
  combined <- .covrel_components(eta, data)
  constant <- e_step_covariates(df, X,
    list(beta0 = params$beta0, beta1 = params$beta1,
         alpha = params$alpha, pi = .03),
    model_type = "symmetric", stationary = FALSE)
  expect_equal(sum(df$weight * combined$loglik_i), constant$loglik,
               tolerance = 1e-8)
})

test_that("combined-model analytic score matches numerical gradient", {
  set.seed(92)
  n <- 55L
  df <- data.frame(y1 = rbinom(n, 1, .5), y2 = rbinom(n, 1, .5),
                   y3 = rbinom(n, 1, .5), weight = runif(n, .5, 2))
  X <- cbind(intercept = 1, x = rnorm(n))
  attr(X, "entry_active") <- attr(X, "persistence_active") <- c(TRUE, TRUE)
  Z <- lapply(seq_len(3L), function(tt)
    cbind(error_intercept = 1, unreliable = rbinom(n, 1, .2)))
  params <- list(beta0 = c(intercept = -1.2, x = .1),
    beta1 = c(intercept = 1.3, x = -.1), alpha = .47,
    delta = c(error_intercept = -2.5, unreliable = .4))
  eta <- pack_covariate_reliability(params, X)
  data <- list(df = df, X = X, Z = Z)
  analytic <- colSums(.covrel_components(eta, data, TRUE)$scores * df$weight)
  fn <- function(z) sum(df$weight * .covrel_components(z, data)$loglik_i)
  numeric <- numeric(length(eta))
  for (j in seq_along(eta)) {
    hh <- 1e-6 * max(1, abs(eta[j])); xp <- xm <- eta
    xp[j] <- xp[j] + hh; xm[j] <- xm[j] - hh
    numeric[j] <- (fn(xp) - fn(xm)) / (2 * hh)
  }
  expect_equal(unname(analytic), unname(numeric), tolerance = 2e-5)
})

test_that("conditional-stationarity score matches numerical gradient", {
  set.seed(93)
  n <- 60L
  df <- data.frame(y1 = rbinom(n, 1, .5), y2 = rbinom(n, 1, .5),
                   y3 = rbinom(n, 1, .5), weight = runif(n, .5, 2))
  X12 <- cbind(intercept = 1, x = rnorm(n))
  X23 <- cbind(intercept = 1, x = rnorm(n))
  for (X in list(X12, X23)) {
    attr(X, "entry_active") <- attr(X, "persistence_active") <- c(TRUE, TRUE)
  }
  attr(X12, "entry_active") <- attr(X12, "persistence_active") <- c(TRUE, TRUE)
  attr(X23, "entry_active") <- attr(X23, "persistence_active") <- c(TRUE, TRUE)
  X <- list(X12 = X12, X23 = X23)
  Z <- lapply(seq_len(3L), function(tt)
    cbind(error_intercept = 1, unreliable = rbinom(n, 1, .2)))
  params <- list(beta0 = c(intercept = -1.1, x = .15),
    beta1 = c(intercept = 1.2, x = -.12),
    delta = c(error_intercept = -2.6, unreliable = .3))
  eta <- pack_covariate_reliability(params, X, stationary = TRUE)
  data <- list(df = df, X = X, Z = Z)
  detail <- .covrel_components(
    eta, data, retain_scores = TRUE, stationary = TRUE)
  analytic <- colSums(detail$scores * df$weight)
  fn <- function(z) sum(df$weight * .covrel_components(
    z, data, stationary = TRUE)$loglik_i)
  numeric <- numeric(length(eta))
  for (j in seq_along(eta)) {
    hh <- 1e-6 * max(1, abs(eta[j])); xp <- xm <- eta
    xp[j] <- xp[j] + hh; xm[j] <- xm[j] - hh
    numeric[j] <- (fn(xp) - fn(xm)) / (2 * hh)
  }
  expect_length(eta, 6L)
  expect_false("initial_employment_logit" %in% names(eta))
  expect_equal(unname(analytic), unname(numeric), tolerance = 2e-5)
  expect_true(all(detail$alpha_i > 0 & detail$alpha_i < 1))
})
