make_ar2r_test_data <- function(n = 18L, zcols = 8L, seed = 91L) {
  set.seed(seed)
  xnames <- c("age", "age_sq", "educ", "female", "log_tenure",
              "log_time_since_work", "never_worked", "tenure_missing",
              "permanent_contract", "informal_sector")
  X <- lapply(1:2, function(tt) {
    x <- matrix(rnorm(n * length(xnames)), n, dimnames = list(NULL, xnames))
    x[, c("never_worked", "tenure_missing", "permanent_contract",
          "informal_sector")] <- rbinom(n * 4, 1, .25)
    x
  })
  znames <- c("error_intercept", "age_inconsistency",
    "education_inconsistency", "race_inconsistency", "gender_inconsistency",
    "two_inconsistencies", "three_inconsistencies", "four_inconsistencies")
  Z <- lapply(1:4, function(tt) cbind(error_intercept = 1,
    matrix(rbinom(n * 7, 1, .12), n, dimnames = list(NULL, znames[-1]))))
  list(y = matrix(rbinom(n * 4, 1, .5), n, 4), weight = runif(n, .5, 1.5),
    weight_sq = rep(1, n), X = X, Z = lapply(Z, function(z) z[,1:zcols,drop=FALSE]),
    entry_active = !xnames %in% c("log_tenure", "tenure_missing",
      "permanent_contract", "informal_sector"),
    persistence_active = !xnames %in% c("log_time_since_work", "never_worked"),
    covariate_names = xnames, n_original = n, n_cells = n, stage = "test")
}

test_that("four-wave inconsistency attribution extends the Table 3 rule", {
  x <- data.frame(age1=30,age2=35,age3=36,age4=42,
                  educ1=10,educ2=10,educ3=7,educ4=8,
                  race1=1,race2=1,race3=2,race4=2,
                  female1=0,female2=1,female3=0,female4=0)
  z <- compute_inconsistency_extent_4w(x)
  expect_equal(as.integer(z[1,paste0("Y_age_",1:4)]), c(1,0,0,1))
  expect_equal(as.integer(z[1,paste0("Y_edu_",1:4)]), c(0,0,1,0))
  expect_equal(as.integer(z[1,paste0("Y_race_",1:4)]), c(0,0,1,0))
  expect_equal(as.integer(z[1,paste0("Y_gender_",1:4)]), c(0,1,0,0))
  expect_equal(as.integer(z[1,paste0("Y_exactly_2_",1:4)]), c(0,0,1,0))
  expect_true(z$extent_age_1 >= 2 && z$extent_age_4 >= 2)
  expect_true(z$extent_edu_3 >= 2)
})

test_that("piecewise duration bins use the stated intervals", {
  x <- c(0, 3, 3.1, 6, 6.1, 12, 12.1, 24, 24.1, 60, 60.1, NA, -1)
  bins <- .ar2_piecewise_duration_bins(x, "tenure")
  expect_equal(colnames(bins), c("tenure_4_6m", "tenure_7_12m",
    "tenure_13_24m", "tenure_25_60m", "tenure_over_60m"))
  expect_equal(rowSums(bins), c(0, 0, rep(1, 9), 0, 0))
  expect_equal(which(bins[, "tenure_4_6m"] == 1), 3:4)
  expect_equal(which(bins[, "tenure_7_12m"] == 1), 5:6)
  expect_equal(which(bins[, "tenure_13_24m"] == 1), 7:8)
  expect_equal(which(bins[, "tenure_25_60m"] == 1), 9:10)
  expect_equal(which(bins[, "tenure_over_60m"] == 1), 11)
})

test_that("AR2 Set 2 analytical score matches numerical derivatives", {
  data <- make_ar2r_test_data()
  eta <- setNames(rnorm(length(ar2_reliability_names(data)), 0, .18),
                  ar2_reliability_names(data))
  eta[grep("entry_intercept", names(eta))] <- -1.5
  eta[grep("persistence_intercept", names(eta))] <- 1.4
  eta["error_intercept"] <- -3
  analytic <- colSums(.ar2r_components(eta, data, TRUE)$scores * data$weight)
  fn <- function(z) sum(data$weight * .ar2r_components(z, data)$loglik_i)
  numeric <- vapply(seq_along(eta), function(j) {
    hh <- 1e-5; zp <- zm <- eta; zp[j] <- zp[j] + hh; zm[j] <- zm[j] - hh
    (fn(zp) - fn(zm)) / (2 * hh)
  }, numeric(1))
  expect_equal(unname(analytic), unname(numeric), tolerance = 4e-4)
})

test_that("reliability stages are exact nested restrictions", {
  full <- make_ar2r_test_data()
  full <- subset_ar2_reliability_stage(full, "table3_column1")
  constant <- subset_ar2_reliability_stage(full, "constant")
  expect_equal(colnames(full$Z[[1]]), c("error_intercept",
    "age_inconsistency", "education_inconsistency", "race_inconsistency",
    "gender_inconsistency", "two_inconsistencies", "three_inconsistencies",
    "four_inconsistencies"))
  e0 <- setNames(rnorm(length(ar2_reliability_names(constant)), 0, .1),
                 ar2_reliability_names(constant))
  e0["error_intercept"] <- -3
  ef <- setNames(rep(0, length(ar2_reliability_names(full))), ar2_reliability_names(full))
  ef[names(e0)] <- e0
  expect_equal(.ar2r_components(e0, constant)$loglik_i,
               .ar2r_components(ef, full)$loglik_i, tolerance = 1e-11)
  q <- ar2_reliability_quantities(ef, full)
  expect_true(all(is.finite(q)) && all(q > 0) && all(q < 1))
})

test_that("prospective apparent-transition events are exhaustive", {
  data <- make_ar2r_test_data()
  data$y[1:2, 1:3] <- rbind(c(0, 1, 0), c(1, 0, 1))
  eta <- setNames(rep(0, length(ar2_reliability_names(data))),
                  ar2_reliability_names(data))
  eta["entry_intercept"] <- -1.4
  eta["persistence_intercept"] <- 1.4
  eta["error_intercept"] <- -3
  fit <- list(eta = eta, params = .ar2r_unpack(eta, data))
  out <- prospective_apparent_transition_ar2(data, fit)
  expect_equal(out$direction,
    c("All apparent transitions", "Apparent entries", "Apparent exits"))
  expect_true(all(is.finite(as.matrix(out[, -1]))))
  expect_equal(rowSums(out[, -1]), rep(1, 3), tolerance = 1e-12)
})
