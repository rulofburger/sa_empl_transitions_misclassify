test_that("four-wave histories enumerate every binary path", {
  histories <- latent_histories_eps_4w()
  expect_equal(dim(histories), c(16L, 4L))
  expect_equal(nrow(unique(histories)), 16L)
  expect_true(all(histories %in% 0:1))
})

test_that("four-wave preparation preserves zero tenure and recodes never-worked", {
  df <- data.frame(weight = c(1, 2, 3))
  for (wave in 1:4) {
    df[[paste0("y", wave)]] <- c(1L, 0L, 0L)
    df[[paste0("tenure", wave)]] <- c(0, -1, -1)
    df[[paste0("timegap", wave)]] <- c(NA, 1, 2)
    df[[paste0("age", wave)]] <- c(25, 30, 30)
    df[[paste0("neverworked", wave)]] <- c(0L, 1L, NA)
  }
  prepared <- prepare_eps_estimation_data_4w(df)
  expect_equal(nrow(prepared), 3L)
  expect_equal(prepared$tenure1, c(0, NA, NA))
  expect_equal(prepared$timegap_cat1, c(NA_integer_, 7L, 2L))
  df$tenure4[1L] <- -1
  expect_equal(nrow(prepare_eps_estimation_data_4w(df)), 2L)
})

test_that("four-wave joint timegap distribution sums to one", {
  categories <- as.matrix(expand.grid(rep(list(1:7), 4L)))
  for (eps_d in c(1e-8, .12, .7)) {
    mass <- exp(log_emission_timegap_spell_joint(categories,
      c(.45, .21, .09, .10, .06), eps_d = eps_d))
    expect_equal(sum(mass), 1, tolerance = 1e-10)
  }
})

test_that("job-segment compression preserves all likelihood outputs", {
  g <- rbind(c(12, 15, 18, 21), c(NA, 1, 4, 7), c(4, 9, NA, 12)) / 12
  g <- g[rep(1:3, each = 3), ]
  arguments <- list(g_mat = g, s_mat = is.finite(g), t_offsets = 0:3,
    lambda_g = c(.25, .28, .17, .15, .13), eps = .2,
    job_change_prob = .03, tenure_heaping_prob = .02,
    tenure_year_revision_prob = .1, tenure_clean_anchor_revision_prob = .2,
    tenure_exact_anchor_retention_prob = .25,
    tenure_local_revision_prob = .3,
    tenure_start_month_probs = (1:12) / sum(1:12),
    interview_month_mat = matrix(rep(c(3, 6, 9, 12), each = nrow(g)),
      nrow(g), 4))
  compressed <- do.call(log_emission_spell_g_monthly, arguments)
  arguments$compress_segments <- FALSE
  reference <- do.call(log_emission_spell_g_monthly, arguments)
  expect_equal(compressed, reference, tolerance = 1e-12)
})

test_that("generic joint timegap spell reproduces the triplet likelihood", {
  categories <- matrix(c(1, 1, 2, 2, 2, 3, 4, 5, 6),
    nrow = 3L, byrow = TRUE)
  hazards <- c(.45, .21, .09, .10, .06)
  generic <- log_emission_timegap_spell_joint(categories, hazards,
    eps_d = .12)
  established <- log_emission_timegap_triplet_joint(
    categories[, 1L], categories[, 2L], categories[, 3L], hazards,
    eps_d = .12)
  expect_equal(generic, established, tolerance = 1e-12)
})

test_that("four-wave posterior histories normalize", {
  df <- data.frame(
    y1 = c(1L, 0L, 1L), y2 = c(1L, 0L, 0L),
    y3 = c(1L, 0L, 0L), y4 = c(1L, 1L, 0L),
    tenure1 = c(2, NA, 1), tenure2 = c(2.25, NA, NA),
    tenure3 = c(2.5, NA, NA), tenure4 = c(2.75, 0, NA),
    timegap_cat1 = c(NA, 2L, NA), timegap_cat2 = c(NA, 3L, 1L),
    timegap_cat3 = c(NA, 3L, 2L), timegap_cat4 = c(NA, NA, 3L),
    interview_month1 = 3L, interview_month2 = 6L,
    interview_month3 = 9L, interview_month4 = 12L,
    weight = 1)
  params <- list(alpha = .55, theta0 = .03, theta1 = .95,
    pi = .03, eps = .15, eps_d = .08, job_change_prob = .02,
    lambda_g = c(.25, .28, .17, .15, .13),
    lambda_d = c(.45, .21, .09, .10, .06), beta_g = 0, beta_d = 0,
    tenure_measurement_model = "monthly",
    timegap_contamination_model = "joint_marginal",
    tenure_heaping_prob = 0, tenure_year_revision_prob = 0,
    tenure_clean_anchor_revision_prob = 0,
    tenure_exact_anchor_retention_prob = 0,
    tenure_local_revision_prob = 0,
    tenure_start_month_probs = rep(1/12, 12L))
  result <- e_step_eps_4w(df, params)
  expect_equal(dim(result$gamma), c(3L, 16L))
  expect_equal(rowSums(result$gamma), rep(1, 3L), tolerance = 1e-12)
  expect_true(is.finite(result$loglik))
})

test_that("four-wave exact-status limit selects the reported history", {
  df <- data.frame(
    y1 = 1L, y2 = 0L, y3 = 0L, y4 = 1L,
    tenure1 = 1, tenure2 = NA, tenure3 = NA, tenure4 = 0,
    timegap_cat1 = NA, timegap_cat2 = 1L, timegap_cat3 = 2L,
    timegap_cat4 = NA,
    interview_month1 = 3L, interview_month2 = 6L,
    interview_month3 = 9L, interview_month4 = 12L,
    weight = 1)
  params <- list(alpha = .5, pi = 0, eps = .2, eps_d = .1,
    job_change_prob = 0, lambda_g = rep(.2, 5L),
    lambda_d = rep(.3, 5L), beta_g = 0, beta_d = 0,
    tenure_measurement_model = "monthly",
    timegap_contamination_model = "joint_marginal",
    tenure_heaping_prob = 0, tenure_year_revision_prob = 0,
    tenure_clean_anchor_revision_prob = 0,
    tenure_exact_anchor_retention_prob = 0,
    tenure_local_revision_prob = 0,
    tenure_start_month_probs = rep(1/12, 12L))
  result <- e_step_eps_4w(df, params)
  expected <- 1L + 1L + 2L * 0L + 4L * 0L + 8L * 1L
  expect_equal(which(result$gamma[1L, ] == 1), expected)
})
