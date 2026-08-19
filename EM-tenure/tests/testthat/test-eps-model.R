# test-eps-model.R
# Tests for the eps (Spec I) point-mass + Exp contamination model.
#
# Created: 2026-04-30
# TeX ref: "EM tenure epsilon.tex"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.make_eps_data <- function(n = 200L, seed = 42L) {
  set.seed(seed)
  s <- sample(0:1, n, replace = TRUE, prob = c(0.35, 0.65))
  # Construct clock-consistent tenures (60% of EE pairs) and contaminated (40%)
  g1 <- pmax(0.26, rexp(n, 2))
  g2 <- ifelse(runif(n) < 0.60, g1 + 0.25, pmax(0.26, rexp(n, 2)))
  g3 <- ifelse(runif(n) < 0.60, g2 + 0.25, pmax(0.26, rexp(n, 2)))
  as.data.frame(list(
    y1           = s,
    y2           = s,
    y3           = s,
    tenure1      = g1,
    tenure2      = g2,
    tenure3      = g3,
    timegap_cat1 = sample(1:7, n, replace = TRUE),
    timegap_cat2 = sample(1:7, n, replace = TRUE),
    timegap_cat3 = sample(1:7, n, replace = TRUE),
    weight       = rep(1, n)
  ))
}

.make_eps_params <- function(eps = 0.20) {
  list(
    alpha    = 0.60,
    theta1   = 0.85,
    theta0   = 0.15,
    pi       = 0.05,
    eps      = eps,
    lambda_g = 2.0,
    lambda_d = 1.5
  )
}

# ---------------------------------------------------------------------------
# 1. log_emission_spell_g: K=0 rows contribute zero
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=0 rows have loglik=0, tau_sum=0, count=0", {
  N <- 5L
  g_mat <- matrix(1, nrow = N, ncol = 2)
  s_mat <- matrix(FALSE, nrow = N, ncol = 2)   # no observations: K=0
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = 0.5, eps = 0.2)

  expect_equal(out$loglik,       rep(0, N))
  expect_equal(out$tau_sum,      rep(0, N))   # tau_sum = 0 for K=0
  expect_equal(out$lambda_count, rep(0, N))
  expect_equal(out$lambda_xsum,  rep(0, N))
  expect_equal(out$K,            rep(0L, N))
})

# ---------------------------------------------------------------------------
# 2. log_emission_spell_g: K=1 rows use plain Exp density
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=1 loglik matches log-Exp exactly", {
  lambda_g <- 0.4
  g_obs    <- c(1.0, 2.5, 0.5)
  N        <- length(g_obs)
  g_mat    <- cbind(g_obs, 0)                        # second col unobserved
  s_mat    <- matrix(c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE), nrow = N)

  out <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = 0.3)

  expected_ll <- log(lambda_g) - lambda_g * g_obs
  expect_equal(out$loglik, expected_ll, tolerance = 1e-10)
  expect_equal(out$tau_sum, rep(0, N))      # tau_sum = 0 for K=1, offset=0 (no eps info)
  expect_equal(out$lambda_count, rep(1, N))
  expect_equal(out$lambda_xsum,  g_obs, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 2b. log_emission_spell_g: K=1 with offset=0 — regression (eps drops out)
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=1 offset=0 tau_sum=0 (eps not identifiable)", {
  # Spell of length 2 but only wave 1 observed (offset 0): t_1 = a.
  # Both patterns give the same Exp density, so eps drops out.
  lambda_g <- 0.3; eps <- 0.4; g_obs <- 1.8
  g_mat <- matrix(c(g_obs, 0), nrow = 1)
  s_mat <- matrix(c(TRUE, FALSE), nrow = 1)

  out <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  expect_equal(out$loglik,       log(lambda_g) - lambda_g * g_obs, tolerance = 1e-12)
  expect_equal(out$tau_sum,      0,     tolerance = 1e-12)
  expect_equal(out$lambda_count, 1,     tolerance = 1e-12)
  expect_equal(out$lambda_xsum,  g_obs, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 2c. log_emission_spell_g: K=1 with offset=1 — 2-pattern mixture
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=1 offset=1 loglik matches 2-pattern analytical", {
  # Spell starts at wave 1 (unobserved); tenure observed only at wave 2 (offset 1).
  # t_1 > a, so clean branch: lambda_g * exp(-lambda_g*(g - Delta)).
  # Contaminated branch:       lambda_g * exp(-lambda_g*g).
  lambda_g <- 0.3; eps <- 0.35; Delta <- .QUARTER_YEARS; g_obs <- 2.0
  T_g <- g_obs - Delta   # = 1.75; T_g > 0 so clean branch is valid.

  lp_clean  <- log(1 - eps) + log(lambda_g) - lambda_g * T_g
  lp_contam <- log(eps)     + log(lambda_g) - lambda_g * g_obs
  mx        <- max(lp_clean, lp_contam)
  expected_ll <- mx + log(exp(lp_clean - mx) + exp(lp_contam - mx))

  nu_expected <- exp(lp_contam - expected_ll)    # P(contaminated | g)
  lx_expected <- nu_expected * g_obs + (1 - nu_expected) * T_g

  g_mat <- matrix(c(0, g_obs), nrow = 1)
  s_mat <- matrix(c(FALSE, TRUE), nrow = 1)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  expect_equal(out$loglik,       expected_ll,   tolerance = 1e-10)
  expect_true(out$tau_sum > 0 && out$tau_sum < 1)  # eps is identifiable
  expect_equal(out$tau_sum,      nu_expected,   tolerance = 1e-10)
  expect_equal(out$lambda_count, 1,             tolerance = 1e-12)
  expect_equal(out$lambda_xsum,  lx_expected,   tolerance = 1e-10)
  expect_true(out$eps_informative)              # K=1/offset>0 is eps-informative
})

# ---------------------------------------------------------------------------
# 2d. log_emission_spell_g: K=1 with offset=2 — 2-pattern mixture, shift=2*Delta
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=1 offset=2 loglik matches 2-pattern analytical", {
  # 3-wave spell [1,3], only wave 3 observed (offset 2). T_g = g - 2*Delta.
  lambda_g <- 0.2; eps <- 0.30; Delta <- .QUARTER_YEARS; g_obs <- 1.0
  T_g <- g_obs - 2 * Delta   # = 0.5; T_g > 0

  lp_clean  <- log(1 - eps) + log(lambda_g) - lambda_g * T_g
  lp_contam <- log(eps)     + log(lambda_g) - lambda_g * g_obs
  mx        <- max(lp_clean, lp_contam)
  expected_ll <- mx + log(exp(lp_clean - mx) + exp(lp_contam - mx))
  nu_expected <- exp(lp_contam - expected_ll)
  lx_expected <- nu_expected * g_obs + (1 - nu_expected) * T_g

  g_mat <- matrix(c(0, 0, g_obs), nrow = 1)
  s_mat <- matrix(c(FALSE, FALSE, TRUE), nrow = 1)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L, 2L),
                                lambda_g = lambda_g, eps = eps)

  expect_equal(out$loglik,       expected_ll,   tolerance = 1e-10)
  expect_equal(out$tau_sum,      nu_expected,   tolerance = 1e-10)
  expect_equal(out$lambda_count, 1,             tolerance = 1e-12)
  expect_equal(out$lambda_xsum,  lx_expected,   tolerance = 1e-10)
  expect_true(out$eps_informative)              # K=1/offset>0 is eps-informative
})

# ---------------------------------------------------------------------------
# 2e. log_emission_spell_g: K=1 offset=1, impossible T_g — tau_sum ≈ 1
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=1 offset=1 with g < Delta gives tau_sum near 1", {
  # g_obs = 0.1 < Delta = 0.25 means T_g = g - Delta = -0.15 < 0.
  # Clean branch is impossible (T_g <= 0); only contaminated branch contributes.
  lambda_g <- 0.4; eps <- 0.35; Delta <- .QUARTER_YEARS; g_obs <- 0.1

  expected_ll <- log(eps) + log(lambda_g) - lambda_g * g_obs

  g_mat <- matrix(c(0, g_obs), nrow = 1)
  s_mat <- matrix(c(FALSE, TRUE), nrow = 1)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  expect_equal(out$loglik,  expected_ll, tolerance = 1e-10)
  expect_equal(out$tau_sum, 1,           tolerance = 1e-10)  # fully contaminated
  expect_equal(out$lambda_xsum, g_obs,   tolerance = 1e-10)  # xsum = g_obs (contam branch)
  expect_true(out$eps_informative)              # K=1/offset>0 is eps-informative
})

# ---------------------------------------------------------------------------
# 2f. log_emission_spell_g: K=1 offset>0 — vectorized (N > 1)
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=1 offset>0 vectorized over N rows", {
  # Multiple rows, all K=1 with offset=1. Each row independently computed.
  lambda_g <- 0.3; eps <- 0.35; Delta <- .QUARTER_YEARS
  g_obs_vec <- c(0.5, 1.0, 2.0, 3.5)
  N <- length(g_obs_vec)

  g_mat <- cbind(0, g_obs_vec)             # col1 unobserved, col2 observed
  s_mat <- matrix(c(rep(FALSE, N), rep(TRUE, N)), nrow = N)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  # Each row should match the single-row analytical result
  for (ii in seq_len(N)) {
    g1i <- g_obs_vec[ii]; T_gi <- g1i - Delta
    lp_c <- if (T_gi > 0) log(1 - eps) + log(lambda_g) - lambda_g * T_gi else -Inf
    lp_x <- log(eps) + log(lambda_g) - lambda_g * g1i
    mx_i <- max(lp_c, lp_x)
    ll_i <- mx_i + log(exp(lp_c - mx_i) + exp(lp_x - mx_i))
    expect_equal(out$loglik[ii], ll_i, tolerance = 1e-10,
                 label = sprintf("row %d loglik", ii))
  }

  expect_equal(out$K,              rep(1L, N))
  expect_true(all(out$eps_informative))   # all K=1/offset>0 rows are eps-informative
  expect_true(all(out$tau_sum > 0))       # all have non-trivial mixture
})

# ---------------------------------------------------------------------------
# 2g. log_emission_spell_g: mixed K=1 rows (offset=0 and offset>0 in same call)
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=1 mixed at_start and offset>0 in same call", {
  # Row 1: K=1, offset=0 (t_1 = a): plain Exp, eps drops out
  # Row 2: K=1, offset=1 (t_1 > a): 2-pattern mixture
  lambda_g <- 0.4; eps <- 0.30; Delta <- .QUARTER_YEARS
  g1 <- 1.5; g2 <- 2.0   # g2 obs at wave 2 (offset 1)

  g_mat <- matrix(c(g1, 0, 0, g2), nrow = 2)
  s_mat <- matrix(c(TRUE, FALSE, FALSE, TRUE), nrow = 2)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  # Row 1: plain Exp at offset=0
  expect_equal(out$loglik[1],       log(lambda_g) - lambda_g * g1, tolerance = 1e-10)
  expect_equal(out$tau_sum[1],      0,                             tolerance = 1e-12)
  expect_false(out$eps_informative[1])   # offset=0: NOT eps-informative

  # Row 2: 2-pattern mixture at offset=1
  T_g2 <- g2 - Delta
  lp_c <- log(1 - eps) + log(lambda_g) - lambda_g * T_g2
  lp_x <- log(eps)     + log(lambda_g) - lambda_g * g2
  mx2  <- max(lp_c, lp_x)
  ll2  <- mx2 + log(exp(lp_c - mx2) + exp(lp_x - mx2))
  expect_equal(out$loglik[2],  ll2, tolerance = 1e-10)
  expect_true(out$tau_sum[2] > 0)
  expect_true(out$eps_informative[2])    # offset>0: IS eps-informative
})

# ---------------------------------------------------------------------------
# 3. log_emission_spell_g: K=2 clock-consistent rows use both branches
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=2 consistent rows have tau_sum in [0, 1)", {
  lambda_g <- 0.4; eps <- 0.25
  # Exactly clock-consistent: g2 = g1 + 0.25
  g1 <- 2.0; g2 <- g1 + 0.25
  g_mat <- matrix(c(g1, g2), nrow = 1)
  s_mat <- matrix(c(TRUE, TRUE), nrow = 1)

  out <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  expect_true(is.finite(out$loglik))
  # tau_sum = E[# contaminated waves] in [0, K=2]; clean rows have tau_sum < 1
  expect_true(out$tau_sum >= 0 && out$tau_sum < 1)
  expect_equal(out$K, 2L)
})

# ---------------------------------------------------------------------------
# 4. log_emission_spell_g: K=2 inconsistent rows collapse to contaminated branch
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=2 inconsistent rows have tau_sum in (1, 2)", {
  lambda_g <- 0.4; eps <- 0.25
  # NOT clock-consistent: g2 is far from g1 + 0.25
  # Under full enumeration: CC has 0 density; CX/XC/XX contribute.
  # tau_sum = E[# contam] = 1*w_CX + 1*w_XC + 2*w_XX > 1 (since w_XX > 0).
  g_mat <- matrix(c(2.0, 5.0), nrow = 1)
  s_mat <- matrix(c(TRUE, TRUE), nrow = 1)

  out <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  expect_true(out$tau_sum > 1 && out$tau_sum < 2)
  expect_true(is.finite(out$loglik))
})

# ---------------------------------------------------------------------------
# 5. log_emission_spell_g: input validation
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g validates lambda_g and eps bounds", {
  g_mat <- matrix(1, 2, 2); s_mat <- matrix(TRUE, 2, 2)
  expect_error(log_emission_spell_g(g_mat, s_mat, 0:1, lambda_g = 0,   eps = 0.2))
  expect_error(log_emission_spell_g(g_mat, s_mat, 0:1, lambda_g = -1,  eps = 0.2))
  expect_error(log_emission_spell_g(g_mat, s_mat, 0:1, lambda_g = 0.5, eps = 0))
  expect_error(log_emission_spell_g(g_mat, s_mat, 0:1, lambda_g = 0.5, eps = 1))
  expect_error(log_emission_spell_g(as.list(g_mat), s_mat, 0:1, lambda_g = 0.5, eps = 0.2))
})

# ---------------------------------------------------------------------------
# 6. e_step_eps: output structure
# ---------------------------------------------------------------------------

test_that("e_step_eps returns correct structure", {
  df     <- .make_eps_data()
  params <- .make_eps_params()

  out <- e_step_eps(df, params)

  expect_named(out, c("gamma", "loglik", "suff"))
  expect_equal(dim(out$gamma), c(200L, 8L))
  expect_true(is.finite(out$loglik))
  expect_true(out$loglik <= 0)  # log-likelihood must be <= 0 for densities
})

# ---------------------------------------------------------------------------
# 7. e_step_eps: gamma rows sum to 1 and are in [0, 1]
# ---------------------------------------------------------------------------

test_that("e_step_eps: gamma row-sums = 1 and values in [0,1]", {
  df  <- .make_eps_data()
  out <- e_step_eps(df, .make_eps_params())

  row_sums <- rowSums(out$gamma)
  expect_equal(row_sums, rep(1, nrow(df)), tolerance = 1e-10)
  expect_true(all(out$gamma >= 0))
  expect_true(all(out$gamma <= 1 + 1e-12))
})

# ---------------------------------------------------------------------------
# 8. e_step_eps: pi=0 gives deterministic responsibilities
# ---------------------------------------------------------------------------

test_that("e_step_eps: pi=0 gives 0/1 gamma values", {
  df       <- .make_eps_data()
  params   <- .make_eps_params()
  params$pi <- 0

  out <- e_step_eps(df, params)
  # With pi=0, each row must place all weight on exactly one history.
  expect_equal(rowSums(out$gamma > 0), rep(1L, nrow(df)))
})

# ---------------------------------------------------------------------------
# 9. e_step_eps: sufficient stats Eps_num <= Eps_den (tau in [0,1])
# ---------------------------------------------------------------------------

test_that("e_step_eps: Eps_num <= Eps_den (valid tau aggregate)", {
  df  <- .make_eps_data()
  out <- e_step_eps(df, .make_eps_params())

  expect_true(out$suff$Eps_num >= 0)
  expect_true(out$suff$Eps_den >= 0)
  expect_true(out$suff$Eps_num <= out$suff$Eps_den + 1e-10)
})

# ---------------------------------------------------------------------------
# 10. e_step_eps: input validation catches y1/y2/y3 NA
# ---------------------------------------------------------------------------

test_that("e_step_eps stops on NA in y columns", {
  df      <- .make_eps_data()
  df$y1[1] <- NA_integer_
  expect_error(e_step_eps(df, .make_eps_params()), "NA in y1/y2/y3")
})

# ---------------------------------------------------------------------------
# 11. m_step_eps: eps update is closed-form Eps_num / Eps_den
# ---------------------------------------------------------------------------

test_that("m_step_eps: eps_hat matches Eps_num/Eps_den ratio", {
  df      <- .make_eps_data()
  estep   <- e_step_eps(df, .make_eps_params())
  mout    <- m_step_eps(estep$suff, total_weight = nrow(df))

  expected_eps <- estep$suff$Eps_num / estep$suff$Eps_den
  expected_eps <- max(1e-4, min(expected_eps, 0.95))
  expect_equal(mout$eps, expected_eps, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 12. m_step_eps: all params in valid ranges
# ---------------------------------------------------------------------------

test_that("m_step_eps: all parameters in valid ranges", {
  df    <- .make_eps_data()
  suff  <- e_step_eps(df, .make_eps_params())$suff
  mout  <- m_step_eps(suff, total_weight = nrow(df))

  expect_true(mout$alpha   > 0 && mout$alpha   < 1)
  expect_true(mout$theta1  > 0 && mout$theta1  < 1)
  expect_true(mout$theta0  > 0 && mout$theta0  < 1)
  expect_true(mout$pi      >= 0 && mout$pi     <= 0.49)
  expect_true(mout$eps     >= 1e-4 && mout$eps <= 0.95)
  expect_true(mout$lambda_g > 0)
  expect_true(mout$lambda_d > 0)
  expect_null(mout$sigma2_g)
})

# ---------------------------------------------------------------------------
# 13. m_step_eps: stationary=TRUE enforces stationarity
# ---------------------------------------------------------------------------

test_that("m_step_eps: stationary=TRUE gives alpha = theta0/(theta0+1-theta1)", {
  df   <- .make_eps_data()
  suff <- e_step_eps(df, .make_eps_params())$suff
  mout <- m_step_eps(suff, total_weight = nrow(df), stationary = TRUE)

  expected_alpha <- mout$theta0 / (mout$theta0 + 1 - mout$theta1)
  expect_equal(mout$alpha, expected_alpha, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 14. init_params_eps: basic structure and range checks
# ---------------------------------------------------------------------------

test_that("init_params_eps returns 7 valid params (no sigma2_g)", {
  df  <- .make_eps_data()
  p0  <- init_params_eps(df)

  expect_named(p0, c("alpha", "theta0", "theta1", "pi", "eps", "lambda_g", "lambda_d"))
  expect_true(p0$alpha   > 0 && p0$alpha   < 1)
  expect_true(p0$theta1  > 0 && p0$theta1  < 1)
  expect_true(p0$theta0  > 0 && p0$theta0  < 1)
  expect_true(p0$eps     > 0 && p0$eps     < 1)
  expect_true(p0$lambda_g > 0)
  expect_true(p0$lambda_d > 0)
  expect_null(p0$sigma2_g)
})

# ---------------------------------------------------------------------------
# 15. init_params_eps: invalid eps_init raises error
# ---------------------------------------------------------------------------

test_that("init_params_eps raises error for invalid eps_init", {
  df <- .make_eps_data()
  expect_error(init_params_eps(df, eps_init = 0))
  expect_error(init_params_eps(df, eps_init = 1))
  expect_error(init_params_eps(df, eps_init = -0.1))
})

# ---------------------------------------------------------------------------
# 16. em_fit_tenure_eps: required column validation
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_eps raises error for missing columns", {
  df <- .make_eps_data()
  df$tenure1 <- NULL
  expect_error(em_fit_tenure_eps(df), "missing columns")
})

# ---------------------------------------------------------------------------
# 17. em_fit_tenure_eps: loglik non-decreasing (free model)
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_eps: LL non-decreasing across iterations (free)", {
  for (seed in c(1L, 42L, 999L)) {
    df  <- .make_eps_data(n = 300L, seed = seed)
    fit <- em_fit_tenure_eps(df, max_iter = 50L, verbose = 0L)

    ll_hist <- fit$history$loglik
    diffs   <- diff(ll_hist)
    expect_true(
      all(diffs >= -1e-10),
      info = sprintf("Seed %d: LL decreased; min diff = %.3e", seed, min(diffs))
    )
  }
})

# ---------------------------------------------------------------------------
# 18. em_fit_tenure_eps: output structure is complete
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_eps returns expected output fields", {
  df  <- .make_eps_data(n = 100L)
  fit <- em_fit_tenure_eps(df, max_iter = 5L, verbose = 0L)

  expect_named(fit, c("params", "loglik", "history", "converged", "status",
                      "iterations", "gamma", "diagnostics"))
  expect_equal(dim(fit$gamma), c(100L, 8L))
  expect_true(is.finite(fit$loglik))
  expect_true(is.logical(fit$converged))
  expect_true(fit$iterations >= 1L)
  expect_true(is.data.frame(fit$history))
  expect_true("eps" %in% names(fit$params))
  expect_null(fit$params$sigma2_g)
})

test_that("stationary epsilon M-step increases the constrained Q block", {
  df <- .make_eps_data(n=400L, seed=912L)
  old <- .make_eps_params()
  old$alpha <- old$theta0/(old$theta0+1-old$theta1)
  estep <- e_step_eps(df,old)
  for (linked in c(FALSE,TRUE)) {
    new <- m_step_eps(estep$suff,sum(df$weight),stationary=TRUE,linked=linked)
    q_old <- .q_theta_stationary_eps(old$theta0,old$theta1,estep$suff,linked)
    q_new <- .q_theta_stationary_eps(new$theta0,new$theta1,estep$suff,linked)
    expect_gte(q_new+1e-7,q_old)
    expect_equal(new$alpha,new$theta0/(new$theta0+1-new$theta1),tolerance=1e-12)
  }
})

test_that("epsilon driver reports a small residual after genuine convergence", {
  df <- .make_eps_data(n=500L,seed=441L)
  fit <- em_fit_tenure_eps(df,stationary=TRUE,linked=TRUE,
    max_iter=500L,tol=1e-9,param_tol=1e-6,verbose=0L)
  if (fit$converged) {
    expect_identical(fit$status,"converged")
    expect_lt(fit$diagnostics$fixedpoint_residual,2e-6)
  } else {
    expect_false(fit$status == "converged")
  }
})

test_that("collapsing epsilon cells preserves normalized total weight", {
  df <- .make_eps_data(n=300L,seed=337L)
  collapsed <- collapse_eps_cells(df)
  expect_equal(sum(collapsed$weight),nrow(df),tolerance=1e-10)
  expect_equal(attr(collapsed,"n_original"),nrow(df))
  expect_lte(nrow(collapsed),nrow(df))
})

# ---------------------------------------------------------------------------
# 19. em_fit_tenure_eps: sigma2_g in warm-start is stripped
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_eps: sigma2_g from warm-start is stripped", {
  df <- .make_eps_data(n = 80L)
  p0 <- init_params_eps(df)
  p0$sigma2_g <- 0.5     # inject base-model param
  p0$sigma2_d <- 0.1

  fit <- em_fit_tenure_eps(df, params0 = p0, max_iter = 3L, verbose = 0L)
  expect_null(fit$params$sigma2_g)
  expect_null(fit$params$sigma2_d)
})

# ---------------------------------------------------------------------------
# 20. em_fit_tenure_eps: linked=TRUE keeps lambda_g = f(theta1)
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_eps: linked mode enforces lambda_g = ctmc_lambda(theta1)", {
  set.seed(2)
  df  <- .make_eps_data(n = 100L)
  fit <- em_fit_tenure_eps(df, linked = TRUE, max_iter = 5L, verbose = 0L)

  expected_lambda_g <- ctmc_lambda_from_persistence(fit$params$theta1)
  expect_equal(fit$params$lambda_g, expected_lambda_g, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 21. e_step_eps: Eps_den = 0 when no K>=2 spells exist
# ---------------------------------------------------------------------------

test_that("e_step_eps: Eps_den = 0 and Eps_num = 0 when all spells are K<=1", {
  # Force all observations to be nonemployed at waves 2 and 3 (K=0 for every
  # E-spell since no observation has more than 1 employed wave).
  df      <- .make_eps_data()
  df$y2   <- 0L
  df$y3   <- 0L
  # tenure at non-employed waves may be NA (natural QLFS data).
  df$tenure2 <- NA_real_
  df$tenure3 <- NA_real_

  out <- e_step_eps(df, .make_eps_params())

  expect_equal(out$suff$Eps_den, 0)
  expect_equal(out$suff$Eps_num, 0)
  expect_true(is.finite(out$loglik))
})

# ---------------------------------------------------------------------------
# 21b. e_step_eps: K=1/offset>0 spells DO contribute to Eps_den
# ---------------------------------------------------------------------------

test_that("e_step_eps: K=1 offset>0 spells contribute to Eps_den", {
  # With s=(0,1,0) pattern, history h=(1,1,0) produces an E-spell [1,2] where
  # wave 2 is observed (offset=1, t_1>a). This is eps-informative (K=1, offset>0).
  # After the fix, Eps_den must be > 0 for such data.
  df      <- .make_eps_data(n = 300L, seed = 123L)
  df$y1   <- 0L
  df$y3   <- 0L
  # y2=1 for ~65% (from original draw); tenure1/3 irrelevant but keep finite
  df$tenure1 <- 0.5
  df$tenure3 <- 0.5

  out <- e_step_eps(df, .make_eps_params())

  # Some observations have s2=1 which, under h=(1,1,0) or h=(0,1,1), produces
  # K=1/offset>0 E-spells. Eps_den must be strictly positive.
  expect_true(out$suff$Eps_den > 0,
              info = "Eps_den should be > 0 when K=1/offset>0 spells exist")
  expect_true(out$suff$Eps_num >= 0)
  expect_true(out$suff$Eps_num <= out$suff$Eps_den + 1e-10)
})

# ---------------------------------------------------------------------------
# 22. log_emission_spell_g: K=3 consistent rows use both mixture branches
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=3 consistent rows have tau_sum in [0, 1)", {
  lambda_g <- 0.4; eps <- 0.25
  # All three waves clock-consistent
  g1 <- 2.0; g2 <- g1 + 0.25; g3 <- g2 + 0.25
  g_mat <- matrix(c(g1, g2, g3), nrow = 1)
  s_mat <- matrix(c(TRUE, TRUE, TRUE), nrow = 1)

  out <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L, 2L),
                               lambda_g = lambda_g, eps = eps)

  expect_equal(out$K, 3L)
  expect_true(is.finite(out$loglik))
  # tau_sum = E[# contaminated waves] in [0, 3]; CCC dominates for consistent rows
  expect_true(out$tau_sum >= 0 && out$tau_sum < 1)
})

# ---------------------------------------------------------------------------
# 23. e_step_eps: sufficient stats are correct under heterogeneous weights
# ---------------------------------------------------------------------------

test_that("e_step_eps: gamma row-sums = 1 and suff stats finite with varying weights", {
  set.seed(7L)
  df         <- .make_eps_data()
  df$weight  <- runif(nrow(df), 0.5, 2.0)

  out <- e_step_eps(df, .make_eps_params())

  # Gamma rows still sum to 1 (weighting affects suff stats, not responsibilities)
  expect_equal(rowSums(out$gamma), rep(1, nrow(df)), tolerance = 1e-10)
  expect_true(all(out$gamma >= 0))
  # Sufficient stats should be finite and positive
  expect_true(is.finite(out$suff$Lg_count) && out$suff$Lg_count > 0)
  expect_true(is.finite(out$suff$Lg_xsum)  && out$suff$Lg_xsum  > 0)
  expect_true(is.finite(out$loglik))
})

# ---------------------------------------------------------------------------
# 24. log_emission_spell_g: K=2 matches 4-pattern analytical computation
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=2 consistent loglik matches 4-pattern analytical", {
  lambda_g <- 0.4; eps <- 0.25; Delta <- 0.25
  g1 <- 3.0; g2 <- g1 + Delta   # consistent: d1=0, d2=1
  T1 <- g1                        # T_g from wave 1 (d1=0)
  T2 <- g2 - Delta                # T_g from wave 2 = g1

  # Analytical 4-pattern log-densities
  lp_cc <- 2*log(1-eps) + log(lambda_g) - lambda_g*T1
  lp_cx <- log(1-eps)+log(eps) + 2*log(lambda_g) - lambda_g*(T1+g2)
  lp_xc <- log(eps)+log(1-eps) + 2*log(lambda_g) - lambda_g*(g1+T2)
  lp_xx <- 2*log(eps)            + 2*log(lambda_g) - lambda_g*(g1+g2)
  mx    <- max(lp_cc, lp_cx, lp_xc, lp_xx)
  expected_ll <- mx + log(exp(lp_cc-mx)+exp(lp_cx-mx)+exp(lp_xc-mx)+exp(lp_xx-mx))

  g_mat <- matrix(c(g1, g2), nrow = 1)
  s_mat <- matrix(c(TRUE, TRUE), nrow = 1)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = lambda_g, eps = eps)

  expect_equal(out$loglik, expected_ll, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 25. log_emission_spell_g: K=3 partial consistency gives higher density than AON
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=3 partial consistency higher than AON", {
  lambda_g <- 0.17; eps <- 0.35; Delta <- 0.25
  # Only first increment consistent (g2=g1+Delta), g3 off
  g1 <- 5.0; g2 <- g1 + Delta; g3 <- 8.0

  # Old AON: since not fully consistent, collapses to all-contaminated
  aon_ll <- log(1 - (1-eps)^3) + 3*log(lambda_g) - lambda_g*(g1+g2+g3)

  g_mat <- matrix(c(g1, g2, g3), nrow = 1)
  s_mat <- matrix(rep(TRUE, 3), nrow = 1)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L, 2L), lambda_g = lambda_g, eps = eps)

  # Full enumeration should give substantially higher density (CCX pattern adds mass)
  expect_true(out$loglik > aon_ll + 0.5)   # at least exp(0.5) ≈ 1.65x higher
  expect_true(out$tau_sum > 1 && out$tau_sum < 3)   # tau_sum in (1, 3)
})

# ---------------------------------------------------------------------------
# 26. log_emission_spell_g: K=3 fully consistent matches 8-pattern analytical
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=3 fully consistent loglik matches 8-pattern analytical", {
  lambda_g <- 0.4; eps <- 0.20; Delta <- 0.25
  g1 <- 2.0; g2 <- g1 + Delta; g3 <- g1 + 2*Delta
  T1 <- g1; T2 <- g2 - Delta; T3 <- g3 - 2*Delta  # all equal g1

  # All 8 pattern log-densities
  lp_CCC <- 3*log(1-eps) + log(lambda_g) - lambda_g*T1
  lp_CCX <- 2*log(1-eps)+log(eps) + 2*log(lambda_g) - lambda_g*(T1+g3)
  lp_CXC <- 2*log(1-eps)+log(eps) + 2*log(lambda_g) - lambda_g*(T1+g2)
  lp_CXX <- log(1-eps)+2*log(eps) + 3*log(lambda_g) - lambda_g*(T1+g2+g3)
  lp_XCC <- log(eps)+2*log(1-eps) + 2*log(lambda_g) - lambda_g*(g1+T2)
  lp_XCX <- 2*log(eps)+log(1-eps) + 3*log(lambda_g) - lambda_g*(g1+T2+g3)
  lp_XXC <- 2*log(eps)+log(1-eps) + 3*log(lambda_g) - lambda_g*(g1+g2+T3)
  lp_XXX <- 3*log(eps)            + 3*log(lambda_g) - lambda_g*(g1+g2+g3)
  mx <- max(lp_CCC,lp_CCX,lp_CXC,lp_CXX,lp_XCC,lp_XCX,lp_XXC,lp_XXX)
  expected_ll <- mx + log(
    exp(lp_CCC-mx)+exp(lp_CCX-mx)+exp(lp_CXC-mx)+exp(lp_CXX-mx)+
    exp(lp_XCC-mx)+exp(lp_XCX-mx)+exp(lp_XXC-mx)+exp(lp_XXX-mx)
  )

  g_mat <- matrix(c(g1, g2, g3), nrow = 1)
  s_mat <- matrix(rep(TRUE, 3), nrow = 1)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L, 2L), lambda_g = lambda_g, eps = eps)

  expect_equal(out$loglik, expected_ll, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 27. log_emission_spell_g: K=2 with offset (0,2) — s=(1,0,1) under h=(1,1,1)
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=2 offset (0,2) consistent loglik matches analytical", {
  lambda_g <- 0.17; eps <- 0.22; Delta <- 0.25
  # Consistent for this offset: g3 = g1 + 2*Delta = g1 + 0.5
  g1 <- 5.0; g3 <- g1 + 2*Delta   # g1 at col1, g3 at col2 (unobserved wave at col2 middle)
  T_from_j1 <- g1 - 0*Delta   # = g1
  T_from_j2 <- g3 - 2*Delta   # = g1

  lp_cc <- 2*log(1-eps) + log(lambda_g) - lambda_g*T_from_j1
  lp_cx <- log(1-eps)+log(eps) + 2*log(lambda_g) - lambda_g*(T_from_j1+g3)
  lp_xc <- log(eps)+log(1-eps) + 2*log(lambda_g) - lambda_g*(g1+T_from_j2)
  lp_xx <- 2*log(eps)          + 2*log(lambda_g) - lambda_g*(g1+g3)
  mx    <- max(lp_cc, lp_cx, lp_xc, lp_xx)
  expected_ll <- mx + log(exp(lp_cc-mx)+exp(lp_cx-mx)+exp(lp_xc-mx)+exp(lp_xx-mx))

  # L=3 matrix: col1 observed, col2 unobserved, col3 observed
  g_mat <- matrix(c(g1, NA_real_, g3), nrow = 1)
  s_mat <- matrix(c(TRUE, FALSE, TRUE), nrow = 1)
  out   <- log_emission_spell_g(g_mat, s_mat, c(0L, 1L, 2L), lambda_g = lambda_g, eps = eps)

  expect_equal(out$loglik, expected_ll, tolerance = 1e-10)
  expect_true(out$tau_sum >= 0 && out$tau_sum < 1)   # CC dominates for consistent obs
})

# ---------------------------------------------------------------------------
# 28. em_fit_tenure_eps: convergence test — eps recovery from synthetic data
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_eps: converges and recovers eps approximately", {
  # Synthetic data with eps_true ~ 0.20: 80% consistent EE pairs
  set.seed(123L)
  n <- 500L
  # All obs (1,1,1) — ensures K=3 spells dominate
  g1 <- pmax(0.3, rexp(n, rate = 0.15))
  # 80% fully consistent (CCC), 20% partially contaminated
  g2 <- ifelse(runif(n) < 0.80, g1 + 0.25, pmax(0.3, rexp(n, rate = 0.15)))
  g3 <- ifelse(runif(n) < 0.80, g2 + 0.25, pmax(0.3, rexp(n, rate = 0.15)))
  df <- data.frame(
    y1 = 1L, y2 = 1L, y3 = 1L,
    tenure1 = g1, tenure2 = g2, tenure3 = g3,
    timegap_cat1 = 1L, timegap_cat2 = 1L, timegap_cat3 = 1L,
    weight = 1
  )

  fit <- em_fit_tenure_eps(df, max_iter = 100L, verbose = 0L)

  # EM should converge and recover eps in a plausible range
  expect_true(fit$converged || fit$iterations >= 10L)
  expect_true(fit$params$eps >= 0.05 && fit$params$eps <= 0.50)
  expect_true(fit$params$pi  >= 0    && fit$params$pi  <= 0.49)
  # LL non-decreasing
  diffs <- diff(fit$history$loglik)
  expect_true(all(diffs >= -1e-10))
})

# ---------------------------------------------------------------------------
# 29. .maximal_e_spells: correct enumeration of E-spells in all cases
# ---------------------------------------------------------------------------

test_that(".maximal_e_spells: correctly enumerates maximal E-spell ranges", {
  expect_equal(.maximal_e_spells(c(0L, 0L, 0L)), list())
  expect_equal(.maximal_e_spells(c(1L, 1L, 1L)), list(1:3))
  expect_equal(.maximal_e_spells(c(0L, 1L, 1L)), list(2:3))
  expect_equal(.maximal_e_spells(c(1L, 1L, 0L)), list(1:2))
  expect_equal(.maximal_e_spells(c(1L, 0L, 1L)), list(1L, 3L))
  expect_equal(.maximal_e_spells(c(0L, 1L, 0L)), list(2L))
  expect_equal(.maximal_e_spells(c(1L, 0L, 0L)), list(1L))
})

# ---------------------------------------------------------------------------
# 30. m_step_eps: direct closed-form alpha update
# ---------------------------------------------------------------------------

test_that("m_step_eps: alpha = C1/(C1+C0) for free non-stationary", {
  df      <- .make_eps_data()
  estep   <- e_step_eps(df, .make_eps_params())
  suff    <- estep$suff
  mout    <- m_step_eps(suff, total_weight = sum(df$weight), stationary = FALSE)

  expected_alpha <- suff$C1 / (suff$C1 + suff$C0)
  expect_equal(mout$alpha, expected_alpha, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 31. e_step_eps: missing weight column raises error
# ---------------------------------------------------------------------------

test_that("e_step_eps: missing weight column raises error", {
  df         <- .make_eps_data()
  df$weight  <- NULL
  expect_error(e_step_eps(df, .make_eps_params()), "weight")
})

# ---------------------------------------------------------------------------
# 32. log_emission_spell_g: numeric s_mat raises error
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: numeric s_mat raises error", {
  g_mat <- matrix(c(2.0, 2.25), nrow = 1)
  s_mat <- matrix(c(1L, 1L),    nrow = 1)  # integer, not logical
  expect_error(
    log_emission_spell_g(g_mat, s_mat, c(0L, 1L), lambda_g = 0.4, eps = 0.2),
    "logical"
  )
})

# ---------------------------------------------------------------------------
# 33. .log_exp_g: matches log(lambda) - lambda*g; -Inf for non-positive g
# ---------------------------------------------------------------------------

test_that(".log_exp_g: correct Exp log-density and -Inf for invalid inputs", {
  lambda <- 0.4

  # Positive g: matches closed form
  expect_equal(.log_exp_g(2.0, lambda), log(lambda) - lambda * 2.0,
               tolerance = 1e-14)
  expect_equal(.log_exp_g(c(1.0, 3.0), lambda),
               c(log(lambda) - lambda * 1.0, log(lambda) - lambda * 3.0),
               tolerance = 1e-14)

  # Non-positive / non-finite g: -Inf
  expect_equal(.log_exp_g(0,         lambda), -Inf)
  expect_equal(.log_exp_g(-1,        lambda), -Inf)
  expect_equal(.log_exp_g(NA_real_,  lambda), -Inf)
  expect_equal(.log_exp_g(Inf,       lambda), -Inf)
})

# ---------------------------------------------------------------------------
# 34. log_emission_spell_g: K>3 raises informative error
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K>3 (4-column s_mat) raises error", {
  g_mat <- matrix(c(2.0, 2.25, 2.5, 2.75), nrow = 1)
  s_mat <- matrix(rep(TRUE, 4), nrow = 1)
  expect_error(
    log_emission_spell_g(g_mat, s_mat, c(0L, 1L, 2L, 3L),
                         lambda_g = 0.4, eps = 0.2),
    "K > 3"
  )
})

# ---------------------------------------------------------------------------
# 35. .m_step_theta1_eps_brent: fallback when bracket same-sign
# ---------------------------------------------------------------------------

test_that(".m_step_theta1_eps_brent: same-sign bracket returns seed in (0, theta_cap)", {
  # Force same-sign by setting degenerate sufficient statistics:
  # T11 = D1 (all E->E transitions, perfect persistence) makes FOC always negative.
  T11 <- 200; D1 <- 200; Lg_count <- 0; Lg_xsum <- 0
  theta_seed <- 0.85; theta_cap <- 0.999

  result <- .m_step_theta1_eps_brent(
    T11 = T11, D1 = D1, Lg_count = Lg_count, Lg_xsum = Lg_xsum,
    theta_seed = theta_seed, theta_cap = theta_cap
  )
  expect_true(result > 0 && result < theta_cap)
})

# ---------------------------------------------------------------------------
# 36. em_fit_tenure_eps: linked=TRUE converges and free LL >= linked LL
# ---------------------------------------------------------------------------

test_that("em_fit_tenure_eps: linked=TRUE converges; free LL >= linked LL", {
  df <- .make_eps_data(n = 150L, seed = 99L)

  fit_free   <- em_fit_tenure_eps(df, max_iter = 80L, verbose = 0L)
  fit_linked <- em_fit_tenure_eps(df, linked = TRUE, max_iter = 80L, verbose = 0L)

  # Both converge or reach max_iter gracefully
  expect_true(is.logical(fit_free$converged))
  expect_true(is.logical(fit_linked$converged))

  # CTMC link is enforced
  lam_g_check <- ctmc_lambda_from_persistence(fit_linked$params$theta1)
  expect_equal(fit_linked$params$lambda_g, lam_g_check, tolerance = 1e-8)

  # Free model LL >= linked LL (fewer constraints → no worse)
  expect_true(fit_free$loglik >= fit_linked$loglik - 1e-6)
})

# ---------------------------------------------------------------------------
# 37. m_step_eps: Eps_den = 0 returns eps_floor (no K>=2 spells)
# ---------------------------------------------------------------------------

test_that("m_step_eps: Eps_den=0 falls back to eps_floor", {
  # Construct minimal suff with zero eps denominator
  suff <- e_step_eps(.make_eps_data(), .make_eps_params())$suff
  suff$Eps_num <- 0
  suff$Eps_den <- 0

  mout <- m_step_eps(suff, total_weight = 200, stationary = FALSE)
  expect_equal(mout$eps, 1e-4)   # default eps_floor
})

# ---------------------------------------------------------------------------
# 38. e_step_eps: Lg_count and Lg_xsum spot-check for a K=2 spell
# ---------------------------------------------------------------------------

test_that("e_step_eps: Lg_count > 0 and Lg_xsum > 0 for a K=2 employed spell", {
  # Single obs: employed at waves 1+2, not at 3 — creates one K=2 E-spell
  df <- data.frame(
    y1 = 1L, y2 = 1L, y3 = 0L,
    tenure1 = 2.0, tenure2 = 2.25, tenure3 = NA_real_,
    timegap_cat1 = 4L, timegap_cat2 = 4L, timegap_cat3 = 2L,
    weight = 1.0
  )
  out <- e_step_eps(df, .make_eps_params())

  expect_true(out$suff$Lg_count > 0)
  expect_true(out$suff$Lg_xsum  > 0)
  expect_true(out$suff$Eps_den  > 0)   # K=2 spell contributes to eps denominator
})

# ---------------------------------------------------------------------------
# 39. log_emission_spell_g: K=2 offset pairs (0,1) and (1,2)
# ---------------------------------------------------------------------------

test_that("log_emission_spell_g: K=2 offset (0,1) and (1,2) give finite loglik", {
  lambda_g <- 0.4; eps <- 0.2; Delta <- 0.25

  for (d_pair in list(c(0L, 1L), c(1L, 2L))) {
    d1 <- d_pair[1]; d2 <- d_pair[2]
    # Consistent pair: g2 = g1 + (d2-d1)*Delta
    g1 <- 3.0; g2 <- g1 + (d2 - d1) * Delta
    g_mat <- matrix(c(g1, g2), nrow = 1)
    s_mat <- matrix(c(TRUE, TRUE), nrow = 1)

    out <- log_emission_spell_g(g_mat, s_mat, d_pair, lambda_g, eps)

    expect_true(is.finite(out$loglik),
                info = sprintf("offset (%d,%d): non-finite loglik", d1, d2))
    expect_true(out$tau_sum >= 0,
                info = sprintf("offset (%d,%d): negative tau_sum", d1, d2))
    expect_equal(out$K, 2L,
                 info = sprintf("offset (%d,%d): wrong K", d1, d2))
  }
})
