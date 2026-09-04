test_that("monthly duration masses form a proper distribution", {
  mass <- exp(.log_duration_month_mass(0:5000, rep(.2, 5), 0))
  expect_equal(sum(mass), 1, tolerance = 1e-10)
  expect_true(all(mass >= 0))
})

test_that("monthly clean clocks and within-quarter starts are enforced", {
  s <- matrix(TRUE, nrow=3, ncol=3)
  g <- rbind(c(0, .25, .50), c(1, 1.25, 1.50),
    c(.25, .50, .75))
  marginal <- log_emission_spell_g_monthly(g, s, 0:2,
    rep(.2, 5), .2, initial_model="marginal")
  expect_true(all(is.finite(marginal$loglik)))
  entry <- log_emission_spell_g_monthly(g, s, 0:2,
    rep(.2, 5), 1e-8, initial_model="within_interval")
  expect_true(is.finite(entry$loglik[1]))
  expect_gt(entry$loglik[1], entry$loglik[2] + 20)
  expect_gt(entry$loglik[1], entry$loglik[3] + 20)
})

test_that("monthly reset emission recovers its generating probability", {
  set.seed(24019)
  n <- 3000L
  q <- .06
  eps <- .15
  m <- matrix(NA_integer_, n, 3)
  m[, 1] <- floor(12 * rexp(n, .2))
  for (j in 2:3) {
    reset <- runif(n) < q
    m[, j] <- ifelse(reset, sample(0:2, n, replace=TRUE),
      m[, j-1] + 3L)
  }
  gross <- matrix(runif(3*n) < eps, n, 3)
  m[gross] <- floor(12 * rexp(sum(gross), .2))
  criterion <- function(value) -sum(log_emission_spell_g_monthly(
    m/12, matrix(TRUE, n, 3), 0:2, rep(.2, 5), eps,
    job_change_prob=value)$loglik)
  q_hat <- optimize(criterion, c(1e-5, .20), tol=1e-6)$minimum
  expect_equal(q_hat, q, tolerance=.015)
})

test_that("monthly simulator preserves zero-month tenure reports", {
  p <- list(alpha=.55, pi=.03, eps=.20, eps_d=.12,
    job_change_prob=.04, lambda_g=rep(.20, 5), lambda_d=rep(.30, 5),
    tenure_measurement_model="monthly")
  d <- simulate_eps_piecewise_job_change(3000, p, seed=24020)
  expect_silent(validate_df_eps(d, allow_zero_tenure=TRUE))
  expect_true(any(unlist(d[paste0("tenure", 1:3)]) == 0, na.rm=TRUE))
  expect_true(all(abs(12 * unlist(d[paste0("tenure", 1:3)]) -
    round(12 * unlist(d[paste0("tenure", 1:3)]))) < 1e-8,
    na.rm=TRUE))
})

test_that("persistent monthly reporting exactly nests independent gross draws", {
  g <- matrix(c(24,27,30, 8,19,22, 3,6,15)/12, nrow=3, byrow=TRUE)
  s <- matrix(TRUE,3,3)
  lambda <- c(.20,.16,.12,.09,.06)
  old <- log_emission_spell_g_monthly(g,s,0:2,lambda,.24,
    job_change_prob=.03,tenure_report_persistence=0)
  nested <- log_emission_spell_g_monthly(g,s,0:2,lambda,.24,
    job_change_prob=.03)
  expect_equal(old$loglik,nested$loglik,tolerance=1e-12)
  expect_equal(old$tau_sum,nested$tau_sum,tolerance=1e-12)
  expect_equal(old$job_changes,nested$job_changes,tolerance=1e-12)
})

test_that("persistence rewards retention of a contaminated start-date anchor", {
  consistent <- matrix(c(24,27,30)/12,1,3)
  revised <- matrix(c(24,41,18)/12,1,3)
  s <- matrix(TRUE,1,3)
  lambda <- c(.20,.16,.12,.09,.06)
  ll0_c <- log_emission_spell_g_monthly(consistent,s,0:2,lambda,.80,
    tenure_report_persistence=0)$loglik
  ll8_c <- log_emission_spell_g_monthly(consistent,s,0:2,lambda,.80,
    tenure_report_persistence=.8)$loglik
  ll0_r <- log_emission_spell_g_monthly(revised,s,0:2,lambda,.80,
    tenure_report_persistence=0)$loglik
  ll8_r <- log_emission_spell_g_monthly(revised,s,0:2,lambda,.80,
    tenure_report_persistence=.8)$loglik
  expect_gt(ll8_c-ll0_c,ll8_r-ll0_r)
  expect_true(all(is.finite(c(ll0_c,ll8_c,ll0_r,ll8_r))))
})

test_that("January-heaped monthly masses normalize for every interview month", {
  lambda <- c(.20,.16,.12,.09,.06)
  log_norm <- .log_january_duration_normalizers(lambda)
  for (interview_month in 1:12) {
    residue <- interview_month-1L
    support <- residue+12L*(0:2500)
    mass <- exp(.log_january_duration_month_mass(support,
      rep(interview_month,length(support)),lambda,0,log_norm))
    expect_equal(sum(mass),1,tolerance=1e-9)
  }
})

test_that("calendar-weighted monthly masses normalize and nest uniform months", {
  lambda <- c(.20,.16,.12,.09,.06)
  support <- 0:30000
  uniform <- rep(1/12,12)
  seasonal <- c(.22,.05,.05,.06,.06,.06,.07,.07,.07,.08,.09,.12)
  for (interview_month in 1:12) {
    base <- .log_duration_month_mass(support,lambda)
    nested <- .log_calendar_duration_month_mass(support,
      rep(interview_month,length(support)),lambda,
      start_month_probs=uniform)
    weighted <- .log_calendar_duration_month_mass(support,
      rep(interview_month,length(support)),lambda,
      start_month_probs=seasonal)
    expect_equal(nested,base,tolerance=1e-12)
    expect_equal(sum(exp(weighted)),1,tolerance=1e-9)
  }
})

test_that("start-month distribution is shared by clean and gross reports", {
  lambda <- c(.20,.16,.12,.09,.06)
  interview <- matrix(3L,1,1)
  s <- matrix(TRUE,1,1)
  january <- matrix(26/12,1,1)
  february <- matrix(25/12,1,1)
  seasonal <- c(.30,rep(.70/11,11))
  clean_j <- log_emission_spell_g_monthly(january,s,0,lambda,1e-8,
    tenure_start_month_probs=seasonal,
    interview_month_mat=interview)$loglik
  clean_f <- log_emission_spell_g_monthly(february,s,0,lambda,1e-8,
    tenure_start_month_probs=seasonal,
    interview_month_mat=interview)$loglik
  gross_j <- log_emission_spell_g_monthly(january,s,0,lambda,1-1e-8,
    tenure_start_month_probs=seasonal,
    interview_month_mat=interview)$loglik
  gross_f <- log_emission_spell_g_monthly(february,s,0,lambda,1-1e-8,
    tenure_start_month_probs=seasonal,
    interview_month_mat=interview)$loglik
  expect_gt(clean_j-clean_f,1)
  expect_gt(gross_j-gross_f,1)
})

test_that("calendar baseline pack and unpack exactly nests uniform months", {
  p <- list(alpha=.55,theta0=.02,theta1=.04,pi=.03,eps=.20,eps_d=.12,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5))
  z <- .piecewise_calendar_revision_monthly_pack(p,.20,.10,.04,.08,
    rep(1/12,12))
  restored <- .piecewise_calendar_revision_monthly_unpack(z)
  expect_equal(restored$tenure_start_month_probs,rep(1/12,12),
    tolerance=1e-12)
})

test_that("calendar-heaping extension exactly nests independent gross reports", {
  g <- matrix(c(2,2.25,2.5,8/12,19/12,22/12),nrow=2,byrow=TRUE)
  s <- matrix(TRUE,2,3)
  interview <- rbind(c(3,6,9),c(6,9,12))
  lambda <- c(.20,.16,.12,.09,.06)
  base <- log_emission_spell_g_monthly(g,s,0:2,lambda,.24,
    job_change_prob=.03)
  nested <- log_emission_spell_g_monthly(g,s,0:2,lambda,.24,
    job_change_prob=.03,tenure_heaping_prob=0,
    interview_month_mat=interview)
  expect_equal(base$loglik,nested$loglik,tolerance=1e-12)
  expect_equal(base$tau_sum,nested$tau_sum,tolerance=1e-12)
})

test_that("calendar heaping favors tenure reports implying January starts", {
  # A March interview with tenure 26 months implies a January start; tenure
  # 25 months implies February.  At high eps, heaping must favor the former.
  january <- matrix(26/12,1,1)
  february <- matrix(25/12,1,1)
  s <- matrix(TRUE,1,1)
  interview <- matrix(3L,1,1)
  lambda <- c(.20,.16,.12,.09,.06)
  ll0_j <- log_emission_spell_g_monthly(january,s,0,lambda,.999,
    tenure_heaping_prob=0,interview_month_mat=interview)$loglik
  ll0_f <- log_emission_spell_g_monthly(february,s,0,lambda,.999,
    tenure_heaping_prob=0,interview_month_mat=interview)$loglik
  ll8_j <- log_emission_spell_g_monthly(january,s,0,lambda,.999,
    tenure_heaping_prob=.8,interview_month_mat=interview)$loglik
  ll8_f <- log_emission_spell_g_monthly(february,s,0,lambda,.999,
    tenure_heaping_prob=.8,interview_month_mat=interview)$loglik
  expect_gt(ll8_j-ll0_j,ll8_f-ll0_f)
})

test_that("monthly simulator generates calendar heaping when requested", {
  p <- list(alpha=.55,pi=.03,eps=.40,eps_d=.12,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5),
    tenure_measurement_model="monthly",tenure_heaping_prob=.70)
  d <- simulate_eps_piecewise_job_change(4000,p,seed=24021)
  employed <- unlist(d[paste0("y",1:3)])==1L
  months <- round(12*unlist(d[paste0("tenure",1:3)]))[employed]
  interviews <- unlist(d[paste0("interview_month",1:3)])[employed]
  january <- months%%12L==interviews-1L
  expect_gt(mean(january),.20)
})

test_that("whole-year revision masses form proper conditional distributions", {
  lambda <- c(.20,.16,.12,.09,.06)
  support <- 0:30000
  log_mass <- .log_whole_year_revision_month_mass(support,
    rep(26L,length(support)),rep(6L,length(support)),
    rep(3L,length(support)),rep(1L,length(support)),lambda)
  expect_equal(sum(exp(log_mass)),1,tolerance=1e-9)
  expect_false(is.finite(log_mass[29L+1L]))
  expect_true(all(support[is.finite(log_mass)]%%12L==5L))
})

test_that("whole-year revision extension nests and rewards revised years", {
  lambda <- c(.20,.16,.12,.09,.06)
  s <- matrix(TRUE,1,2)
  interview <- matrix(c(3L,6L),1,2)
  yearly <- matrix(c(26L,41L)/12,1,2)
  other_month <- matrix(c(26L,40L)/12,1,2)
  base <- log_emission_spell_g_monthly(yearly,s,0:1,lambda,.999,
    tenure_heaping_prob=.2,interview_month_mat=interview)
  nested <- log_emission_spell_g_monthly(yearly,s,0:1,lambda,.999,
    tenure_heaping_prob=.2,tenure_year_revision_prob=0,
    interview_month_mat=interview)
  expect_equal(base$loglik,nested$loglik,tolerance=1e-12)
  ll0_y <- base$loglik
  ll7_y <- log_emission_spell_g_monthly(yearly,s,0:1,lambda,.999,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.7,
    interview_month_mat=interview)$loglik
  ll0_o <- log_emission_spell_g_monthly(other_month,s,0:1,lambda,.999,
    tenure_heaping_prob=.2,interview_month_mat=interview)$loglik
  ll7_o <- log_emission_spell_g_monthly(other_month,s,0:1,lambda,.999,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.7,
    interview_month_mat=interview)$loglik
  expect_gt(ll7_y-ll0_y,ll7_o-ll0_o)
})

test_that("clean-anchor revision extension nests the gross-anchor model", {
  lambda <- c(.20,.16,.12,.09,.06)
  s <- matrix(TRUE,2,2)
  interview <- matrix(c(3L,6L,3L,6L),2,2,byrow=TRUE)
  g <- matrix(c(26L,41L,26L,40L)/12,2,2,byrow=TRUE)
  base <- log_emission_spell_g_monthly(g,s,0:1,lambda,.20,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.3,
    interview_month_mat=interview)
  nested <- log_emission_spell_g_monthly(g,s,0:1,lambda,.20,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.3,
    tenure_clean_anchor_revision_prob=0,
    interview_month_mat=interview)
  expect_equal(base$loglik,nested$loglik,tolerance=1e-12)
  expect_equal(base$tau_sum,nested$tau_sum,tolerance=1e-12)
})

test_that("clean-anchor revisions favor same-month different-year reports", {
  lambda <- c(.20,.16,.12,.09,.06)
  s <- matrix(TRUE,1,2)
  interview <- matrix(c(3L,6L),1,2)
  yearly <- matrix(c(26L,41L)/12,1,2)
  other_month <- matrix(c(26L,40L)/12,1,2)
  ll0_y <- log_emission_spell_g_monthly(yearly,s,0:1,lambda,.05,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=0,
    interview_month_mat=interview)$loglik
  ll7_y <- log_emission_spell_g_monthly(yearly,s,0:1,lambda,.05,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.7,
    interview_month_mat=interview)$loglik
  ll0_o <- log_emission_spell_g_monthly(other_month,s,0:1,lambda,.05,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=0,
    interview_month_mat=interview)$loglik
  ll7_o <- log_emission_spell_g_monthly(other_month,s,0:1,lambda,.05,
    tenure_heaping_prob=.2,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.7,
    interview_month_mat=interview)$loglik
  expect_gt(ll7_y-ll0_y,ll7_o-ll0_o)
})

test_that("monthly simulator generates whole-year revisions when requested", {
  p <- list(alpha=.95,pi=.001,eps=.95,eps_d=.12,
    job_change_prob=.001,lambda_g=rep(.20,5),lambda_d=rep(.30,5),
    tenure_measurement_model="monthly",tenure_heaping_prob=.15,
    tenure_year_revision_prob=.80)
  d <- simulate_eps_piecewise_job_change(4000,p,seed=24022)
  m1 <- round(12*d$tenure1); m2 <- round(12*d$tenure2)
  eligible <- d$y1==1L & d$y2==1L & is.finite(m1) & is.finite(m2)
  whole_year_revision <- (m2-m1-3L)%%12L==0L & m2!=m1+3L
  expect_gt(mean(whole_year_revision[eligible]),.30)
})

test_that("exact-anchor retention extension nests calendar revisions", {
  lambda <- c(.20,.16,.12,.09,.06)
  s <- matrix(TRUE,2,3)
  interview <- matrix(c(3L,6L,9L,6L,9L,12L),2,3,byrow=TRUE)
  g <- matrix(c(26L,29L,44L,9L,20L,23L)/12,2,3,byrow=TRUE)
  base <- log_emission_spell_g_monthly(g,s,0:2,lambda,.30,
    tenure_heaping_prob=.1,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.15,
    interview_month_mat=interview)
  nested <- log_emission_spell_g_monthly(g,s,0:2,lambda,.30,
    tenure_heaping_prob=.1,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.15,
    tenure_exact_anchor_retention_prob=0,
    interview_month_mat=interview)
  expect_equal(base$loglik,nested$loglik,tolerance=1e-12)
  expect_equal(base$tau_sum,nested$tau_sum,tolerance=1e-12)
})

test_that("exact-anchor retention favors exact continuations", {
  lambda <- c(.20,.16,.12,.09,.06)
  s <- matrix(TRUE,1,2)
  interview <- matrix(c(3L,6L),1,2)
  exact <- matrix(c(26L,29L)/12,1,2)
  revised_year <- matrix(c(26L,41L)/12,1,2)
  other <- matrix(c(26L,40L)/12,1,2)
  evaluate <- function(g,rho) log_emission_spell_g_monthly(
    g,s,0:1,lambda,.999,tenure_heaping_prob=.1,
    tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.2,
    tenure_exact_anchor_retention_prob=rho,
    interview_month_mat=interview)$loglik
  gain_exact <- evaluate(exact,.7)-evaluate(exact,0)
  gain_year <- evaluate(revised_year,.7)-evaluate(revised_year,0)
  gain_other <- evaluate(other,.7)-evaluate(other,0)
  expect_gt(gain_exact,gain_year)
  expect_gt(gain_exact,gain_other)
})

test_that("monthly simulator generates exact-anchor retention", {
  p <- list(alpha=.95,pi=.001,eps=.95,eps_d=.12,
    job_change_prob=.001,lambda_g=rep(.20,5),lambda_d=rep(.30,5),
    tenure_measurement_model="monthly",tenure_heaping_prob=.10,
    tenure_year_revision_prob=.10,
    tenure_clean_anchor_revision_prob=.10,
    tenure_exact_anchor_retention_prob=.80)
  d <- simulate_eps_piecewise_job_change(4000,p,seed=24023)
  m1 <- round(12*d$tenure1); m2 <- round(12*d$tenure2)
  eligible <- d$y1==1L & d$y2==1L & is.finite(m1) & is.finite(m2)
  expect_gt(mean(m2[eligible]==m1[eligible]+3L),.60)
})

test_that("calendar pack and unpack preserves exact retention", {
  p <- list(alpha=.55,theta0=.02,theta1=.04,pi=.03,eps=.20,eps_d=.12,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5))
  z <- .piecewise_calendar_revision_monthly_pack(p,.20,.10,.04,.08,
    rep(1/12,12),.25)
  restored <- .piecewise_calendar_revision_monthly_unpack(z)
  expect_equal(restored$tenure_exact_anchor_retention_prob,.25,
    tolerance=1e-12)
})

test_that("local revision kernel is normalized and excludes exact retention", {
  support <- .TENURE_LOCAL_MONTHS
  expected <- c(2L,24L)
  for (anchor in expected) {
    current <- anchor+support
    keep <- current>=0L
    lp <- .log_local_revision_month_mass(current[keep],
      rep(anchor,sum(keep)),0L)
    expect_equal(sum(exp(lp)),1,tolerance=1e-12)
    expect_equal(.log_local_revision_month_mass(anchor,anchor,0L),-Inf)
  }
})

test_that("local revision extension nests exact-anchor model at zero", {
  lambda <- c(.20,.16,.12,.09,.06)
  s <- matrix(TRUE,2,3)
  interview <- matrix(c(3L,6L,9L,6L,9L,12L),2,3,byrow=TRUE)
  g <- matrix(c(26L,31L,34L,9L,12L,17L)/12,2,3,byrow=TRUE)
  base <- log_emission_spell_g_monthly(g,s,0:2,lambda,.30,
    tenure_heaping_prob=.1,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.15,
    tenure_exact_anchor_retention_prob=.25,
    interview_month_mat=interview)
  nested <- log_emission_spell_g_monthly(g,s,0:2,lambda,.30,
    tenure_heaping_prob=.1,tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.15,
    tenure_exact_anchor_retention_prob=.25,
    tenure_local_revision_prob=0,
    interview_month_mat=interview)
  expect_equal(base$loglik,nested$loglik,tolerance=1e-12)
  expect_equal(base$tau_sum,nested$tau_sum,tolerance=1e-12)
})

test_that("local revision favors nearby non-exact changes", {
  lambda <- c(.20,.16,.12,.09,.06)
  s <- matrix(TRUE,1,2)
  interview <- matrix(c(3L,6L),1,2)
  nearby <- matrix(c(26L,31L)/12,1,2)
  distant <- matrix(c(26L,39L)/12,1,2)
  evaluate <- function(g,kappa) log_emission_spell_g_monthly(
    g,s,0:1,lambda,.999,tenure_heaping_prob=.1,
    tenure_year_revision_prob=.2,
    tenure_clean_anchor_revision_prob=.2,
    tenure_exact_anchor_retention_prob=.2,
    tenure_local_revision_prob=kappa,
    interview_month_mat=interview)$loglik
  expect_gt(evaluate(nearby,.7)-evaluate(nearby,0),
    evaluate(distant,.7)-evaluate(distant,0))
})

test_that("calendar pack and unpack preserves local revision", {
  p <- list(alpha=.55,theta0=.02,theta1=.04,pi=.03,eps=.20,eps_d=.12,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5))
  z <- .piecewise_calendar_revision_monthly_pack(p,.20,.10,.04,.08,
    rep(1/12,12),.25,.35)
  restored <- .piecewise_calendar_revision_monthly_unpack(z)
  expect_equal(restored$tenure_exact_anchor_retention_prob,.25,
    tolerance=1e-12)
  expect_equal(restored$tenure_local_revision_prob,.35,tolerance=1e-12)
})

test_that("shared reliability components move both contamination probabilities", {
  p <- list(eps=.25,eps_d=.15,duration_reliability_shift=1.2)
  reliable <- .duration_reliability_component_params(p,"reliable")
  unreliable <- .duration_reliability_component_params(p,"unreliable")
  expect_lt(reliable$eps,p$eps)
  expect_lt(reliable$eps_d,p$eps_d)
  expect_gt(unreliable$eps,p$eps)
  expect_gt(unreliable$eps_d,p$eps_d)
  expect_equal(reliable$duration_reliability_shift,0)
  expect_equal(unreliable$duration_reliability_shift,0)
})

test_that("zero shared-reliability shift exactly nests the local model", {
  p <- list(alpha=.55,pi=.03,eps=.25,eps_d=.15,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5),
    beta_g=0,beta_d=0,
    tenure_measurement_model="monthly",tenure_heaping_prob=.01,
    tenure_year_revision_prob=.10,
    tenure_clean_anchor_revision_prob=.12,
    tenure_exact_anchor_retention_prob=.20,
    tenure_local_revision_prob=.30)
  d <- simulate_eps_piecewise_job_change(80,p,seed=24024)
  base <- e_step_eps(d,p,suff_stats=FALSE)
  p$duration_reliability_shift <- 0
  nested <- e_step_eps(d,p,suff_stats=FALSE)
  expect_equal(base$loglik,nested$loglik,tolerance=1e-12)
  expect_equal(base$row_loglik,nested$row_loglik,tolerance=1e-12)
  expect_equal(base$gamma,nested$gamma,tolerance=1e-12)
})

test_that("shared-reliability likelihood equals its two-point mixture", {
  p <- list(alpha=.55,pi=.03,eps=.25,eps_d=.15,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5),
    beta_g=0,beta_d=0,
    tenure_measurement_model="monthly",tenure_heaping_prob=.01,
    tenure_year_revision_prob=.10,
    tenure_clean_anchor_revision_prob=.12,
    tenure_exact_anchor_retention_prob=.20,
    tenure_local_revision_prob=.30)
  d <- simulate_eps_piecewise_job_change(60,p,seed=24025)
  p$duration_reliability_shift <- .9
  mix <- e_step_eps(d,p,suff_stats=FALSE)
  reliable <- e_step_eps(d,.duration_reliability_component_params(p,
    "reliable"),suff_stats=FALSE)
  unreliable <- e_step_eps(d,.duration_reliability_component_params(p,
    "unreliable"),suff_stats=FALSE)
  maximum <- pmax(reliable$row_loglik,unreliable$row_loglik)
  expected <- maximum+log(.5*exp(reliable$row_loglik-maximum)+
    .5*exp(unreliable$row_loglik-maximum))
  expect_equal(mix$row_loglik,expected,tolerance=1e-12)
  expect_equal(rowSums(mix$gamma),rep(1,nrow(d)),tolerance=1e-12)
  expect_true(all(mix$duration_reliability_posterior>0 &
    mix$duration_reliability_posterior<1))
  expect_error(e_step_eps(d,p,suff_stats=TRUE),"suff_stats=FALSE")
})

test_that("calendar pack and unpack preserves shared reliability", {
  p <- list(alpha=.55,theta0=.02,theta1=.04,pi=.03,eps=.20,eps_d=.12,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5))
  z <- .piecewise_calendar_revision_monthly_pack(p,.20,.10,.04,.08,
    rep(1/12,12),.25,.35,1.1)
  restored <- .piecewise_calendar_revision_monthly_unpack(z)
  expect_equal(restored$duration_reliability_shift,1.1,tolerance=1e-12)
})

test_that("shared-reliability simulator records both latent classes", {
  p <- list(alpha=.55,pi=.03,eps=.25,eps_d=.15,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5),
    tenure_measurement_model="monthly",duration_reliability_shift=1)
  d <- simulate_eps_piecewise_job_change(500,p,seed=24026)
  expect_true(all(c("reliable","unreliable") %in%
    unique(d$duration_reliability_class)))
  expect_equal(mean(d$duration_reliability_class=="unreliable"),.5,
    tolerance=.08)
})

test_that("separate reliability shifts move the intended clocks", {
  p <- list(eps=.25,eps_d=.15,tenure_reliability_shift=1.2,
    timegap_reliability_shift=.4)
  reliable <- .duration_reliability_component_params(p,"reliable")
  unreliable <- .duration_reliability_component_params(p,"unreliable")
  expect_equal(qlogis(unreliable$eps)-qlogis(p$eps),1.2,
    tolerance=1e-12)
  expect_equal(qlogis(unreliable$eps_d)-qlogis(p$eps_d),.4,
    tolerance=1e-12)
  expect_equal(qlogis(p$eps)-qlogis(reliable$eps),1.2,
    tolerance=1e-12)
  expect_equal(qlogis(p$eps_d)-qlogis(reliable$eps_d),.4,
    tolerance=1e-12)
  expect_equal(reliable$tenure_reliability_shift,0)
  expect_equal(reliable$timegap_reliability_shift,0)
})

test_that("equal separate shifts reproduce common reliability", {
  p <- list(alpha=.55,pi=.03,eps=.25,eps_d=.15,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5),
    beta_g=0,beta_d=0,tenure_measurement_model="monthly",
    tenure_heaping_prob=.01,tenure_year_revision_prob=.10,
    tenure_clean_anchor_revision_prob=.12,
    tenure_exact_anchor_retention_prob=.20,
    tenure_local_revision_prob=.30,duration_reliability_shift=.9)
  d <- simulate_eps_piecewise_job_change(70,p,seed=24027)
  common <- e_step_eps(d,p,suff_stats=FALSE)
  p$duration_reliability_shift <- 0
  p$tenure_reliability_shift <- .9
  p$timegap_reliability_shift <- .9
  separate <- e_step_eps(d,p,suff_stats=FALSE)
  expect_equal(separate$loglik,common$loglik,tolerance=1e-12)
  expect_equal(separate$row_loglik,common$row_loglik,tolerance=1e-12)
  expect_equal(separate$duration_reliability_posterior,
    common$duration_reliability_posterior,tolerance=1e-12)
})

test_that("calendar pack and unpack preserves separate reliability shifts", {
  p <- list(alpha=.55,theta0=.02,theta1=.04,pi=.03,eps=.20,eps_d=.12,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5))
  z <- .piecewise_calendar_revision_monthly_pack(p,.20,.10,.04,.08,
    rep(1/12,12),.25,.35,NULL,1.1,.6)
  restored <- .piecewise_calendar_revision_monthly_unpack(z)
  expect_equal(restored$tenure_reliability_shift,1.1,tolerance=1e-12)
  expect_equal(restored$timegap_reliability_shift,.6,tolerance=1e-12)
  expect_equal(restored$duration_reliability_shift,0,tolerance=1e-12)
  expect_error(.piecewise_calendar_revision_monthly_pack(p,.20,.10,.04,
    .08,rep(1/12,12),.25,.35,.8,1.1,.6),"cannot both")
})

test_that("fixed-difference reduction and expansion are exact", {
  p <- list(alpha=.55,theta0=.02,theta1=.04,pi=.03,eps=.20,eps_d=.12,
    job_change_prob=.04,lambda_g=rep(.20,5),lambda_d=rep(.30,5))
  z <- .piecewise_calendar_revision_monthly_pack(p,.20,.10,.04,.08,
    rep(1/12,12),.25,.35,NULL,1.1,.6)
  difference <- unname(z["timegap_reliability_shift"]-
    z["tenure_reliability_shift"])
  reduced <- .piecewise_calendar_revision_monthly_difference_reduce(z)
  expanded <- .piecewise_calendar_revision_monthly_difference_expand(
    reduced,difference)[names(z)]
  expect_equal(expanded,z,tolerance=1e-14)
  equal <- .piecewise_calendar_revision_monthly_difference_expand(
    reduced,0)
  expect_equal(unname(equal["tenure_reliability_shift"]),
    unname(equal["timegap_reliability_shift"]),tolerance=1e-14)
  expect_error(.piecewise_calendar_revision_monthly_difference_expand(
    reduced,7),"[-6,6]",fixed=TRUE)
})
