test_that("the sequential simulator preserves defaults and exposes four waves", {
  p <- list(alpha=.5, pi=.04, eps=.2, eps_d=.1, job_change_prob=.02,
    lambda_g=rep(.2,5),lambda_d=rep(.3,5),
    tenure_measurement_model="monthly")
  expect_identical(simulate_eps_piecewise_job_change(100,p,seed=32),
    simulate_eps_piecewise_job_change(100,p,seed=32,waves=3))
  d <- simulate_eps_piecewise_job_change(100,p,seed=32,waves=4,exact_risk=TRUE)
  expect_silent(validate_df_eps_4w(d))
  expect_true(all(is.na(d$tenure4[d$y4==0])))
  expect_true(all(is.na(d$timegap_cat4[d$y4==1])))
  expect_true(all(abs(12*d$tenure4-round(12*d$tenure4))<1e-10,na.rm=TRUE))
  expect_true(all(d$interview_month4==(d$interview_month3+2L)%%12L+1L))
  p$pi <- 0
  p$tenure_reliability_shift <- 1
  p$timegap_reliability_shift <- 2
  d <- simulate_eps_piecewise_job_change(200,p,seed=33,waves=4,exact_risk=TRUE)
  expect_equal(unname(as.matrix(d[paste0("y",1:4)])),
    unname(as.matrix(d[paste0("h",1:4)])))
  expect_setequal(unique(d$duration_reliability_class),c("reliable","unreliable"))
  expect_error(simulate_eps_piecewise_job_change(2,p,waves=5),"waves")
})
