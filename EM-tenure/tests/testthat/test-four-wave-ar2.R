ar2_test_parameters <- function() list(alpha=.49,pi=.04,eps=.2,eps_d=.15,
  job_change_prob=.02,lambda_g=c(.25,.3,.2,.15,.1),lambda_d=c(.5,.3,.2,.1,.08),
  beta_g=0,beta_d=0,tenure_measurement_model="monthly",timegap_contamination_model="joint_marginal",
  tenure_heaping_prob=.03,tenure_year_revision_prob=.1,tenure_clean_anchor_revision_prob=.15,
  tenure_exact_anchor_retention_prob=.2,tenure_local_revision_prob=.1,
  tenure_start_month_probs=(1:12)/78,tenure_reliability_shift=.6,timegap_reliability_shift=.8,
  timegap_clock_model="continuous_joint",ar2_entry=0,ar2_exit=0)

test_that("zero AR2 effects preserve the seeded sequential simulator", {
  p <- ar2_test_parameters()
  expect_identical(simulate_four_wave_ar2(200,p,seed=904,exact_risk=TRUE),
    simulate_eps_piecewise_job_change(200,p,seed=904,waves=4,exact_risk=TRUE))
})
test_that("AR2 likelihood nests AR1 and agrees with direct-prior reference", {
  withr::local_dir(testthat::test_path("..","..",".."))
  p <- ar2_test_parameters()
  d <- simulate_four_wave_ar2(100,p,seed=905,exact_risk=TRUE)
  d$weight <- seq_len(nrow(d))/mean(seq_len(nrow(d)))
  data <- prepare_four_wave_kernel_data(d)
  old <- evaluate_four_wave_monthly_fast(data,p,threads=1)
  new <- evaluate_four_wave_ar2(data,p,threads=1)
  expect_equal(new$loglik,old$loglik,tolerance=1e-10)
  expect_equal(new$gamma,old$gamma,tolerance=1e-12)
  expect_equal(duration_weighted_transition_rates_ar2(d,list(params=p,gamma=new$gamma))$entry_rate,
    duration_weighted_transition_rates_4w(d,list(params=p,gamma=old$gamma,
      integration_method="exact_piecewise"))$entry_rate,tolerance=1e-12)
  expect_equal(duration_weighted_transition_rates_ar2(d,list(params=p,gamma=new$gamma))$exit_rate,
    duration_weighted_transition_rates_4w(d,list(params=p,gamma=old$gamma,
      integration_method="exact_piecewise"))$exit_rate,tolerance=1e-12)
  expect_equal(new$loglik,sum(d$weight*new$row_loglik),tolerance=1e-10)
  for (b in list(c(.8,-.6),c(-1,1.2),c(3,3))) {
    p$ar2_entry<-b[1];p$ar2_exit<-b[2]
    fast <- evaluate_four_wave_ar2(data,p,threads=1)
    ref <- e_step_four_wave_ar2_reference(d,p)
    expect_equal(fast$row_loglik,ref$row_loglik,tolerance=1e-10)
    expect_equal(fast$gamma,ref$gamma,tolerance=1e-10)
    expect_equal(rowSums(fast$gamma),rep(1,nrow(d)),tolerance=1e-12)
  }
})
test_that("only recent transitions receive the appropriate AR2 odds shift", {
  p <- ar2_test_parameters();p$ar2_entry<-log(2);p$ar2_exit<-log(3)
  d <- simulate_four_wave_ar2(2,p,seed=906,exact_risk=TRUE)
  data <- prepare_four_wave_kernel_data(d); r <- .ar2_baseline_risks(data,p)
  h <- latent_histories_eps_4w();ratio <- .ar2_logratio(data,p)
  for (path in list(c(0,0,0,0),c(1,1,1,1),c(0,0,0,1),c(1,1,1,0))) {
    j <- which(apply(h,1,function(x) all(x==path)))
    expect_equal(ratio[,j],c(0,0))
  }
  j <- which(apply(h,1,function(x) all(x==c(1,0,1,1))))
  expected <- log(plogis(qlogis(r$entry[,2])+log(2))/r$entry[,2])+
    log((1-plogis(qlogis(r$exit[,3])+log(3)))/(1-r$exit[,3]))
  expect_equal(ratio[,j],expected,tolerance=1e-12)
  gamma <- matrix(0,nrow(d),16);gamma[,j] <- 1
  rate <- duration_weighted_transition_rates_ar2(d,list(params=p,gamma=gamma))[1,]
  expect_equal(rate$entry_rate,weighted.mean(.ar2_shift(r$entry[,2],log(2)),d$weight))
  expect_equal(rate$exit_rate,weighted.mean((r$exit[,1]+.ar2_shift(r$exit[,3],log(3)))/2,d$weight))
  expect_equal(.ar2_shift(.2,log(2))/(1-.ar2_shift(.2,log(2))),.5)
  expect_error(.ar2_effects(list(ar2_entry=NA)),"Invalid")
})
