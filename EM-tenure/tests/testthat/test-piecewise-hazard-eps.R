test_that("piecewise cumulative hazard and inverse are exact", {
  h <- c(1,2,3,4,5)
  x <- c(.1,.25,1,3,5,8)
  expected <- c(.1,.25,1.75,7.75,15.75,30.75)
  H <- .duration_cumhaz(x,h,0)
  expect_equal(H,expected,tolerance=1e-12)
  expect_equal(.duration_inverse_cumhaz(H,h,0),x,tolerance=1e-12)
})

test_that("equal piecewise hazards nest the exponential model", {
  rate <- .31; h <- rep(rate,5)
  x <- c(.05,.25,.7,2,4,8)
  expect_equal(.duration_cumhaz(x,h,0),.duration_cumhaz(x,rate,0),
    tolerance=1e-12)
  expect_equal(.log_duration_density(x,h,0),
    .log_duration_density(x,rate,0),tolerance=1e-12)
  expect_equal(.duration_transition_probability(x,h,0),
    .duration_transition_probability(x,rate,0),tolerance=1e-12)
  expect_equal(.duration_category_transition_probability(1:7,h,0),
    .duration_category_transition_probability(1:7,rate,0),tolerance=1e-8)
})

test_that("positive final piecewise hazard gives a proper finite-mean duration", {
  h <- c(.5,.3,.2,.1,.04)
  mass <- integrate(function(x) exp(.log_duration_density(x,h,0)),
    0,Inf,subdivisions=500L)$value
  expect_equal(mass,1,tolerance=1e-6)
  expect_true(is.finite(duration_mean_years(h,0)))
  expect_true(.duration_marginal_transition_probability(h,0) > 0)
})

test_that("piecewise epsilon likelihood is finite with missing wrong-state clocks", {
  set.seed(501)
  n <- 60L; y <- replicate(3,sample(0:1,n,TRUE))
  df <- data.frame(y1=y[,1],y2=y[,2],y3=y[,3],
    tenure1=rexp(n,.3),tenure2=rexp(n,.3),tenure3=rexp(n,.3),
    timegap_cat1=sample(1:7,n,TRUE),timegap_cat2=sample(1:7,n,TRUE),
    timegap_cat3=sample(1:7,n,TRUE),weight=1)
  df <- prepare_eps_estimation_data(df)
  p <- list(alpha=.55,theta0=.05,theta1=.95,pi=.1,eps=.25,
    lambda_g=c(.35,.25,.18,.12,.08),beta_g=0,
    lambda_d=c(.4,.25,.15,.1,.06),beta_d=0)
  out <- e_step_eps(df,p,suff_stats=FALSE)
  expect_true(is.finite(out$loglik))
  expect_true(all(is.finite(out$gamma)))
  expect_equal(rowSums(out$gamma),rep(1,n),tolerance=1e-10)
})

test_that("zero timegap contamination exactly nests clock-consistent emissions", {
  hg <- c(.4,.3,.2,.15,.1)
  curr <- rep(1:7,each=7); prev <- rep(1:7,7)
  expect_equal(
    log_emission_transition_d_contaminated(curr,prev,hg,eps_d=0),
    log_emission_transition_d(curr,prev,hg),tolerance=0)
})

test_that("timegap contamination gives impossible clock changes finite mass", {
  hd <- c(.5,.3,.2,.15,.1)
  clean <- log_emission_transition_d(7L,1L,hd)
  mixed <- log_emission_transition_d_contaminated(7L,1L,hd,eps_d=.1)
  expect_true(is.infinite(clean) && clean < 0)
  expect_true(is.finite(mixed))
})

test_that("piecewise parameter transform retains timegap contamination", {
  p <- list(alpha=.5,pi=.08,eps=.22,eps_d=.28,
    lambda_g=c(.4,.3,.2,.15,.1),lambda_d=c(.5,.3,.2,.15,.1))
  z <- .piecewise_eps_pack(p,timegap_contamination=TRUE)
  out <- .piecewise_eps_unpack(z,timegap_contamination=TRUE)
  expect_equal(out$eps_d,p$eps_d,tolerance=1e-12)
  expect_equal(out$lambda_g,p$lambda_g,tolerance=1e-12)
  expect_equal(out$lambda_d,p$lambda_d,tolerance=1e-12)
})

test_that("piecewise likelihood with timegap contamination is finite", {
  set.seed(502)
  n <- 40L; y <- replicate(3,sample(0:1,n,TRUE))
  df <- data.frame(y1=y[,1],y2=y[,2],y3=y[,3],
    tenure1=rexp(n,.3),tenure2=rexp(n,.3),tenure3=rexp(n,.3),
    timegap_cat1=sample(1:7,n,TRUE),timegap_cat2=sample(1:7,n,TRUE),
    timegap_cat3=sample(1:7,n,TRUE),weight=1)
  df <- prepare_eps_estimation_data(df)
  p <- list(alpha=.55,theta0=.05,theta1=.95,pi=.1,eps=.25,eps_d=.1,
    lambda_g=c(.35,.25,.18,.12,.08),beta_g=0,
    lambda_d=c(.4,.25,.15,.1,.06),beta_d=0)
  out <- e_step_eps(df,p,suff_stats=FALSE)
  expect_true(is.finite(out$loglik))
  expect_equal(rowSums(out$gamma),rep(1,n),tolerance=1e-10)
})

test_that("timegap decomposition separates reset-compatible contradictions", {
  df <- data.frame(y1=c(0L,0L,0L),y2=c(0L,0L,0L),y3=c(1L,1L,1L),
    tenure1=NA_real_,tenure2=NA_real_,tenure3=NA_real_,
    timegap_cat1=c(2L,2L,2L),timegap_cat2=c(3L,1L,7L),
    timegap_cat3=NA_integer_,weight=1)
  p <- list(alpha=.5,theta0=.1,theta1=.9,pi=.05,eps=.2,eps_d=.25,
    lambda_g=rep(.2,5),beta_g=0,lambda_d=rep(.3,5),beta_d=0)
  fit <- list(params=p,gamma=matrix(1/8,nrow=3,ncol=8))
  out <- timegap_contamination_decomposition(df,fit)$summary
  expect_setequal(out$mechanism,
    c("Clock-feasible","Reset-compatible","Not reset-compatible"))
  expect_equal(sum(out$share_of_latent_continuations),1,tolerance=1e-12)
})

test_that("local timegap contamination is a proper conditional distribution", {
  hd <- c(.5,.3,.2,.15,.1)
  for (j in 1:7) {
    lp <- .log_timegap_local_contamination(1:7,rep(j,7),hd,decay=1)
    expect_equal(sum(exp(lp)),1,tolerance=1e-12)
  }
  expect_true(.log_timegap_local_contamination(7L,1L,hd,decay=1) <
    .log_timegap_local_contamination(3L,1L,hd,decay=1))
})

test_that("joint timegap contamination nests a clean three-wave clock", {
  hd <- c(.5,.3,.2,.15,.1)
  c1 <- c(1L,5L,6L); c2 <- c(2L,5L,7L); c3 <- c(3L,6L,7L)
  clean <- log_emission_interval_d(c1,hd)+
    log_emission_transition_d(c2,c1,hd)+
    log_emission_transition_d(c3,c2,hd)
  joint <- log_emission_timegap_triplet_joint(c1,c2,c3,hd,eps_d=1e-9)
  expect_equal(joint,clean,tolerance=1e-7)
})

test_that("local and joint timegap likelihoods are finite", {
  set.seed(503); n <- 30L
  df <- data.frame(y1=sample(0:1,n,TRUE),y2=sample(0:1,n,TRUE),
    y3=sample(0:1,n,TRUE),tenure1=rexp(n,.3),tenure2=rexp(n,.3),
    tenure3=rexp(n,.3),timegap_cat1=sample(1:7,n,TRUE),
    timegap_cat2=sample(1:7,n,TRUE),timegap_cat3=sample(1:7,n,TRUE),weight=1)
  df <- prepare_eps_estimation_data(df)
  base <- list(alpha=.5,theta0=.1,theta1=.9,pi=.08,eps=.2,eps_d=.2,
    lambda_g=c(.4,.3,.2,.15,.1),beta_g=0,
    lambda_d=c(.5,.3,.2,.15,.1),beta_d=0)
  for (model in c("local","joint_marginal")) {
    p <- base; p$timegap_contamination_model <- model
    out <- e_step_eps(df,p,suff_stats=FALSE)
    expect_true(is.finite(out$loglik))
    expect_equal(rowSums(out$gamma),rep(1,n),tolerance=1e-10)
  }
})
