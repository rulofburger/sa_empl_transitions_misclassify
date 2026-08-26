simulate_ar1_4w <- function(n,theta0=.08,theta1=.92,alpha=.5,pi=.02,seed=1) {
  set.seed(seed); h <- matrix(0L,n,4); h[,1] <- rbinom(n,1,alpha)
  for(tt in 2:4) h[,tt] <- rbinom(n,1,ifelse(h[,tt-1]==1,theta1,theta0))
  y <- abs(h-matrix(rbinom(4*n,1,pi),n,4))
  data.frame(y1=y[,1],y2=y[,2],y3=y[,3],y4=y[,4],weight=runif(n,.5,1.5))
}

test_that("cell probabilities are valid", {
  p <- list(theta0=.08,theta1=.92,alpha=.48,pi=.02)
  x <- ar1_4w_cell_probabilities(p,"symmetric")
  expect_length(x,16); expect_true(all(x>0)); expect_equal(sum(x),1,tolerance=1e-12)
})

test_that("cell collapse preserves weights and squared weights", {
  d <- simulate_ar1_4w(300,seed=2); c <- collapse_ar1_4w_cells(d)
  expect_equal(sum(c$count),nrow(d)); expect_equal(sum(c$weight),sum(d$weight),tolerance=1e-12)
  expect_equal(sum(c$weight_sq),sum(d$weight^2),tolerance=1e-12)
})

test_that("four-wave symmetric model is identified and recovers parameters", {
  d <- simulate_ar1_4w(10000,theta0=.07,theta1=.93,alpha=.46,pi=.025,seed=3)
  fit <- fit_ar1_4w_mle(d,"symmetric",FALSE,random_starts=4,verbose=0)
  expect_true(fit$converged); expect_true(fit$identified)
  expect_equal(fit$params$theta0,.07,tolerance=.015)
  expect_equal(fit$params$theta1,.93,tolerance=.015)
  expect_equal(fit$params$alpha,.46,tolerance=.02)
  expect_equal(fit$params$pi,.025,tolerance=.01)
  inf <- analytical_se_ar1_4w(d,fit)
  expect_true(all(is.finite(inf$summary$se))); expect_true(all(inf$summary$se>0))
})

test_that("free-alpha likelihood weakly dominates stationary likelihood", {
  d <- simulate_ar1_4w(3000,alpha=.35,seed=4)
  stat <- fit_ar1_4w_mle(d,"symmetric",TRUE,random_starts=2,verbose=0)
  free <- fit_ar1_4w_mle(d,"symmetric",FALSE,random_starts=2,verbose=0)
  expect_gte(free$loglik,stat$loglik-1e-5)
})
