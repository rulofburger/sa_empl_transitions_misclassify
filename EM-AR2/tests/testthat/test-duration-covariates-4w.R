test_that("duration AR(2) likelihood scores agree with numerical derivatives", {
  set.seed(19)
  n <- 24
  data <- list(y=matrix(rbinom(4*n,1,.5),n,4), weight=runif(n,.5,1.5),
    weight_sq=rep(1,n), n_original=n,
    X_entry=lapply(1:2,function(i)cbind(intercept=1,lag2=rbinom(n,1,.5),
      log_time_since_work=rnorm(n),never_worked=rbinom(n,1,.2),timegap_missing=0)),
    X_persistence=lapply(1:2,function(i)cbind(intercept=1,lag2=rbinom(n,1,.5),
      log_tenure=rnorm(n),tenure_missing=0)))
  z <- c(.1,-.2,.3,-1.6,.2,-.1,.15,.05,1.5,.2,.3,-.1,-2.5)
  analytic <- colSums(.ar2_duration_log_components(z,data,TRUE)$scores*data$weight)
  fn <- function(x) sum(data$weight*.ar2_duration_log_components(x,data)$loglik_i)
  numeric <- vapply(seq_along(z),function(j){
    step <- 1e-5; zp<-zm<-z; zp[j]<-zp[j]+step; zm[j]<-zm[j]-step
    (fn(zp)-fn(zm))/(2*step)
  },numeric(1))
  expect_equal(analytic,numeric,tolerance=2e-4)
})

test_that("duration AR(2) probabilities are valid", {
  n <- 10
  data <- list(X_entry=list(matrix(0,n,5)),X_persistence=list(matrix(0,n,4)))
  colnames(data$X_entry[[1]]) <- c("intercept","lag2","log_time_since_work","never_worked","timegap_missing")
  colnames(data$X_persistence[[1]]) <- c("intercept","lag2","log_tenure","tenure_missing")
  p <- .ar2_duration_unpack(rep(0,13),data)
  expect_equal(sum(p$alpha),1)
  expect_true(all(p$alpha>0))
  expect_true(p$pi>0 && p$pi<.5)
})
