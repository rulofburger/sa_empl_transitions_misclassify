make_covinc_4w_test_data <- function(n=80L) {
  set.seed(41)
  common <- cbind(age=rnorm(n),age_sq=rnorm(n),educ=rnorm(n),
    race_2=rbinom(n,1,.2),race_3=rbinom(n,1,.1),race_4=rbinom(n,1,.1),
    female=rbinom(n,1,.5))
  X <- lapply(1:3,function(tt)cbind(common,log_tenure=rnorm(n),
    log_time_since_work=rnorm(n),never_worked=rbinom(n,1,.1),
    tenure_missing=rbinom(n,1,.01),timegap_missing=rbinom(n,1,.02),
    permanent_contract=rbinom(n,1,.35),informal_sector=rbinom(n,1,.25)))
  Z <- lapply(1:4,function(tt)cbind(intercept=1,
    age_inconsistent=rbinom(n,1,.08),education_inconsistent=rbinom(n,1,.06)))
  list(y=matrix(rbinom(4*n,1,.5),n,4),weight=runif(n,.5,2),
    weight_sq=runif(n,.5,2)^2,X=X,Z=Z,
    entry_active=!colnames(X[[1]])%in%
      c("log_tenure","tenure_missing","permanent_contract","informal_sector"),
    persistence_active=!colnames(X[[1]])%in%
      c("log_time_since_work","never_worked","timegap_missing"),
    covariate_names=colnames(X[[1]]),n_original=n,n_cells=n)
}

testthat::test_that("four-wave type and state posteriors are normalized", {
  d <- make_covinc_4w_test_data(); p <- .initial_fmm_covinc_4w(d)
  e <- e_step_fmm_covinc_4w(d,p)
  testthat::expect_equal(rowSums(e$posterior_type),rep(1,nrow(d$y)),tolerance=1e-11)
  for(k in 1:2)for(tt in 1:4)
    testthat::expect_equal(rowSums(e$state[[k]][[tt]]),e$posterior_type[,k],tolerance=1e-11)
})

testthat::test_that("duration variables are active only in their relevant transition equation", {
  d <- make_covinc_4w_test_data(); p <- .initial_fmm_covinc_4w(d)
  testthat::expect_true("log_time_since_work" %in% names(p$beta0))
  testthat::expect_true("never_worked" %in% names(p$beta0))
  testthat::expect_false("log_tenure" %in% names(p$beta0))
  testthat::expect_true("log_tenure" %in% names(p$beta1))
  testthat::expect_false("log_time_since_work" %in% names(p$beta1))
})

testthat::test_that("analytical observed-likelihood score matches a direction derivative", {
  d <- make_covinc_4w_test_data(); p <- .initial_fmm_covinc_4w(d)
  o <- .direct_fmm_covinc_objective(d,p); z <- .pack_fmm_covinc_4w(p)
  set.seed(9); direction <- rnorm(length(z)); direction <- direction/sqrt(sum(direction^2))
  h <- 1e-5
  numerical <- (o$fn(z+h*direction)-o$fn(z-h*direction))/(2*h)
  analytical <- sum(o$gr(z)*direction)
  testthat::expect_equal(analytical,numerical,tolerance=2e-7)
})

testthat::test_that("type-specific misclassification score matches a direction derivative", {
  d <- make_covinc_4w_test_data(); p <- .initial_fmm_covinc_4w(d)
  p$delta <- setNames(c(p$delta[1],p$delta[1]+.2,p$delta[2:3]),
    c("intercept_A","intercept_B","age_inconsistent","education_inconsistent"))
  o <- .direct_fmm_covinc_objective(d,p); z <- .pack_fmm_covinc_4w(p)
  set.seed(19); direction <- rnorm(length(z)); direction <- direction/sqrt(sum(direction^2))
  h <- 1e-5
  numerical <- (o$fn(z+h*direction)-o$fn(z-h*direction))/(2*h)
  analytical <- sum(o$gr(z)*direction)
  testthat::expect_equal(analytical,numerical,tolerance=2e-7)
})

testthat::test_that("four-wave inconsistency attribution extends the three-wave rule", {
  d <- data.frame(age1=c(30,30),age2=c(35,31),age3=c(36,37),age4=c(37,38),
    educ1=c(5,5),educ2=c(5,5),educ3=c(5,5),educ4=c(5,5))
  z <- compute_inconsistencies_4w(d)
  testthat::expect_equal(unname(z[1,paste0("Y_age_",1:4)]),c(1,0,0,0))
  testthat::expect_equal(unname(z[2,paste0("Y_age_",1:4)]),c(0,0,1,0))
})
