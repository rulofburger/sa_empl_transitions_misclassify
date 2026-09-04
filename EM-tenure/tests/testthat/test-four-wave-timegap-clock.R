test_that("joint continuous timegap clocks normalize for every observed subset", {
  h <- c(4,1,.5,.3,.2)
  for (mask in 0:15) {
    use <- which(as.logical(intToBits(mask)[1:4]))
    for (eps in c(0,.3,1)) {
      cats <- matrix(NA_integer_,if(length(use)) 7^length(use) else 1,4)
      if (length(use)) cats[,use] <- as.matrix(expand.grid(rep(list(1:7),length(use))))
      expect_equal(sum(exp(log_emission_timegap_clock_joint(cats,h,eps_d=eps))),
        1,tolerance=1e-12)
    }
  }
})

test_that("marginalizing any report leaves the same joint clock for the remainder", {
  h <- c(4,1,.5,.3,.2)
  cats <- as.matrix(expand.grid(rep(list(1:7),4)))
  mass <- exp(log_emission_timegap_clock_joint(cats,h,eps_d=.3))
  for (drop in 1:4) {
    key <- as.vector((cats[,-drop]-1)%*%7^(0:2))+1L
    marginal <- as.vector(rowsum(mass,key,reorder=TRUE))
    reduced <- matrix(NA_integer_,343,4)
    reduced[,-drop] <- as.matrix(expand.grid(rep(list(1:7),3)))
    expect_equal(marginal,exp(log_emission_timegap_clock_joint(reduced,h,eps_d=.3)),
      tolerance=1e-12)
  }
})

test_that("clean clock masses use the full interval intersection", {
  h <- c(4,1,.5,.3,.2)
  cats <- as.matrix(expand.grid(rep(list(1:7),4)))
  independent <- apply(cats,1,function(path) {
    intervals <- vapply(path,.timegap_interval,numeric(2))
    lower <- max(intervals[1,]-(0:3)/4)
    upper <- min(intervals[2,]-(0:3)/4)
    if (lower>=upper) 0 else exp(.log_duration_interval_prob(lower,upper,h))
  })
  expect_equal(exp(log_emission_timegap_clock_joint(cats,h,eps_d=0)),
    independent,tolerance=1e-12)
  expect_equal(log_emission_timegap_clock_joint(matrix(c(NA,NA,1),1),h,eps_d=0),-Inf)
  expect_true(is.finite(log_emission_timegap_clock_joint(
    matrix(c(NA,NA,1),1),h,eps_d=.2)))
  expect_error(log_emission_timegap_clock_joint(matrix(8,1),h,eps_d=.2),"categories")
})
