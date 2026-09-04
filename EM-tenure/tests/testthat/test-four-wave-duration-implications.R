source(testthat::test_path("..","..","R","four_wave_duration_implications.R"))

test_that("Table 8 posterior decomposition partitions all latent paths", {
  h <- latent_histories_eps_4w()
  df <- as.data.frame(h); names(df) <- paste0("y",1:4)
  df$weight <- 1
  result <- posterior_duration_implications_4w(df,list(gamma=diag(16)),"test")
  expect_equal(result$classification_only,rep(0,3))
  expect_equal(result$true_reversal,rep(.5,3))
  expect_equal(result$true_persistent,rep(.5,3))
  expect_equal(unname(attr(result,"apparent_weight")),c(16,8,8))
})

test_that("opposite-direction true changes are not classification-only", {
  h <- latent_histories_eps_4w()
  df <- data.frame(y1=0,y2=1,y3=0,y4=1,weight=2)
  gamma <- matrix(as.numeric(apply(h,1,function(x) all(x==c(1,0,1,0)))),1)
  result <- posterior_duration_implications_4w(df,list(gamma=gamma),"test")
  expect_equal(result$classification_only,rep(0,3))
  expect_equal(result$true_reversal,rep(1,3))
  expect_equal(result$true_persistent,rep(0,3))
  gamma <- matrix(as.numeric(rowSums(h)==0),1)
  result <- posterior_duration_implications_4w(df,list(gamma=gamma),"test")
  expect_equal(result$classification_only,rep(1,3))
  expect_error(posterior_duration_implications_4w(df,list(gamma=gamma*.5),"test"))
})
