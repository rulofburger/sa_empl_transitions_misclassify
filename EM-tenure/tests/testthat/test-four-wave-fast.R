test_that("compiled input preparation retains the exact monthly grid", {
  df <- data.frame(weight = 1)
  for (j in 1:4) {
    df[[paste0("y", j)]] <- 1L
    df[[paste0("tenure", j)]] <- (j - 1) / 4
    df[[paste0("timegap_cat", j)]] <- NA_integer_
    df[[paste0("interview_month", j)]] <- j * 3L
  }
  prepared <- prepare_four_wave_kernel_data(df)
  expect_equal(as.vector(prepared$month), c(0L, 3L, 6L, 9L))
  expect_equal(as.vector(prepared$category), rep(0L, 4))
  df$tenure1 <- .01
  expect_error(prepare_four_wave_kernel_data(df), "monthly grid")
  df$tenure1 <- 0
  df$interview_month4 <- 13
  expect_error(prepare_four_wave_kernel_data(df), "Invalid interview")
})
