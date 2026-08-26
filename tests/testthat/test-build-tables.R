# Created: 2026-05-07
# Tests for helper functions defined in build_tables.R.
#
# These helpers are pure functions (no file I/O, no estimation) so they can be
# tested by sourcing a minimal preamble that defines the constants they depend on.
#
# Run from project root:
#   testthat::test_file("tests/testthat/test-build-tables.R")

library(testthat)

# ---------------------------------------------------------------------------
# Minimal preamble: define helpers in isolation (no full script sourcing).
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

.CRIT_p01 <- qnorm(0.995)
.CRIT_p05 <- qnorm(0.975)
.CRIT_p10 <- qnorm(0.95)

.stars <- function(est, se) {
  if (is.na(est) || is.na(se) || se <= 0) return("")
  z <- abs(est / se)
  if (z > .CRIT_p01) return("$^{**}$")
  if (z > .CRIT_p05) return("$^{*}$")
  if (z > .CRIT_p10) return("$^{.}$")
  return("")
}

.fmt_estimate <- function(est, se, factor = 1, digits = 4L) {
  if (is.na(est)) return(c("---", "(---)"))
  star    <- .stars(est, se)
  est_str <- paste0(formatC(est * factor, format = "f", digits = digits), star)
  se_str  <- if (!is.na(se)) {
    sprintf("(%s)", formatC(se * factor, format = "f", digits = digits))
  } else {
    "(---)"
  }
  c(est_str, se_str)
}

.fmt_param <- function(est, se, digits = 4L) .fmt_estimate(est, se, factor = 1,   digits = digits)
.fmt_pct   <- function(est, se, digits = 2L) .fmt_estimate(est, se, factor = 100, digits = digits)
.fmt_pp    <- function(est, se, digits = 2L) .fmt_pct(est, se, digits)

.fmt_ll <- function(ll) {
  if (is.na(ll) || is.null(ll)) return("---")
  sprintf("%.1fM", ll / 1e6)
}

.fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

.detect_B <- function(boot_dir) {
  fls <- sort(list.files(boot_dir, pattern = "boot_.*_B[0-9]+\\.rds$"))
  if (length(fls) == 0L) return(NULL)
  as.integer(sub(".*_B([0-9]+)\\.rds$", "\\1", fls[length(fls)]))
}

.get_se <- function(summary_or_map, quantity) {
  if (is.null(summary_or_map)) return(NA_real_)
  if (is.data.frame(summary_or_map)) {
    idx <- which(summary_or_map$quantity == quantity)
    if (length(idx) == 0L) NA_real_ else summary_or_map$se[idx[1L]]
  } else {
    v <- summary_or_map[[quantity]]
    if (is.null(v)) NA_real_ else v
  }
}

.se_map <- function(boot_obj) {
  if (is.null(boot_obj) || is.null(boot_obj$summary)) return(NULL)
  setNames(boot_obj$summary$se, boot_obj$summary$quantity)
}

# ---------------------------------------------------------------------------
# Tests: .CRIT_p10 correctness
# ---------------------------------------------------------------------------

test_that(".CRIT_p10 is qnorm(0.95), not qnorm(0.945)", {
  expect_equal(.CRIT_p10, qnorm(0.95), tolerance = 1e-10)
  expect_gt(.CRIT_p10, 1.644)
  expect_lt(.CRIT_p10, 1.645 + 0.001)
})

# ---------------------------------------------------------------------------
# Tests: .stars()
# ---------------------------------------------------------------------------

test_that(".stars returns empty string for NA est or se", {
  expect_equal(.stars(NA_real_, 0.1), "")
  expect_equal(.stars(1.0, NA_real_), "")
})

test_that(".stars returns empty string when se <= 0", {
  expect_equal(.stars(1.0, 0), "")
  expect_equal(.stars(1.0, -0.1), "")
})

test_that(".stars returns ** just above p<0.01 threshold", {
  # z just above 2.576
  z_crit <- .CRIT_p01 + 1e-8
  expect_equal(.stars(z_crit, 1), "$^{**}$")
})

test_that(".stars returns * just above p<0.05 threshold (but below p<0.01)", {
  z_crit <- .CRIT_p05 + 1e-8
  expect_equal(.stars(z_crit, 1), "$^{*}$")
})

test_that(".stars returns . just above p<0.10 threshold (but below p<0.05)", {
  z_crit <- .CRIT_p10 + 1e-8
  expect_equal(.stars(z_crit, 1), "$^{.}$")
})

test_that(".stars returns empty just below p<0.10 threshold", {
  z_crit <- .CRIT_p10 - 1e-8
  expect_equal(.stars(z_crit, 1), "")
})

test_that(".stars boundary: exactly at .CRIT_p01 gets ** (strictly greater required)", {
  # z == .CRIT_p01 exactly — NOT strictly greater → should NOT be **
  expect_equal(.stars(.CRIT_p01, 1), "$^{*}$")
})

# ---------------------------------------------------------------------------
# Tests: .fmt_estimate()
# ---------------------------------------------------------------------------

test_that(".fmt_estimate returns --- for NA estimate", {
  out <- .fmt_estimate(NA_real_, 0.05)
  expect_equal(out[1], "---")
  expect_equal(out[2], "(---)")
})

test_that(".fmt_estimate returns (---) when se is NA", {
  out <- .fmt_estimate(0.5, NA_real_)
  expect_match(out[2], "---")
})

test_that(".fmt_estimate scales by factor=100 (percentage)", {
  out <- .fmt_estimate(0.12345, 0.01, factor = 100, digits = 2L)
  expect_match(out[1], "12.35")
  expect_match(out[2], "1.00")
})

test_that(".fmt_estimate respects digits argument", {
  out <- .fmt_estimate(0.5, 0.05, factor = 1, digits = 2L)
  expect_match(out[1], "^0\\.50")
  out4 <- .fmt_estimate(0.5, 0.05, factor = 1, digits = 4L)
  expect_match(out4[1], "^0\\.5000")
})

test_that(".fmt_param formats to 4 d.p. at raw scale", {
  out <- .fmt_param(0.1234567, 0.01)
  expect_match(out[1], "^0\\.1235")
})

test_that(".fmt_pct formats to 2 d.p. at percent scale", {
  out <- .fmt_pct(0.1234, 0.01)
  expect_match(out[1], "^12\\.34")
})

# ---------------------------------------------------------------------------
# Tests: .fmt_ll()
# ---------------------------------------------------------------------------

test_that(".fmt_ll formats in millions", {
  expect_equal(.fmt_ll(-357200000), "-357.2M")
})

test_that(".fmt_ll returns --- for NA", {
  expect_equal(.fmt_ll(NA_real_), "---")
})

test_that(".fmt_ll returns --- for NULL", {
  expect_equal(.fmt_ll(NULL), "---")
})

# ---------------------------------------------------------------------------
# Tests: .fmt_n()
# ---------------------------------------------------------------------------

test_that(".fmt_n inserts comma separator", {
  expect_equal(.fmt_n(12345), "12,345")
})

test_that(".fmt_n handles small numbers without comma", {
  expect_equal(.fmt_n(999), "999")
})

# ---------------------------------------------------------------------------
# Tests: .detect_B()
# ---------------------------------------------------------------------------

test_that(".detect_B returns NULL on empty directory", {
  tmp <- withr::local_tempdir()
  expect_null(.detect_B(tmp))
})

test_that(".detect_B parses B from single file", {
  tmp <- withr::local_tempdir()
  file.create(file.path(tmp, "boot_none_stat_B200.rds"))
  expect_equal(.detect_B(tmp), 200L)
})

test_that(".detect_B returns largest B when multiple files present", {
  tmp <- withr::local_tempdir()
  file.create(file.path(tmp, "boot_none_stat_B100.rds"))
  file.create(file.path(tmp, "boot_none_stat_B200.rds"))
  expect_equal(.detect_B(tmp), 200L)
})

test_that(".detect_B uses sorted order (deterministic, not OS-order)", {
  tmp <- withr::local_tempdir()
  file.create(file.path(tmp, "boot_a_B50.rds"))
  file.create(file.path(tmp, "boot_b_B200.rds"))
  result1 <- .detect_B(tmp)
  result2 <- .detect_B(tmp)
  expect_equal(result1, result2)
  expect_equal(result1, 200L)
})

# ---------------------------------------------------------------------------
# Tests: .get_se()
# ---------------------------------------------------------------------------

test_that(".get_se returns NA_real_ for NULL input", {
  expect_identical(.get_se(NULL, "theta0"), NA_real_)
})

test_that(".get_se looks up from data.frame path", {
  df <- data.frame(quantity = c("theta0", "theta1"), se = c(0.01, 0.02),
                   stringsAsFactors = FALSE)
  expect_equal(.get_se(df, "theta0"), 0.01)
  expect_equal(.get_se(df, "theta1"), 0.02)
})

test_that(".get_se returns NA_real_ for missing quantity in data.frame", {
  df <- data.frame(quantity = "theta0", se = 0.01, stringsAsFactors = FALSE)
  expect_identical(.get_se(df, "theta1"), NA_real_)
})

test_that(".get_se looks up from named vector path (O(1) map)", {
  m <- c(theta0 = 0.01, theta1 = 0.02)
  expect_equal(.get_se(m, "theta0"), 0.01)
})

test_that(".get_se returns NA_real_ for missing key in named vector", {
  m <- c(theta0 = 0.01)
  expect_identical(.get_se(m, "theta1"), NA_real_)
})

# ---------------------------------------------------------------------------
# Tests: .se_map()
# ---------------------------------------------------------------------------

test_that(".se_map returns NULL for NULL input", {
  expect_null(.se_map(NULL))
})

test_that(".se_map returns NULL when $summary is NULL", {
  expect_null(.se_map(list(summary = NULL)))
})

test_that(".se_map converts data.frame to named vector correctly", {
  obj <- list(summary = data.frame(quantity = c("a", "b"), se = c(0.1, 0.2),
                                   stringsAsFactors = FALSE))
  m <- .se_map(obj)
  expect_equal(m[["a"]], 0.1)
  expect_equal(m[["b"]], 0.2)
  expect_true(is.numeric(m))
})

test_that(".se_map handles empty summary data.frame", {
  obj <- list(summary = data.frame(quantity = character(0), se = numeric(0)))
  m <- .se_map(obj)
  expect_true(is.numeric(m))
  expect_equal(length(m), 0L)
})
