# EM-tenure test runner
# Usage: source("EM-tenure/tests/testthat.R") or Rscript EM-tenure/tests/testthat.R
library(testthat)
test_dir(file.path(dirname(sys.frame(1)$ofile), "testthat"))
