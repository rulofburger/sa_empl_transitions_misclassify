# ==============================================================================
# EM-AR2 tests: source all module files before running tests
# ==============================================================================
# This helper is automatically run by testthat before each test file.
# It sources all EM-AR2 R modules so functions are available in tests.
# ==============================================================================

# Determine project root (two levels up from tests/testthat/)
.test_root <- tryCatch(
  here::here(),
  error = function(e) {
    # Fallback: navigate up from current file location
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", "..", ".."))
  }
)

source(file.path(.test_root, "EM-AR2", "R", "source_all.R"))
rm(.test_root)
