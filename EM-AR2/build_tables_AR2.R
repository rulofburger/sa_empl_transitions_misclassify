# ==============================================================================
# EM-AR2: LaTeX tables with analytical delta-method standard errors
# ==============================================================================

if (!file.exists("EM-AR2/R/source_all.R"))
  stop("Source EM-AR2/build_tables_AR2.R from the project root")
source("EM-AR2/R/source_all.R")

results_dir <- "EM-AR2/output/results"
tables_dir <- "EM-AR2/output/tables"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
labels <- c("nopi", "sym", "asym")
headers <- c("No miscl.", "Symmetric $\\pi$", "Asymmetric $(\\pi_0,\\pi_1)$")

latest_fit <- function(label) {
  files <- list.files(results_dir,
    pattern = paste0("^em_ar2_", label, "_[0-9]{8}_[0-9]{6}\\.rds$"),
    full.names = TRUE)
  if (!length(files)) stop("No fit found for model ", label)
  readRDS(files[which.max(file.info(files)$mtime)])
}
fits <- setNames(lapply(labels, latest_fit), labels)
if (any(vapply(fits, function(x) is.null(x$analytical_inference), logical(1))))
  stop("Latest fits lack analytical inference; rerun EM-AR2/estimate_pipeline.R")

lookup <- function(label, quantity, field) {
  tab <- fits[[label]]$analytical_inference$summary
  value <- tab[tab$quantity == quantity, field]
  if (length(value)) value[[1L]] else NA_real_
}
stars <- function(est, se) {
  if (!is.finite(est) || !is.finite(se) || se <= 0) return("")
  z <- abs(est / se)
  if (z >= qnorm(.995)) "$^{***}$" else if (z >= qnorm(.975)) "$^{**}$" else
    if (z >= qnorm(.95)) "$^{*}$" else ""
}
cell <- function(label, quantity, scale = 1, digits = 4L) {
  est <- lookup(label, quantity, "estimate")
  se <- lookup(label, quantity, "se")
  if (!is.finite(est)) return(c("---", ""))
  c(paste0(formatC(scale * est, format = "f", digits = digits), stars(est, se)),
    paste0("(", formatC(scale * se, format = "f", digits = digits), ")"))
}

write_table <- function(rows, caption, label, path, percentage = FALSE) {
  scale <- if (percentage) 100 else 1
  digits <- if (percentage) 2L else 4L
  lines <- c("\\begin{table}[htbp]", "\\centering",
             paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
             "\\begin{tabular}{lccc}", "\\toprule",
             paste0(paste(c("", headers), collapse = " & "), " \\\\"), "\\midrule")
  for (quantity in names(rows)) {
    values <- lapply(labels, cell, quantity = quantity, scale = scale, digits = digits)
    lines <- c(lines,
      paste0(paste(c(rows[[quantity]], vapply(values, `[`, character(1), 1L)), collapse = " & "), " \\\\"),
      paste0(paste(c("", vapply(values, `[`, character(1), 2L)), collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\midrule",
    paste0(paste(c("Log-likelihood", vapply(fits, function(x) formatC(x$loglik, format="f", digits=1), character(1))), collapse=" & "), " \\\\"),
    paste0(paste(c("$N$", rep(format(fits[[1]]$n_obs, big.mark=","), 3)), collapse=" & "), " \\\\"),
    "\\bottomrule",
    "\\multicolumn{4}{l}{\\footnotesize Analytical survey-weighted sandwich standard errors in parentheses.} \\\\",
    "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
  message("Written: ", path)
}

write_table(
  c(theta0 = "$\\theta_0$", theta01 = "$\\theta_{01}$",
    theta1 = "$\\theta_1$", theta10 = "$\\theta_{10}$",
    pi = "$\\pi$", pi0 = "$\\pi_0$", pi1 = "$\\pi_1$"),
  "Four-wave AR(2) EM model: parameter estimates",
  "tab:ar2_params", file.path(tables_dir, "table_AR2_params.tex"))

write_table(
  c(p_00 = "$p_{00}$ (\\%)", p_10 = "$p_{10}$ (\\%)",
    p_01 = "$p_{01}$ (\\%)", p_11 = "$p_{11}$ (\\%)",
    employment_rate = "Steady-state employment rate (\\%)",
    pi = "$\\pi$ (\\%)", pi0 = "$\\pi_0$ (\\%)", pi1 = "$\\pi_1$ (\\%)"),
  "Four-wave AR(2) EM model: implied probabilities",
  "tab:ar2_implied", file.path(tables_dir, "table_AR2_implied.tex"),
  percentage = TRUE)
