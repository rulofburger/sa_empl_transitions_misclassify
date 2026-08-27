# ==============================================================================
# Matching-rule composition diagnostics and covariate-adjusted sensitivity
#
# This is a supplementary pipeline for Table 2. It deliberately does not write
# table_matching_implied.tex or any existing appendix-table input.
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2) # the shared ingestion script constructs diagnostic plots

source(here::here("EM-baseline", "R", "source_all.R"))
source(here::here("EM-baseline-ext", "R", "source_all.R"))

set.seed(20260826L)
results_dir <- here::here("EM-baseline", "output", "results")
tables_dir <- here::here("EM-baseline", "output", "tables")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

num <- function(x) as.numeric(unclass(x))

load_matching_panel <- function(panel) {
  old <- getOption("sa_empl_transitions.qlfs_3wave_panel")
  on.exit(options(sa_empl_transitions.qlfs_3wave_panel = old), add = TRUE)
  options(sa_empl_transitions.qlfs_3wave_panel = paste0("df_qlfs_", panel, ".rds"))
  ingest <- new.env(parent = globalenv())
  sys.source(here::here("scripts", "ingest_data_3waves_SA.R"), envir = ingest)
  d <- ingest$df_qlfs
  baseline <- complete.cases(d[c("y1", "y2", "y3", "weight")]) &
    is.finite(d$weight) & d$weight > 0
  covariate <- baseline & complete.cases(d[c(paste0("age", 1:3),
    paste0("educ", 1:3), "race1", "female1")])
  d <- d[covariate, , drop = FALSE]
  out <- data.frame(
    hhnr = d$hhnr, pnr = d$pnr,
    period1 = d$period1, period2 = d$period2, period3 = d$period3,
    y1 = as.integer(d$y1), y2 = as.integer(d$y2), y3 = as.integer(d$y3),
    weight = as.numeric(d$weight),
    age1 = num(d$age1), educ1 = num(d$educ1), race1 = num(d$race1),
    female1 = num(d$female1), stringsAsFactors = FALSE)
  attr(out, "baseline_n") <- sum(baseline)
  out
}

panels <- setNames(lapply(c("A", "B", "C"), load_matching_panel), c("A", "B", "C"))
cat("Complete-covariate samples:",
    paste(names(panels), vapply(panels, nrow, integer(1L)), collapse = "; "), "\n")

history_key <- function(d)
  paste(d$hhnr, d$pnr, d$period1, d$period2, d$period3, sep = "|")
key_B <- history_key(panels$B); key_C <- history_key(panels$C)
panels$A$group_B_not_C <- as.integer(history_key(panels$A) %in% key_B &
                                      !history_key(panels$A) %in% key_C)
panels$A$group_A_not_B <- as.integer(!history_key(panels$A) %in% key_B)
panels$A$group_C <- 1L - panels$A$group_B_not_C - panels$A$group_A_not_B
stopifnot(all(panels$A$group_C %in% 0:1),
          sum(panels$A$group_B_not_C * panels$A$group_A_not_B) == 0L)

# A common coding and scaling makes coefficients and predictive margins directly
# comparable across the three panels. Continuous variables use Panel A moments.
make_reference <- function(d) list(
  age = c(mean(d$age1), sd(d$age1)),
  age_sq = c(mean(d$age1^2), sd(d$age1^2)),
  educ = c(mean(d$educ1), sd(d$educ1)))
reference <- make_reference(panels$A)

make_set2 <- function(d, add_groups = FALSE) {
  X <- cbind(
    intercept = 1,
    age = (d$age1 - reference$age[1L]) / reference$age[2L],
    age_sq = (d$age1^2 - reference$age_sq[1L]) / reference$age_sq[2L],
    educ = (d$educ1 - reference$educ[1L]) / reference$educ[2L],
    race_2 = as.integer(d$race1 == 2), race_3 = as.integer(d$race1 == 3),
    race_4 = as.integer(d$race1 == 4), female = as.integer(d$female1))
  if (add_groups) X <- cbind(X, panel_B_not_C = d$group_B_not_C,
                             panel_A_not_B = d$group_A_not_B)
  attr(X, "entry_active") <- rep(TRUE, ncol(X))
  attr(X, "persistence_active") <- rep(TRUE, ncol(X))
  X
}
designs <- lapply(panels, make_set2)

collapse_cov <- function(df, X) {
  z <- .collapse_covariate_information_cells(df, X)
  z$n_original <- nrow(df)
  z
}

fit_cov_multistart <- function(df, X, model_type, nested = NULL, starts = 4L) {
  cells <- collapse_cov(df, X)
  seeds <- vector("list", starts)
  if (!is.null(nested)) {
    seeds[[1L]] <- nested$params; seeds[[1L]]$pi <- 0.01
  } else {
    seeds[[1L]] <- init_params_covariates(ncol(X), model_type)
    seeds[[1L]]$alpha <- sum(df$weight * df$y1) / sum(df$weight)
  }
  for (s in seq.int(2L, starts)) {
    p <- seeds[[1L]]
    p$beta0 <- p$beta0 + rnorm(length(p$beta0), 0, .20)
    p$beta1 <- p$beta1 + rnorm(length(p$beta1), 0, .20)
    p$alpha <- plogis(qlogis(p$alpha) + rnorm(1L, 0, .20))
    if (model_type == "symmetric")
      p$pi <- min(.45, plogis(qlogis(p$pi) + rnorm(1L, 0, .30)))
    seeds[[s]] <- p
  }
  candidates <- lapply(seeds, function(p) em_fit_covariates(
    cells$df, cells$X, model_type = model_type, stationary = FALSE,
    params0 = p, max_iter = 1000L, tol = 1e-8, verbose = 0L))
  fit <- candidates[[which.max(vapply(candidates, `[[`, numeric(1L), "loglik"))]]
  full <- e_step_covariates(df, X, fit$params, model_type, stationary = FALSE)
  fit$gamma <- full$gamma; fit$loglik <- full$loglik
  fit$starts <- data.frame(start = seq_along(candidates),
    loglik = vapply(candidates, `[[`, numeric(1L), "loglik"),
    converged = vapply(candidates, `[[`, logical(1L), "converged"))
  if (!fit$converged) stop("Covariate model failed to converge")
  fit
}

fits <- inference <- setNames(vector("list", 3L), names(panels))
for (P in names(panels)) {
  cat("\nEstimating adjusted models for panel", P, "\n")
  fit_path <- file.path(results_dir,
    paste0("matching_composition_fits_", tolower(P), ".rds"))
  inf_path <- file.path(results_dir,
    paste0("matching_composition_inference_", tolower(P), ".rds"))
  if (file.exists(fit_path) && file.exists(inf_path) &&
      !identical(Sys.getenv("MATCHING_COMPOSITION_REFIT"), "1")) {
    fits[[P]] <- readRDS(fit_path); inference[[P]] <- readRDS(inf_path)
    cat("  resumed saved point estimates and analytical inference\n")
    next
  }
  none <- fit_cov_multistart(panels[[P]], designs[[P]], "none")
  sym <- fit_cov_multistart(panels[[P]], designs[[P]], "symmetric", none)
  if (sym$loglik < none$loglik - 1e-4) stop("Nesting check failed in panel ", P)
  fits[[P]] <- list(none = none, symmetric = sym)
  inference[[P]] <- list(
    none = analytical_se_covariates(panels[[P]], designs[[P]], none, "none"),
    symmetric = analytical_se_covariates(panels[[P]], designs[[P]], sym, "symmetric"))
  saveRDS(fits[[P]], fit_path)
  saveRDS(inference[[P]], inf_path)
}

weighted_mean <- function(x, w) sum(w * x) / sum(w)
composition <- do.call(rbind, lapply(names(panels), function(P) {
  d <- panels[[P]]; w <- d$weight
  data.frame(panel = P, baseline_N = attr(d, "baseline_n"), N = nrow(d),
    age = weighted_mean(d$age1, w), education = weighted_mean(d$educ1, w),
    female = weighted_mean(d$female1, w),
    race_1 = weighted_mean(d$race1 == 1, w), race_2 = weighted_mean(d$race1 == 2, w),
    race_3 = weighted_mean(d$race1 == 3, w), race_4 = weighted_mean(d$race1 == 4, w),
    initial_employment = weighted_mean(d$y1, w))
}))
write.csv(composition, file.path(results_dir, "matching_composition_summary.csv"),
          row.names = FALSE)

predictive_margin <- function(params, X, w) c(
  entry = weighted_mean(pnorm(as.vector(X %*% params$beta0)), w),
  exit = weighted_mean(1 - pnorm(as.vector(X %*% params$beta1)), w),
  pi = unname(params$pi %||% NA_real_),
  initial_employment = unname(params$alpha))

margin_delta <- function(fit, inf, model_type, X, w) {
  eta <- .pack_covariate_params(fit$params, X, model_type)
  fn <- function(z) predictive_margin(
    .unpack_covariate_params(z, fit$params, X, model_type), X, w)
  q <- fn(eta); J <- matrix(NA_real_, length(q), length(eta),
                            dimnames = list(names(q), names(eta)))
  for (j in seq_along(eta)) {
    h <- 1e-5 * (1 + abs(eta[j])); ep <- em <- eta
    ep[j] <- ep[j] + h; em[j] <- em[j] - h
    J[, j] <- (fn(ep) - fn(em)) / (2 * h)
  }
  se <- sqrt(pmax(diag(J %*% inf$vcov_eta %*% t(J)), 0))
  data.frame(quantity = names(q), estimate = unname(q), se = unname(se))
}

margins <- list()
for (P in names(panels)) for (M in c("none", "symmetric")) {
  for (target_name in c("Own sample", "Panel C covariates")) {
    target <- if (target_name == "Own sample") P else "C"
    key <- paste(P, M, if (target_name == "Own sample") "own" else "C_standardized", sep = "_")
    margins[[key]] <- transform(margin_delta(fits[[P]][[M]], inference[[P]][[M]],
      M, designs[[target]], panels[[target]]$weight), panel = P, model = M,
      target = target_name)
  }
}
margins <- do.call(rbind, margins); rownames(margins) <- NULL
write.csv(margins, file.path(results_dir, "matching_standardized_margins.csv"),
          row.names = FALSE)

get_margin <- function(P, M, target, q)
  margins$estimate[margins$panel == P & margins$model == M &
    margins$target == target & margins$quantity == q]
decomposition <- do.call(rbind, lapply(c("A", "B"), function(P) {
  do.call(rbind, lapply(c("none", "symmetric"), function(M) {
    do.call(rbind, lapply(c("entry", "exit"), function(q) {
      own_P <- get_margin(P, M, "Own sample", q)
      std_P <- get_margin(P, M, "Panel C covariates", q)
      own_C <- get_margin("C", M, "Own sample", q)
      total <- own_P - own_C; comp <- own_P - std_P
      data.frame(comparison = paste0(P, "--C"), model = M, rate = q,
        total_gap = total, composition_component = comp,
        residual_component = std_P - own_C,
        share_explained = if (abs(total) < 1e-10) NA_real_ else comp / total)
    }))
  }))
}))
write.csv(decomposition, file.path(results_dir, "matching_composition_decomposition.csv"),
          row.names = FALSE)

# Pooled staged models: group indicators in neither equation, the transition
# equations only, the error equation only, and both equations.
make_stage <- function(transition_groups, error_groups) {
  X <- make_set2(panels$A, add_groups = transition_groups)
  Z0 <- cbind(error_intercept = rep(1, nrow(panels$A)))
  if (error_groups) Z0 <- cbind(Z0,
    panel_B_not_C = panels$A$group_B_not_C,
    panel_A_not_B = panels$A$group_A_not_B)
  Z <- replicate(3L, Z0, simplify = FALSE)
  collapse_covariate_reliability(panels$A, X, Z)
}

map_start <- function(data, previous = NULL) {
  Xnames <- colnames(.as_transition_design(data$X)$X12)
  delta_names <- colnames(data$Z[[1L]])
  p <- list(beta0 = setNames(rep(0, length(Xnames)), Xnames),
            beta1 = setNames(rep(0, length(Xnames)), Xnames),
            alpha = fits$A$symmetric$params$alpha,
            delta = setNames(rep(0, length(delta_names)), delta_names))
  for (block in c("beta0", "beta1")) {
    src <- if (is.null(previous)) fits$A$symmetric$params[[block]] else previous$params[[block]]
    common <- intersect(names(p[[block]]), names(src)); p[[block]][common] <- src[common]
  }
  p$delta["error_intercept"] <- qlogis(2 * fits$A$symmetric$params$pi)
  if (!is.null(previous)) {
    common <- intersect(names(p$delta), names(previous$params$delta))
    p$delta[common] <- previous$params$delta[common]
    p$alpha <- previous$params$alpha
  }
  p
}

fit_stage <- function(data, start, label) {
  eta <- pack_covariate_reliability(start, data$X)
  starts <- list(eta, eta + rnorm(length(eta), 0, .08),
                 eta + rnorm(length(eta), 0, .16))
  candidates <- lapply(starts, function(z)
    fit_covariate_reliability(data, z, maxit = 1500L, reltol = 1e-9))
  fit <- candidates[[which.max(vapply(candidates, `[[`, numeric(1L), "loglik"))]]
  if (!fit$converged) stop(label, " failed to converge")
  fit$candidate_table <- data.frame(start = seq_along(candidates),
    loglik = vapply(candidates, `[[`, numeric(1L), "loglik"),
    code = vapply(candidates, `[[`, integer(1L), "optimizer_code"))
  fit
}

stage_specs <- list(
  baseline = c(FALSE, FALSE), transition = c(TRUE, FALSE),
  error = c(FALSE, TRUE), both = c(TRUE, TRUE))
stage_data <- stage_fits <- stage_inf <- list()
stage_path <- file.path(results_dir, "matching_pooled_group_models.rds")
saved_stages <- if (file.exists(stage_path) &&
  !identical(Sys.getenv("MATCHING_COMPOSITION_REFIT"), "1")) readRDS(stage_path) else NULL
for (nm in names(stage_specs)) {
  cat("\nEstimating pooled stage:", nm, "\n")
  spec <- stage_specs[[nm]]
  stage_data[[nm]] <- make_stage(spec[1L], spec[2L])
  if (!is.null(saved_stages)) {
    stage_fits[[nm]] <- saved_stages$fits[[nm]]
    stage_inf[[nm]] <- saved_stages$inference[[nm]]
    cat("  resumed saved point estimates and analytical inference\n")
    next
  }
  previous <- if (nm == "baseline") NULL else stage_fits$baseline
  stage_fits[[nm]] <- fit_stage(stage_data[[nm]], map_start(stage_data[[nm]], previous), nm)
  stage_inf[[nm]] <- analytical_se_covariate_reliability(stage_data[[nm]], stage_fits[[nm]])
}
saveRDS(list(fits = stage_fits, inference = stage_inf, specs = stage_specs), stage_path)

scenario_margin <- function(fit, inf, data, group) {
  eta <- fit$eta
  fn <- function(z) {
    p <- unpack_covariate_reliability(z, data$X, colnames(data$Z[[1L]]))
    base <- make_set2(panels$A, add_groups = "panel_B_not_C" %in% names(p$beta0))
    if ("panel_B_not_C" %in% colnames(base)) {
      base[, "panel_B_not_C"] <- as.integer(group == "B_not_C")
      base[, "panel_A_not_B"] <- as.integer(group == "A_not_B")
    }
    z <- setNames(rep(0, length(p$delta)), names(p$delta)); z["error_intercept"] <- 1
    if ("panel_B_not_C" %in% names(z)) z["panel_B_not_C"] <- as.integer(group == "B_not_C")
    if ("panel_A_not_B" %in% names(z)) z["panel_A_not_B"] <- as.integer(group == "A_not_B")
    c(entry = weighted_mean(pnorm(as.vector(base %*% p$beta0)), panels$A$weight),
      exit = weighted_mean(1 - pnorm(as.vector(base %*% p$beta1)), panels$A$weight),
      pi = .5 * plogis(sum(z * p$delta)))
  }
  q <- fn(eta); J <- matrix(NA_real_, length(q), length(eta))
  for (j in seq_along(eta)) {
    h <- 1e-5 * (1 + abs(eta[j])); ep <- em <- eta
    ep[j] <- ep[j] + h; em[j] <- em[j] - h
    J[, j] <- (fn(ep) - fn(em)) / (2 * h)
  }
  data.frame(quantity = names(q), estimate = unname(q),
    se = sqrt(pmax(diag(J %*% inf$vcov_eta %*% t(J)), 0)))
}

stage_margins <- do.call(rbind, lapply(names(stage_fits), function(nm)
  do.call(rbind, lapply(c("C", "B_not_C", "A_not_B"), function(g)
    transform(scenario_margin(stage_fits[[nm]], stage_inf[[nm]], stage_data[[nm]], g),
              model = nm, group = g)))))
rownames(stage_margins) <- NULL
write.csv(stage_margins, file.path(results_dir, "matching_pooled_group_margins.csv"),
          row.names = FALSE)

# ---- Compact LaTeX writers --------------------------------------------------
fmt <- function(x, digits = 2L) ifelse(is.na(x), "---", sprintf(paste0("%.", digits, "f"), x))
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
write_table <- function(path, caption, label, headers, rows, note, align = NULL) {
  if (is.null(align)) align <- paste0("l", paste(rep("c", length(headers) - 1L), collapse = ""))
  lines <- c("\\begin{table}[htbp]", "\\centering", paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"), paste0("\\begin{tabular}{", align, "}"),
    "\\toprule", paste(headers, collapse = " & "), "\\\\", "\\midrule")
  for (r in rows) {
    lines <- c(lines, paste(r$est, collapse = " & "), "\\\\")
    if (!is.null(r$se) && any(nzchar(r$se))) lines <- c(lines, paste(r$se, collapse = " & "), "\\\\")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\begin{minipage}{0.98\\linewidth}",
    paste0("\\footnotesize \\textit{Note:} ", note), "\\end{minipage}", "\\end{table}")
  writeLines(lines, path)
}

comp_rows <- lapply(c("baseline_N", "N", "age", "education", "female", paste0("race_", 1:4),
                      "initial_employment"), function(q) {
  lab <- c(baseline_N = "$N$, Table 2 sample", N = "$N$, complete-covariate sample",
    age = "Age (years)", education = "Education (years)", female = "Female (\\%)",
    race_1 = "Black African (\\%)", race_2 = "Coloured (\\%)",
    race_3 = "Indian/Asian (\\%)", race_4 = "White (\\%)",
    initial_employment = "Initially employed (\\%)")[[q]]
  vals <- composition[[q]]
  if (q %in% c("baseline_N", "N")) vals <- vapply(vals, fmt_n, character(1L))
  else if (q %in% c("female", paste0("race_", 1:4), "initial_employment")) vals <- fmt(100 * vals)
  else vals <- fmt(vals)
  list(est = c(lab, vals), se = rep("", 4L))
})
write_table(file.path(tables_dir, "table_matching_composition.tex"),
  "Observable composition under alternative matching rules", "tab:matching_composition",
  c("", "Panel A", "Panel B", "Panel C"), comp_rows,
  "Weighted means use each panel's survey weights. The complete-covariate samples require nonmissing age and education in all three waves and nonmissing baseline race and sex; the first row reports the larger Table 2 estimation samples.")

model_order <- expand.grid(panel = c("A", "B", "C"), model = c("none", "symmetric"),
                           KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
model_order <- model_order[order(match(model_order$panel, c("A", "B", "C")),
                                 match(model_order$model, c("none", "symmetric"))), ]
margin_rows <- list()
for (q in c("entry", "exit", "pi")) for (target in c("Own sample", "Panel C covariates")) {
  est <- se <- numeric(nrow(model_order))
  for (j in seq_len(nrow(model_order))) {
    z <- margins[margins$panel == model_order$panel[j] & margins$model == model_order$model[j] &
      margins$target == target & margins$quantity == q, ]
    est[j] <- z$estimate; se[j] <- z$se
  }
  lab <- paste(if (q == "pi") "Misclassification rate" else paste0(tools::toTitleCase(q), " rate"),
               if (target == "Own sample") "-- own composition" else "-- Panel C composition", "(\\%)")
  margin_rows[[length(margin_rows) + 1L]] <- list(est = c(lab, fmt(100 * est)),
    se = c("", paste0("(", fmt(100 * se), ")")))
}
margin_rows <- c(margin_rows, list(
  list(est = c("Initial employment rate (\\%)", vapply(seq_len(nrow(model_order)), function(j)
    fmt(100 * fits[[model_order$panel[j]]][[model_order$model[j]]]$params$alpha), character(1L))), se = rep("", 7L)),
  list(est = c("$N$", vapply(model_order$panel, function(P) fmt_n(nrow(panels[[P]])), character(1L))), se = rep("", 7L))))
write_table(file.path(tables_dir, "table_matching_standardized.tex"),
  "Covariate-adjusted matching-rule sensitivity", "tab:matching_standardized",
  c("", "A: No error", "A: Symmetric", "B: No error", "B: Symmetric", "C: No error", "C: Symmetric"),
  margin_rows,
  "All models include Set 2 covariates (age, age squared, education, race, and sex), estimate initial employment freely, and do not impose stationarity. Panel C composition reports survey-weighted predictive margins evaluated on the same Panel C covariate distribution. Parentheses contain analytical sandwich/delta standard errors conditional on that target distribution.",
  align = "lcccccc")

dec_rows <- lapply(seq_len(nrow(decomposition)), function(i) {
  z <- decomposition[i, ]
  list(est = c(paste(z$comparison, if (z$model == "none") "No error" else "Symmetric",
                     tools::toTitleCase(z$rate), sep = ": "),
               fmt(100 * z$total_gap), fmt(100 * z$composition_component),
               fmt(100 * z$residual_component), fmt(100 * z$share_explained)),
       se = rep("", 5L))
})
write_table(file.path(tables_dir, "table_matching_decomposition.tex"),
  "Decomposition of matching-rule differences", "tab:matching_decomposition",
  c("", "Total gap (p.p.)", "Composition (p.p.)", "Residual (p.p.)", "Composition share (\\%)"),
  dec_rows,
  "The composition component is the own-sample predictive margin minus the same panel model evaluated on Panel C covariates. The residual is the Panel C-standardized prediction minus Panel C's own-sample prediction. These are descriptive decompositions; the final share can be unstable when the total gap is small.",
  align = "lcccc")

stage_headers <- c("", "(1) Baseline", "(2) Transitions", "(3) Error", "(4) Both")
stage_rows <- list()
for (g in c("C", "B_not_C", "A_not_B")) for (q in c("entry", "exit", "pi")) {
  vals <- ses <- numeric(4L)
  for (j in seq_along(stage_specs)) {
    z <- stage_margins[stage_margins$model == names(stage_specs)[j] &
      stage_margins$group == g & stage_margins$quantity == q, ]
    vals[j] <- z$estimate; ses[j] <- z$se
  }
  glab <- c(C = "Retained in C", B_not_C = "B but not C", A_not_B = "A but not B")[[g]]
  qlab <- c(entry = "entry", exit = "exit", pi = "misclassification")[[q]]
  stage_rows[[length(stage_rows) + 1L]] <- list(
    est = c(paste0(glab, ": ", qlab, " (\\%)"), fmt(100 * vals)),
    se = c("", paste0("(", fmt(100 * ses), ")")))
}
stage_rows <- c(stage_rows, list(
  list(est = c("Log-likelihood (millions)", fmt(vapply(stage_fits, `[[`, numeric(1L), "loglik") / 1e6, 3L)), se = rep("", 5L)),
  list(est = c("$N$", rep(fmt_n(nrow(panels$A)), 4L)), se = rep("", 5L)),
  list(est = c("Group indicators in transitions", "No", "Yes", "No", "Yes"), se = rep("", 5L)),
  list(est = c("Group indicators in error rate", "No", "No", "Yes", "Yes"), se = rep("", 5L))))
write_table(file.path(tables_dir, "table_matching_pooled_groups.tex"),
  "Pooled matching-group models", "tab:matching_pooled_groups", stage_headers, stage_rows,
  "All columns use Panel A, Set 2 controls, a symmetric error model, a free initial employment rate, and no stationarity restriction. Rates are predictive margins standardized to Panel A's covariate distribution while changing only the matching-group indicators. The error probability uses the bounded link $0.5\\,\\mathrm{logit}^{-1}(z'\\delta)$. Parentheses contain analytical sandwich/delta standard errors.",
  align = "lcccc")

cat("\nObservable composition:\n"); print(composition, row.names = FALSE, digits = 5)
cat("\nPanel-C standardized predictive margins:\n")
print(margins[margins$target == "Panel C covariates", ], row.names = FALSE, digits = 5)
cat("\nComposition decomposition:\n"); print(decomposition, row.names = FALSE, digits = 5)
cat("\nPooled matching-group predictive margins:\n"); print(stage_margins, row.names = FALSE, digits = 5)
cat("\nSupplementary Tables A2a--A2d written without changing Table 2 or Appendix Table A2.\n")
