# Transition-only no-misclassification benchmarks and duration diagnostics.
# Destination-wave durations are excluded from estimation and used only for
# posterior-predictive-style validation of clock-consistent spell histories.
if (!file.exists("data/raw/df_qlfs_A.rds")) stop("Run from the project root")
required_packages <- c("dplyr", "ggplot2", "splines", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ",
  paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

.gap_bounds <- c(0, .25, .5, .75, 1, 3, 5, Inf)
.duration_labels <- c("0-3 months", "3-6 months", "6-9 months",
  "9-12 months", "1-3 years", "3-5 years", "5+ years", "Never worked")
.gap_mid_years <- setNames(c(.125, .375, .625, .875, 2, 4, 7.5), 1:7)

duration_category <- function(x, never = FALSE) {
  # Treat an exact upper endpoint as belonging to the preceding bin. This is
  # especially important for a reported tenure of exactly three months.
  x_for_cut <- ifelse(is.finite(x) & x > 0, x - 1e-10, x)
  ans <- as.character(cut(x_for_cut, breaks = .gap_bounds, right = FALSE,
    labels = .duration_labels[1:7], include.lowest = TRUE))
  ans[never %in% TRUE] <- "Never worked"
  ans
}

prepare_panel <- function(waves) {
  raw <- readRDS("data/raw/df_qlfs_A.rds")
  wanted <- unique(c(paste0("employed", 1:waves), paste0("tenure", 1:waves),
    paste0("timegap", 1:waves), paste0("neverworked", 1:waves),
    paste0("age", 1:waves), paste0("educ", 1:waves),
    paste0("weight", 1:waves)))
  absent <- setdiff(wanted, names(raw))
  if (length(absent)) stop("Panel is missing: ", paste(absent, collapse = ", "))
  raw <- raw[, wanted, drop = FALSE]
  raw[] <- lapply(raw, function(x) as.numeric(unclass(x)))
  y_names <- paste0("employed", 1:waves)
  keep <- raw$age1 > 17 & raw$age1 < 56 & complete.cases(raw[y_names])
  raw <- raw[keep, , drop = FALSE]
  weight_mat <- as.matrix(raw[paste0("weight", 1:waves)])
  raw$weight <- exp(rowMeans(log(weight_mat)))
  keep_w <- is.finite(raw$weight) & raw$weight > 0
  raw <- raw[keep_w, , drop = FALSE]
  raw$weight <- nrow(raw) * raw$weight / sum(raw$weight)
  raw
}

make_transition_rows <- function(panel, model = c("AR(1)", "AR(2)")) {
  model <- match.arg(model)
  origin_waves <- if (model == "AR(1)") 1:2 else 2:3
  out <- lapply(origin_waves, function(t) {
    origin <- panel[[paste0("employed", t)]]
    destination <- panel[[paste0("employed", t + 1L)]]
    tenure_origin <- panel[[paste0("tenure", t)]] / 12
    tenure_destination <- panel[[paste0("tenure", t + 1L)]] / 12
    tenure_origin[tenure_origin <= 0] <- NA_real_
    tenure_destination[tenure_destination <= 0] <- NA_real_
    gap_origin_code <- panel[[paste0("timegap", t)]]
    gap_destination_code <- panel[[paste0("timegap", t + 1L)]]
    never_origin_raw <- panel[[paste0("neverworked", t)]]
    never_destination_raw <- panel[[paste0("neverworked", t + 1L)]]
    never_origin <- !is.na(never_origin_raw) & never_origin_raw == 1
    never_destination <- !is.na(never_destination_raw) & never_destination_raw == 1
    timegap_origin <- unname(.gap_mid_years[as.character(gap_origin_code)])
    # A never-worked respondent has no literal time since last job. Its age
    # since 16 is used only as an entry-risk covariate and is separately flagged.
    timegap_origin[never_origin] <- pmax(panel[[paste0("age", t)]][never_origin] - 16, 0)
    data.frame(id = seq_len(nrow(panel)), model = model, interval = t,
      lag2 = if (model == "AR(2)") panel[[paste0("employed", t - 1L)]] else 0,
      origin = origin, destination = destination,
      tenure_origin = tenure_origin, tenure_destination = tenure_destination,
      gap_origin_code = gap_origin_code, gap_destination_code = gap_destination_code,
      timegap_origin = timegap_origin,
      never_origin = as.integer(never_origin),
      never_destination = as.integer(never_destination), weight = panel$weight)
  })
  bind_rows(out)
}

fit_transition_model <- function(rows) {
  entry <- rows %>% filter(origin == 0, is.finite(timegap_origin))
  exit <- rows %>% filter(origin == 1, is.finite(tenure_origin))
  lag_term <- if (unique(rows$model) == "AR(2)") " + lag2" else ""
  entry_formula <- as.formula(paste0(
    "destination ~ splines::ns(log1p(timegap_origin), df=4) + never_origin + factor(interval)",
    lag_term))
  exit_formula <- as.formula(paste0(
    "I(1-destination) ~ splines::ns(log1p(tenure_origin), df=4) + factor(interval)",
    lag_term))
  entry_fit <- suppressWarnings(glm(entry_formula, data = entry,
    family = binomial(), weights = weight))
  exit_fit <- suppressWarnings(glm(exit_formula, data = exit,
    family = binomial(), weights = weight))
  rows$p_change <- NA_real_
  use_entry <- rows$origin == 0 & is.finite(rows$timegap_origin)
  use_exit <- rows$origin == 1 & is.finite(rows$tenure_origin)
  rows$p_change[use_entry] <- predict(entry_fit, rows[use_entry,], type = "response")
  rows$p_change[use_exit] <- predict(exit_fit, rows[use_exit,], type = "response")
  list(rows = rows, entry = entry_fit, exit = exit_fit)
}

weighted_distribution <- function(data, groups) {
  data %>% filter(is.finite(weight), weight > 0, !is.na(duration_bin)) %>%
    group_by(across(all_of(groups)), duration_bin, .drop = FALSE) %>%
    summarise(weight = sum(weight), .groups = "drop") %>%
    group_by(across(all_of(groups))) %>%
    mutate(share = weight / sum(weight)) %>% ungroup()
}

counterfactual_distributions <- function(fitted) {
  r <- fitted$rows %>% filter(is.finite(p_change))
  # Deterministic low-discrepancy draws approximate event time conditional on a
  # transition. They affect only the distribution within the first quarter.
  r$event_age <- ((r$id * .61803398875 + r$interval * .41421356237) %% 1) * .25
  observed <- bind_rows(
    r %>% filter(destination == 1, is.finite(tenure_destination)) %>%
      transmute(model, state = "Employment tenure", source = "Observed",
        duration_bin = duration_category(tenure_destination), weight),
    r %>% filter(destination == 0, gap_destination_code %in% 1:7 | never_destination == 1) %>%
      transmute(model, state = "Time since employment", source = "Observed",
        duration_bin = duration_category(unname(.gap_mid_years[as.character(gap_destination_code)]),
                                         never_destination == 1), weight))
  counterfactual <- bind_rows(
    r %>% filter(origin == 1, is.finite(tenure_origin)) %>%
      transmute(model, state = "Employment tenure", source = "No-error model",
        duration_bin = duration_category(tenure_origin + .25),
        weight = weight * (1 - p_change)),
    r %>% filter(origin == 0) %>%
      transmute(model, state = "Employment tenure", source = "No-error model",
        duration_bin = duration_category(event_age),
        weight = weight * p_change),
    r %>% filter(origin == 1) %>%
      transmute(model, state = "Time since employment", source = "No-error model",
        duration_bin = duration_category(event_age),
        weight = weight * p_change),
    r %>% filter(origin == 0, gap_origin_code %in% 1:7 | never_origin == 1) %>%
      transmute(model, state = "Time since employment", source = "No-error model",
        duration_bin = duration_category(timegap_origin + .25, never_origin == 1),
        weight = weight * (1 - p_change))
  )
  weighted_distribution(bind_rows(observed, counterfactual),
                        c("model", "state", "source"))
}

conditional_transition_distributions <- function(fitted) {
  r <- fitted$rows
  observed <- bind_rows(
    r %>% filter(origin == 0, destination == 1, is.finite(tenure_destination)) %>%
      transmute(model, transition = "Apparent entry: destination tenure",
        source = "Observed", duration_bin = duration_category(tenure_destination), weight),
    r %>% filter(origin == 1, destination == 0,
                 gap_destination_code %in% 1:7 | never_destination == 1) %>%
      transmute(model, transition = "Apparent exit: destination timegap",
        source = "Observed",
        duration_bin = duration_category(unname(.gap_mid_years[as.character(gap_destination_code)]),
                                         never_destination == 1), weight))
  counterfactual <- bind_rows(
    r %>% filter(origin == 0, destination == 1) %>%
      transmute(model, transition = "Apparent entry: destination tenure",
        source = "No-error clock", duration_bin = duration_category(.125), weight),
    r %>% filter(origin == 1, destination == 0) %>%
      transmute(model, transition = "Apparent exit: destination timegap",
        source = "No-error clock", duration_bin = duration_category(.125), weight))
  weighted_distribution(bind_rows(observed, counterfactual),
                        c("model", "transition", "source"))
}

gap_reachable <- function(prev, curr) {
  valid <- prev %in% 1:7 & curr %in% 1:7
  ans <- rep(FALSE, length(prev))
  for (i in which(valid)) {
    a0 <- .gap_bounds[prev[i]]; b0 <- .gap_bounds[prev[i] + 1L]
    a1 <- .gap_bounds[curr[i]]; b1 <- .gap_bounds[curr[i] + 1L]
    L <- max(a0, a1 - .25)
    U <- min(b0, if (is.infinite(b1)) Inf else b1 - .25)
    ans[i] <- L < U
  }
  ans
}

clock_diagnostics <- function(fitted) {
  r <- fitted$rows
  metric <- function(use, inconsistent, label) {
    ok <- use & !is.na(inconsistent)
    data.frame(model = unique(r$model), diagnostic = label,
      inconsistent_share = weighted.mean(inconsistent[ok], r$weight[ok]),
      risk_weight = sum(r$weight[ok]))
  }
  entry_use <- r$origin == 0 & r$destination == 1 & is.finite(r$tenure_destination)
  exit_use <- r$origin == 1 & r$destination == 0 &
    (r$gap_destination_code %in% 1:7 | r$never_destination == 1)
  ee_use <- r$origin == 1 & r$destination == 1 &
    is.finite(r$tenure_origin) & is.finite(r$tenure_destination)
  nn_use <- r$origin == 0 & r$destination == 0 &
    (r$gap_origin_code %in% 1:7 | r$never_origin == 1) &
    (r$gap_destination_code %in% 1:7 | r$never_destination == 1)
  nn_consistent <- (r$never_origin == 1 & r$never_destination == 1) |
    (r$never_origin == 0 & r$never_destination == 0 &
       gap_reachable(r$gap_origin_code, r$gap_destination_code))
  bind_rows(
    metric(entry_use, r$tenure_destination > .25,
      "Entry followed by tenure >3 months"),
    metric(exit_use, r$never_destination == 1 | r$gap_destination_code != 1,
      "Exit followed by timegap >3 months"),
    metric(ee_use, abs(r$tenure_destination-r$tenure_origin-.25) > 1/12,
      "Employment continuation with inconsistent clock"),
    metric(nn_use, !nn_consistent,
      "Nonemployment continuation with inconsistent clock"))
}

transition_summary <- function(fitted) {
  r <- fitted$rows %>% filter(is.finite(p_change))
  r %>% group_by(model, origin) %>% summarise(
    observed_change = weighted.mean(destination != origin, weight),
    fitted_change = weighted.mean(p_change, weight), n = n(), .groups = "drop") %>%
    mutate(equation = ifelse(origin == 0, "Entry", "Exit"))
}

weighted_cdf <- function(x, weight, grid) {
  use <- is.finite(x) & is.finite(weight) & weight > 0
  x <- x[use]; weight <- weight[use]
  ord <- order(x); x <- x[ord]; weight <- weight[ord]
  cw <- cumsum(weight); idx <- findInterval(grid, x)
  ifelse(idx == 0, 0, cw[pmax(idx,1)] / sum(weight))
}

tenure_cdf_comparison <- function(fitted, grid=seq(0,15,by=.05)) {
  r <- fitted$rows %>% filter(is.finite(p_change))
  r$event_age <- ((r$id * .61803398875 + r$interval * .41421356237) %% 1) * .25
  curve <- function(x,w,context,source) data.frame(model=unique(r$model),
    context=context,source=source,tenure_years=grid,
    cumulative_share=weighted_cdf(x,w,grid))
  observed_all <- r %>% filter(destination==1,is.finite(tenure_destination))
  simulated_stay <- r %>% filter(origin==1,is.finite(tenure_origin))
  simulated_enter <- r %>% filter(origin==0)
  observed_entry <- r %>% filter(origin==0,destination==1,
                                 is.finite(tenure_destination))
  bind_rows(
    curve(observed_all$tenure_destination,observed_all$weight,
          "All employed destination waves","Empirical"),
    curve(c(simulated_stay$tenure_origin+.25,simulated_enter$event_age),
          c(simulated_stay$weight*(1-simulated_stay$p_change),
            simulated_enter$weight*simulated_enter$p_change),
          "All employed destination waves","Simulated no-error"),
    curve(observed_entry$tenure_destination,observed_entry$weight,
          "Apparent entries only","Empirical"),
    curve(observed_entry$event_age,observed_entry$weight,
          "Apparent entries only","Simulated no-error"))
}

weighted_pdf <- function(x, weight, grid, bandwidth=.25) {
  use <- is.finite(x) & is.finite(weight) & weight > 0
  x <- x[use]; weight <- weight[use]
  # Reflect observations around zero to correct kernel-density boundary bias
  # for this nonnegative duration variable.
  reflected_x <- c(x, -x)
  reflected_weight <- c(weight, weight) / (2 * sum(weight))
  estimate <- density(reflected_x, weights=reflected_weight, bw=bandwidth,
    from=min(grid), to=max(grid), n=2048)
  2 * approx(estimate$x, estimate$y, xout=grid, rule=2)$y
}

tenure_pdf_comparison <- function(fitted, grid=seq(0,15,by=.025),
                                  bandwidth=.25) {
  r <- fitted$rows %>% filter(is.finite(p_change))
  r$event_age <- ((r$id * .61803398875 + r$interval * .41421356237) %% 1) * .25
  curve <- function(x,w,context,source) data.frame(model=unique(r$model),
    context=context,source=source,tenure_years=grid,
    density=weighted_pdf(x,w,grid,bandwidth))
  observed_all <- r %>% filter(destination==1,is.finite(tenure_destination))
  simulated_stay <- r %>% filter(origin==1,is.finite(tenure_origin))
  simulated_enter <- r %>% filter(origin==0)
  observed_entry <- r %>% filter(origin==0,destination==1,
                                 is.finite(tenure_destination))
  bind_rows(
    curve(observed_all$tenure_destination,observed_all$weight,
          "All employed destination waves","Empirical"),
    curve(c(simulated_stay$tenure_origin+.25,simulated_enter$event_age),
          c(simulated_stay$weight*(1-simulated_stay$p_change),
            simulated_enter$weight*simulated_enter$p_change),
          "All employed destination waves","Simulated no-error"),
    curve(observed_entry$tenure_destination,observed_entry$weight,
          "Apparent entries only","Empirical"),
    curve(observed_entry$event_age,observed_entry$weight,
          "Apparent entries only","Simulated no-error"))
}

panel3 <- prepare_panel(3)
panel4 <- prepare_panel(4)
fit_ar1 <- fit_transition_model(make_transition_rows(panel3, "AR(1)"))
fit_ar2 <- fit_transition_model(make_transition_rows(panel4, "AR(2)"))
fits <- list(fit_ar1, fit_ar2)

unconditional <- bind_rows(lapply(fits, counterfactual_distributions))
conditional <- bind_rows(lapply(fits, conditional_transition_distributions))
unconditional$duration_bin <- factor(unconditional$duration_bin,
                                     levels=.duration_labels)
conditional$duration_bin <- factor(conditional$duration_bin,
                                   levels=.duration_labels)
diagnostics <- bind_rows(lapply(fits, clock_diagnostics))
transitions <- bind_rows(lapply(fits, transition_summary))
tenure_cdf <- bind_rows(lapply(fits, tenure_cdf_comparison))
tenure_pdf <- bind_rows(lapply(fits, tenure_pdf_comparison))

stopifnot(all(is.finite(unconditional$share)), all(is.finite(conditional$share)),
  all(is.finite(diagnostics$inconsistent_share)),
  all(is.finite(tenure_cdf$cumulative_share)),
  all(is.finite(tenure_pdf$density)),
  max(abs(unconditional %>% group_by(model,state,source) %>%
    summarise(x=sum(share),.groups="drop") %>% pull(x) - 1)) < 1e-10,
  max(abs(conditional %>% group_by(model,transition,source) %>%
    summarise(x=sum(share),.groups="drop") %>% pull(x) - 1)) < 1e-10)

outdir <- "diagnostics/output/no_error_duration_counterfactual"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(unconditional, file.path(outdir, "unconditional_distributions.csv"), row.names=FALSE)
write.csv(conditional, file.path(outdir, "transition_distributions.csv"), row.names=FALSE)
write.csv(diagnostics, file.path(outdir, "clock_inconsistency.csv"), row.names=FALSE)
write.csv(transitions, file.path(outdir, "transition_fit.csv"), row.names=FALSE)
write.csv(tenure_cdf, file.path(outdir, "tenure_cdf_comparison.csv"), row.names=FALSE)
write.csv(tenure_pdf, file.path(outdir, "tenure_pdf_comparison.csv"), row.names=FALSE)

palette <- c("Observed" = "#172A3A", "No-error model" = "#F28E2B",
             "No-error clock" = "#F28E2B")
theme_diag <- theme_minimal(base_size=12) + theme(
  panel.grid.major.x = element_blank(), legend.position = "top",
  axis.text.x = element_text(angle=35,hjust=1),
  strip.text = element_text(face="bold"), plot.title.position="plot",
  plot.margin=margin(10,15,10,25))

p_unconditional <- ggplot(unconditional,
  aes(duration_bin, share, fill=source)) +
  geom_col(position=position_dodge(width=.78),width=.7) +
  facet_grid(model ~ state, scales="free_x") +
  scale_fill_manual(values=palette) + scale_y_continuous(labels=scales::percent) +
  labs(title="Observed durations versus transition-only no-error counterfactual",
    subtitle="Counterfactual clocks are propagated or reset using fitted flexible transition risks",
    x=NULL,y="Survey-weighted share",fill=NULL) + theme_diag
ggsave(file.path(outdir,"unconditional_duration_distributions.png"),
       p_unconditional,width=12,height=7.5,dpi=180)

p_conditional <- ggplot(conditional,
  aes(duration_bin,share,fill=source)) +
  geom_col(position=position_dodge(width=.78),width=.7) +
  facet_grid(model ~ transition,scales="free_x") +
  scale_fill_manual(values=palette) + scale_y_continuous(labels=scales::percent) +
  labs(title="Duration reports following apparent employment transitions",
    subtitle="With no status error, a new tenure or timegap cannot exceed the three-month interview interval",
    x=NULL,y="Survey-weighted share",fill=NULL) + theme_diag
ggsave(file.path(outdir,"transition_duration_distributions.png"),
       p_conditional,width=12,height=7.5,dpi=180)

p_diagnostics <- ggplot(diagnostics,
  aes(inconsistent_share,reorder(diagnostic,inconsistent_share),fill=model)) +
  geom_col(position=position_dodge(width=.75),width=.65) +
  geom_text(aes(label=scales::percent(inconsistent_share,accuracy=.1)),
    position=position_dodge(width=.75),hjust=-.08,size=3.5) +
  scale_x_continuous(labels=scales::percent,limits=c(0,max(diagnostics$inconsistent_share)*1.18)) +
  scale_fill_manual(values=c("AR(1)"="#4E79A7","AR(2)"="#59A14F")) +
  labs(title="Share of reported histories that violate a no-error spell clock",
    subtitle="One-month tolerance is allowed for tenure increments; timegap uses reported category intervals",
    x="Survey-weighted inconsistent share",y=NULL,fill=NULL) +
  theme_minimal(base_size=12) + theme(panel.grid.major.y=element_blank(),
    legend.position="top",plot.title.position="plot")
ggsave(file.path(outdir,"clock_inconsistency_shares.png"),
       p_diagnostics,width=11,height=6.5,dpi=180)

p_tenure_cdf <- ggplot(tenure_cdf,
  aes(tenure_years,cumulative_share,color=source,linetype=source)) +
  geom_line(linewidth=1.05) + facet_grid(model ~ context) +
  scale_color_manual(values=c("Empirical"="#172A3A",
    "Simulated no-error"="#F28E2B")) +
  scale_linetype_manual(values=c("Empirical"="solid",
    "Simulated no-error"="22")) +
  scale_x_continuous(breaks=c(0,.25,1,3,5,10,15),
    labels=c("","3m","1y","3y","5y","10y","15y")) +
  scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
  labs(title="Empirical and simulated employment-tenure distributions",
    subtitle="Flexible transition models use origin tenure only; destination tenure is held out",
    x="Tenure",y="Survey-weighted cumulative share",color=NULL,linetype=NULL) +
  theme_minimal(base_size=12) + theme(legend.position="top",
    panel.grid.minor=element_blank(),strip.text=element_text(face="bold"),
    plot.title.position="plot",plot.margin=margin(10,15,10,20))
ggsave(file.path(outdir,"tenure_empirical_vs_simulated_cdf.png"),
       p_tenure_cdf,width=12,height=7.5,dpi=180)

p_tenure_pdf <- ggplot(tenure_pdf,
  aes(tenure_years + 1/12,density,color=source,linetype=source)) +
  geom_line(linewidth=1.05) +
  facet_wrap(vars(model,context),ncol=2,scales="free_y") +
  scale_color_manual(values=c("Empirical"="#172A3A",
    "Simulated no-error"="#F28E2B")) +
  scale_linetype_manual(values=c("Empirical"="solid",
    "Simulated no-error"="22")) +
  scale_x_log10(breaks=c(0,.25,1,3,5,10,15) + 1/12,
    labels=c("0","3m","1y","3y","5y","10y","15y")) +
  scale_y_continuous(expand=expansion(mult=c(0,.05))) +
  labs(title="Empirical and simulated employment-tenure densities",
    subtitle="Survey-weighted densities; log(tenure + one month) axis and common three-month bandwidth",
    x="Tenure (shifted log scale)",y="Probability density",color=NULL,linetype=NULL) +
  theme_minimal(base_size=12) + theme(legend.position="top",
    panel.grid.minor=element_blank(),strip.text=element_text(face="bold"),
    plot.title.position="plot",plot.margin=margin(10,15,10,20))
ggsave(file.path(outdir,"tenure_empirical_vs_simulated_pdf.png"),
       p_tenure_pdf,width=12,height=7.5,dpi=180)
ggsave(file.path(outdir,"tenure_empirical_vs_simulated_pdf.pdf"),
       p_tenure_pdf,width=12,height=7.5,device=cairo_pdf)

timegap_plot_data <- bind_rows(
  unconditional %>% filter(state=="Time since employment") %>%
    transmute(model,context="All nonemployed destination waves",
      source=ifelse(source=="Observed","Empirical","Simulated no-error"),
      duration_bin,share),
  conditional %>% filter(transition=="Apparent exit: destination timegap") %>%
    transmute(model,context="Apparent exits only",
      source=ifelse(source=="Observed","Empirical","Simulated no-error"),
      duration_bin,share))
p_timegap <- ggplot(timegap_plot_data,
  aes(duration_bin,share,fill=source)) +
  geom_col(position=position_dodge(width=.78),width=.7) +
  facet_grid(model ~ context,scales="free_x") +
  scale_fill_manual(values=c("Empirical"="#172A3A",
    "Simulated no-error"="#F28E2B")) +
  scale_y_continuous(labels=scales::percent) +
  labs(title="Empirical and simulated time-since-employment distributions",
    subtitle="Reported categories are compared with clocks propagated or reset by the no-error model",
    x=NULL,y="Survey-weighted share",fill=NULL) + theme_diag
ggsave(file.path(outdir,"timegap_empirical_vs_simulated.png"),
       p_timegap,width=12,height=7.5,dpi=180)

cat("\nTransition-only model fit\n")
print(transitions,row.names=FALSE,digits=5)
cat("\nNo-error clock inconsistencies\n")
print(diagnostics,row.names=FALSE,digits=5)
cat("\nOutputs written to ",outdir,"\n",sep="")
