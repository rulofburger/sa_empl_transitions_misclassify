
# Import
df_qlfs <- readRDS(
  here::here("data", "raw", "df_qlfs_A.rds")
  ) %>%
  filter(age1 > 17 & age1 < 56) %>%  
  filter(!is.na(employed1)) %>%
  filter(!is.na(employed2)) %>%
  filter(!is.na(employed3)) %>%
  select(contains(c("employed", "age", "educ", "race", "female", "weight", "tenure", "timegap", "formal", "period", "contracttype", "neverworked"))) %>%  # keep age, educ, empl status
  select(-contains(c("4"))) %>% 
  rename(
    y1 = employed1,
    y2 = employed2,
    y3 = employed3
  ) %>% 
  select(-contains(c("4"))) %>% 
  mutate(weight = (weight1*weight2*weight3)^0.333) %>% 
  as.data.frame

# Identify haven_labelled columns
haven_cols <- names(df_qlfs)[sapply(df_qlfs, haven::is.labelled)]

# a) Duplicate with _haven suffix (preserve originals)
df_qlfs[paste0(haven_cols, "_haven")] <- df_qlfs[haven_cols]

# b) Convert originals to numeric in-place
df_qlfs[haven_cols] <- lapply(df_qlfs[haven_cols], \(x) as.numeric(unclass(x)))

# Make Employment Binary
df_qlfs <- df_qlfs %>%
  mutate(
    tenure1 = if_else(y1 == 0, 0, tenure1),
    tenure2 = if_else(y2 == 0, 0, tenure2),
    tenure3 = if_else(y3 == 0, 0, tenure3),
    timegap1 = if_else(y1 == 1, 0, timegap1),
    timegap2 = if_else(y2 == 1, 0, timegap2),
    timegap3 = if_else(y3 == 1, 0, timegap3)
  ) %>% 
  mutate(
    tenure1 = if_else(tenure1 < 0, NA_real_, tenure1),
    tenure2 = if_else(tenure2 < 0, NA_real_, tenure2),
    tenure3 = if_else(tenure3 < 0, NA_real_, tenure3)
    )

# ---------------------------------------------------------------------------
# STEP 4a: Preserve raw timegap category codes BEFORE midpoint mapping
# ---------------------------------------------------------------------------
# Store the raw integer codes (1–7) as timegap_cat1–3.
#   Code 0  = "Never worked" (not applicable) → NA
#   Codes 1–7 = valid QLFS categories → integer 1L–7L
#   Codes 8, 99 = unspecified / missing → NA (dropped later)
#
# These category codes are used by the discrete interval-censored timegap
# model (discrete_timegap = TRUE in the EM). The continuous midpoints below
# are retained for descriptive use and backward compatibility.
df_qlfs <- df_qlfs %>%
  mutate(
    timegap_cat1 = case_when(
      timegap1 %in% 1:7 ~ as.integer(timegap1),
      TRUE               ~ NA_integer_
    ),
    timegap_cat2 = case_when(
      timegap2 %in% 1:7 ~ as.integer(timegap2),
      TRUE               ~ NA_integer_
    ),
    timegap_cat3 = case_when(
      timegap3 %in% 1:7 ~ as.integer(timegap3),
      TRUE               ~ NA_integer_
    )
  )


# Create timegap variables
df_qlfs <- df_qlfs %>%
  mutate(
    timegap1 = case_when(
      timegap1 == 0 ~ 0,
      timegap1 == 1 ~ 1.5,
      timegap1 == 2 ~ 4.5,
      timegap1 == 3 ~ 7.5,
      timegap1 == 4 ~ 10.5,
      timegap1 == 5 ~ 24,
      timegap1 == 6 ~ 48,
      timegap1 == 7 ~ 90,
      timegap1 == 8 ~ NA_real_,
      timegap1 == 99 ~ NA_real_
    )
  )%>%
  mutate(
    timegap2 = case_when(
      timegap2 == 0 ~ 0,
      timegap2 == 1 ~ 1.5,
      timegap2 == 2 ~ 4.5,
      timegap2 == 3 ~ 7.5,
      timegap2 == 4 ~ 10.5,
      timegap2 == 5 ~ 24,
      timegap2 == 6 ~ 48,
      timegap2 == 7 ~ 90,
      timegap2 == 8 ~ NA_real_,
      timegap2 == 99 ~ NA_real_
    )
  )%>%
  mutate(
    timegap3 = case_when(
      timegap3 == 0 ~ 0,
      timegap3 == 1 ~ 1.5,
      timegap3 == 2 ~ 4.5,
      timegap3 == 3 ~ 7.5,
      timegap3 == 4 ~ 10.5,
      timegap3 == 5 ~ 24,
      timegap3 == 6 ~ 48,
      timegap3 == 7 ~ 90,
      timegap3 == 8 ~ NA_real_,
      timegap3 == 99 ~ NA_real_
    )
  ) %>% 
  mutate(
    timegap1 = if_else(y1 == 0 & unclass(neverworked1) == 1, age1 - educ1 - 6, timegap1),
    timegap2 = if_else(y2 == 0 & unclass(neverworked2) == 1, age2 - educ2 - 6, timegap2),
    timegap3 = if_else(y3 == 0 & unclass(neverworked3) == 1, age3 - educ3 - 6, timegap3),
    timegap1 = if_else(timegap1 < 0, 0, timegap1),
    timegap2 = if_else(timegap2 < 0, 0, timegap2),
    timegap3 = if_else(timegap3 < 0, 0, timegap3)
  ) %>% 
  mutate(
    timegap1 = timegap1/12,
    timegap2 = timegap2/12,
    timegap3 = timegap3/12,
    tenure1 = tenure1/12,
    tenure2 = tenure2/12,
    tenure3 = tenure3/12
  )

df_qlfs <- df_qlfs %>% 
  mutate(
    race1 = as.numeric(unclass(race1)),
    race2 = as.numeric(unclass(race2)),
    race3 = as.numeric(unclass(race3)),
    female1 = as.numeric(unclass(female1)),
    female2 = as.numeric(unclass(female2)),
    female3 = as.numeric(unclass(female3)),
  )

# ---------------------------------------------------------------------------
# STEP 4b: Nearest-non-zero imputation for wrong-state durations
# ---------------------------------------------------------------------------
# BACKGROUND: The EM-tenure model must evaluate duration emissions for both
# the observed state and potential misclassified states. The ingest script
# sets tenure = 0 for nonemployed observations and timegap = 0 for employed
# observations. These zeros cause log_emg(0) = -Inf, which structurally
# rules out misclassification for those observations (DIAGNOSIS.md Issue 1).
#
# REMEDY: For each individual and each wave where the "wrong-state" duration
# is zero, impute a plausible non-zero value from another wave of the same
# individual.
#
# TIMEGAP imputation (for employed waves, timegap_cat = NA):
#   - Donor: the nearest wave t' where timegap_cat_{t'} is in 1:7 (valid
#     category, non-NA). Nearest means |t' - t| = 1 preferred over 2.
#   - IMPORTANT ASSUMPTION (Q6 from brainstorm): Only valid category codes
#     1–7 are eligible donors. Never-worked waves have timegap_cat = NA and
#     are EXCLUDED from the donor pool. This means an employed person at
#     wave t whose only nonemployed waves are never-worked receives the floor
#     category (1), not a never-worked duration. After the never-worked filter
#     in estimate_pipeline.R, never-worked individuals are removed entirely
#     from estimation, so this assumption does not affect EM estimates.
#   - Floor fallback: If all three waves are employed (no valid donor),
#     assign category 1 ([0, 0.25) years = [0, 3) months). This is the
#     smallest possible nonemployment duration — deliberately conservative
#     to avoid over-attributing misclassification.
#
# TENURE imputation (for nonemployed waves, tenure = 0):
#   - Donor: the nearest wave t' where tenure_{t'} > 0.
#   - Floor fallback: If all three waves are nonemployed (no valid donor),
#     assign .DURATION_FLOOR = 0.25/12 ≈ 0.021 years ≈ 1 week.
#
# NOTE: This imputation is performed on the FULL sample (before the
# never-worked filter in estimate_pipeline.R). Never-worked waves have
# timegap_cat = NA and are not eligible donors, so they do not contaminate
# the imputed values for other observations.
# ---------------------------------------------------------------------------

# Helper: nearest-non-zero donor index for a length-3 vector of values,
# where valid = non-NA and satisfies the validity condition.
# Returns the value from the nearest valid donor, or floor_val if none.
.impute_nearest <- function(v1, v2, v3, valid1, valid2, valid3, floor_val) {
  result <- c(
    if (!valid1) NA else v1,
    if (!valid2) NA else v2,
    if (!valid3) NA else v3
  )
  # For each position, find nearest valid donor
  out <- result
  for (t in 1:3) {
    if (!is.na(out[t])) next  # already valid
    # Prefer nearest: check t±1 before t±2
    donors_by_distance <- switch(t,
      `1` = c(2L, 3L),
      `2` = c(1L, 3L),
      `3` = c(2L, 1L)
    )
    found <- FALSE
    for (d in donors_by_distance) {
      if (!is.na(result[d])) {
        out[t] <- result[d]
        found <- TRUE
        break
      }
    }
    if (!found) out[t] <- floor_val
  }
  out
}

# Apply nearest-non-zero imputation row-wise using mapply
.floor_cat   <- 1L        # category 1 = [0, 0.25) years (floor for timegap)
.floor_dur   <- 0.25 / 12 # ≈ 0.021 years (floor for tenure)

# --- Timegap category imputation (for employed waves: y_t = 1) ---
# At employed waves, timegap_cat is currently NA (codes 0/8/99 map to NA,
# and employed waves get code 0 from the earlier if_else). We impute from
# the nearest wave with a valid category (1–7).
imputed_cats <- mapply(
  function(y1, y2, y3, c1, c2, c3) {
    .impute_nearest(
      v1 = c1, v2 = c2, v3 = c3,
      valid1 = !is.na(c1) & c1 >= 1L & c1 <= 7L,
      valid2 = !is.na(c2) & c2 >= 1L & c2 <= 7L,
      valid3 = !is.na(c3) & c3 >= 1L & c3 <= 7L,
      floor_val = .floor_cat
    )
  },
  df_qlfs$y1, df_qlfs$y2, df_qlfs$y3,
  df_qlfs$timegap_cat1, df_qlfs$timegap_cat2, df_qlfs$timegap_cat3,
  SIMPLIFY = TRUE  # returns 3 x N matrix
)
# Overwrite timegap_cat columns with imputed values (applied to ALL rows;
# for rows where the original was already valid, the value is unchanged)
df_qlfs$timegap_cat1 <- as.integer(imputed_cats[1L, ])
df_qlfs$timegap_cat2 <- as.integer(imputed_cats[2L, ])
df_qlfs$timegap_cat3 <- as.integer(imputed_cats[3L, ])

# --- Tenure imputation (for nonemployed waves: y_t = 0, tenure = 0) ---
# Nonemployed waves currently have tenure = 0 (set above by if_else).
# Impute from the nearest wave with tenure > 0.
imputed_ten <- mapply(
  function(y1, y2, y3, g1, g2, g3) {
    .impute_nearest(
      v1 = g1, v2 = g2, v3 = g3,
      valid1 = !is.na(g1) & g1 > 0,
      valid2 = !is.na(g2) & g2 > 0,
      valid3 = !is.na(g3) & g3 > 0,
      floor_val = .floor_dur
    )
  },
  df_qlfs$y1, df_qlfs$y2, df_qlfs$y3,
  df_qlfs$tenure1, df_qlfs$tenure2, df_qlfs$tenure3,
  SIMPLIFY = TRUE
)
df_qlfs$tenure1 <- imputed_ten[1L, ]
df_qlfs$tenure2 <- imputed_ten[2L, ]
df_qlfs$tenure3 <- imputed_ten[3L, ]

message(sprintf(
  "Nearest-non-zero imputation: timegap_cat and tenure imputed for wrong-state waves."
))

# ---------------------------------------------------------------------------
# ASSERTIONS: verify imputation produced valid values
# ---------------------------------------------------------------------------
stopifnot(all(df_qlfs$timegap_cat1 %in% 1:7, na.rm = TRUE))
stopifnot(all(df_qlfs$timegap_cat2 %in% 1:7, na.rm = TRUE))
stopifnot(all(df_qlfs$timegap_cat3 %in% 1:7, na.rm = TRUE))
# All employed obs should now have positive tenure (after imputation of nonemployed waves)
stopifnot(all(df_qlfs$tenure1[df_qlfs$y1 == 1] > 0, na.rm = TRUE))
stopifnot(all(df_qlfs$tenure2[df_qlfs$y2 == 1] > 0, na.rm = TRUE))
stopifnot(all(df_qlfs$tenure3[df_qlfs$y3 == 1] > 0, na.rm = TRUE))

# ---------------------------------------------------------------------------
# DROP ROWS WITH MISSING TENURE OR TIMEGAP (for EM-tenure estimation)
# ---------------------------------------------------------------------------
# The EM-tenure model requires complete duration data at all 3 waves.
# Tenure can be NA when an employed respondent has a negative raw value.
# Timegap can be NA from coded categories 8 (unspecified) or 99 (missing).
# We drop these rows here so all downstream estimation uses complete cases.
# This is a deliberate design choice: the EM emission densities need
# observed durations at every wave and cannot handle missing values.
n_before <- nrow(df_qlfs)
df_qlfs <- df_qlfs %>%
  filter(
    !is.na(tenure1) & !is.na(tenure2) & !is.na(tenure3) &
    !is.na(timegap1) & !is.na(timegap2) & !is.na(timegap3)
  )
n_after <- nrow(df_qlfs)
message(sprintf(
  "Dropped %d rows with NA tenure/timegap (%d -> %d, %.1f%% retained)",
  n_before - n_after, n_before, n_after, 100 * n_after / n_before
))

# model_tenure <- df_qlfs %>% 
#   filter(y1 == 1) %>% 
#   lm(tenure1 ~ age1 + educ1)

df_qlfs %>%
  filter(y1 == 1) %>%
  ggplot(aes(y = y2, x = log(tenure1 + 0.125))) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 5)
  )

df_qlfs %>%
  filter(y1 == 0) %>%
  ggplot(aes(y = y2, x = log(timegap1 + 0.125), group = neverworked1, color = as.factor(neverworked1))) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 5)
  )
