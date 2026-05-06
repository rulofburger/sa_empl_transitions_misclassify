---
date: 2026-04-30
title: "NA validation must distinguish structural missingness from data errors"
category: "data-quality"
language: "R"
tags: [na-validation, qlfs, survey-data, employment-status, tenure, estep, input-validation]
root-cause: "Tenure NA check tested all three waves unconditionally; in QLFS data NA tenure at a non-employed wave is correct by design, not a data error."
severity: "P2"
---

# NA Validation Must Distinguish Structural Missingness from Data Errors

## Problem

The `e_step_eps()` input validation block contained:

```r
na_tenure <- is.na(df$tenure1) | is.na(df$tenure2) | is.na(df$tenure3)
if (any(na_tenure)) stop(...)
```

In QLFS panel data, workers who are **not employed** at a given wave legitimately
carry `NA` tenure (the variable is not applicable). This check therefore rejected
**every observation with a non-employed wave**, making the function unusable on real
QLFS data without a fragile upstream pre-imputation step.

The same pattern appeared (subtly) in downstream validation: `bad_cats` used
`na.rm = TRUE` in `all(..., na.rm = TRUE)` to check timegap categories, but the
`na_timegap` guard ran *after* `bad_cats`. If a timegap column had NAs *and*
out-of-range values, only the out-of-range error was reported; the NA was silently
swallowed by the first check.

## Root Cause

**Structural missingness** (NA because the variable does not apply at that
combination of values) was conflated with **data errors** (NA because of a coding
error, merge failure, or corrupted record).

The correct invariant for QLFS tenure data is:
- `y_t = 1` (employed) → `tenure_t` must be non-NA and positive
- `y_t = 0` (not employed) → `tenure_t` is naturally `NA` and must not be validated

## Solution

Replace the unconditional NA check with an employment-conditional check:

```r
# NA tenure is only an error at EMPLOYED waves. Non-employed waves
# naturally carry NA tenure in QLFS data.
na_tenure_emp <- (df$y1 == 1L & is.na(df$tenure1)) |
                 (df$y2 == 1L & is.na(df$tenure2)) |
                 (df$y3 == 1L & is.na(df$tenure3))
if (any(na_tenure_emp)) {
  stop(sprintf("e_step_eps: %d obs have NA tenure at an employed wave.",
               sum(na_tenure_emp)))
}
```

Also fix check ordering: move `na_timegap` guard **before** `bad_cats`, and remove
`na.rm = TRUE` from `bad_cats` so NA-containing columns trigger `bad_cats = TRUE`
rather than being silently passed:

```r
# Check NAs first so bad_cats cannot silently pass NA-containing columns.
na_timegap <- is.na(df$timegap_cat1) | is.na(df$timegap_cat2) |
              is.na(df$timegap_cat3)
if (any(na_timegap)) stop(sprintf("... %d obs have NA in timegap_cat.", sum(na_timegap)))

bad_cats <- !all(df$timegap_cat1 %in% 1:7) ||   # no na.rm = TRUE
            !all(df$timegap_cat2 %in% 1:7) ||
            !all(df$timegap_cat3 %in% 1:7)
```

## Prevention

When writing input validation for a model that operates on panel data:

1. **Map every variable to its applicability condition** before writing NA checks.
   Ask: "Under what rows is this variable expected to be defined?"
2. **Structural NA** = variable is inapplicable for this observation → allow it.
3. **Error NA** = variable should be defined but isn't → stop with a clear message.
4. **Ordering**: always check NA *before* range/value checks on the same column.
   Never rely on `na.rm = TRUE` inside a range check to handle upstream NAs —
   that creates silent pass-through.

Template for survey-data E-step validation:

```r
# 1. Required columns present
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

# 2. NA in always-defined columns (employment status, weights, timegap)
na_always <- is.na(df$y1) | is.na(df$y2) | is.na(df$y3) |
             is.na(df$timegap_cat1) | ...
if (any(na_always)) stop(sprintf("... %d obs have NA in required columns.", sum(na_always)))

# 3. Range checks on always-defined columns (now safe without na.rm = TRUE)
bad_cats <- !all(df$timegap_cat1 %in% 1:7) || ...
if (bad_cats) stop("timegap_cat values must be integers 1-7.")

# 4. Conditional NA check: tenure only at employed waves
na_tenure_emp <- (df$y1 == 1L & is.na(df$tenure1)) | ...
if (any(na_tenure_emp)) stop(sprintf("... %d employed obs have NA tenure.", sum(na_tenure_emp)))

# 5. Conditional range check: tenure > 0 only at employed waves
bad_tenure <- (df$y1 == 1L & df$tenure1 <= 0) | ...
if (any(bad_tenure)) stop(sprintf("... %d employed obs have non-positive tenure.", sum(bad_tenure)))
```

## Related

- [2026-03-10-em-gem-linked-ctmc.md](../data-quality/2026-03-10-em-gem-linked-ctmc.md)
- EM-tenure/R/estep_eps.R — fix applied here
