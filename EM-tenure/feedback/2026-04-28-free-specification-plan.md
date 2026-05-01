# 2026-04-28: Free Specification — What's Changing and Why

## The Problem (recap)

The EM-tenure model uses the CTMC link to force a deterministic relationship
between the transition probabilities (θ₁, θ₀) and the exponential spell-length
rates (λ_g, λ_d):

    λ_g = −log(1 − θ₁) / 3

This means θ₁ must simultaneously:
1. Match the **transition counts** (which want θ₁ ≈ 0.97)
2. Match the **tenure distribution** (mean 6.5 years → wants λ_g very small → wants θ₁ very small)

These two signals are in irreconcilable conflict. The EM "compromises" at
θ₁ ≈ 0.50, which satisfies neither. The simple binary model (without durations)
gets θ₁ ≈ 0.97, θ₀ ≈ 0.03, π ≈ 3% — all sensible. The EM-tenure model should
match these transition estimates while *also* fitting the duration data.

## The Solution: Free Specification

We **decouple λ from θ**. The transition probabilities are estimated from
the transition counts alone (just like the simple model). The exponential rates
are estimated from the duration data alone. They no longer interfere with each
other.

Concretely, the M-step changes from:

- **Old (linked)**: θ₁ solves a joint equation combining transition counts
  AND tenure distributions → compromise at ~0.50
- **New (free)**: θ₁ = T₁₁/D₁ (closed form, from transitions only);
  λ_g solved separately from EMG emissions only

This means:
- **θ₁ should now match the simple model** (~0.97 with misclassification,
  ~0.91 without)
- **λ_g is free to fit the observed tenure distribution** (long average tenure
  → small λ_g)
- The two are no longer in conflict

## Identification (why this works)

A natural worry: if λ and σ² are both estimated from the duration data,
are they separately identified? Yes, because they come from **different** parts
of the data:

1. **σ²_g** comes from continuation increments: Δg = (g_t − g_{t−1}) − 0.25.
   This involves *only* σ²_g, never λ_g.
2. **λ_g** comes from EMG observations (wave 1, some misclassified cases).
   Given σ²_g from step 1, the EMG has only one free parameter (λ_g), which
   is uniquely identified.
3. **θ₁** comes from transition counts, completely independent of durations.

No ridge, no confounding. We can verify this empirically by running from
multiple starting values and checking they converge to the same answer.

## The CTMC Link Becomes a Testable Restriction

The linked specification isn't gone — it becomes a hypothesis we can test:

    H₀: λ_g = −log(1 − θ₁) / 3 and λ_d = −log(1 − θ₀) / 3

After fitting both the free and linked models, a likelihood ratio test
(2 degrees of freedom) tells us whether the CTMC memoryless property holds.
We expect it to be rejected (tenure has strong duration dependence), which
is itself an interesting empirical finding.

## Never-Worked Population

Previously excluded from the EM sample. Now included by treating them as
continuously nonemployed since age 16:

- Duration = (age − 16) × 12 months
- Mapped to a timegap category using standard binning
- Ages ≥ 21 → category 7 ([60, ∞) months)
- This is conservative: the true nonemployment spell is at least this long

## What Changes in the Code

All changes are confined to `EM-tenure/`:

1. **mstep.R**: Add `linked` flag. When `linked = FALSE`:
   - θ₁ = T₁₁/D₁ (closed form)
   - θ₀ = T₀₁/D₀ (closed form)
   - λ_g via Brent on EMG-only score
   - λ_d via Brent on interval-censored Exp

2. **em_driver.R**: Add `linked` parameter to `em_fit_tenure()`.
   Default is `FALSE` (free specification). `TRUE` gives the old behaviour.

3. **simulate.R**: Allow user-supplied λ_g, λ_d independent of θ.

4. **estimate_pipeline.R**: Remove never-worked filter; add age-based
   timegap imputation.

5. **diagnostics.R** (new): Optional post-estimation functions for
   multi-start, profile likelihood, and numerical Hessian.

## Expected Results

After these changes, running on QLFS data should produce:

| Parameter | With miscl. | Without miscl. |
|-----------|------------|----------------|
| θ₁        | ~0.96–0.97 | ~0.91          |
| θ₀        | ~0.03–0.04 | ~0.08          |
| π         | ~2.5–3%    | 0 (fixed)      |
| λ_g       | small (long tenure) | small   |
| λ_d       | from timegap data   | from data |

These should be consistent with the simple binary-sequence MLE for the
transition parameters, with the added information from duration data
going into the λ estimates.
