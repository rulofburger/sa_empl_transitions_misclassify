# EM-tenure Troubleshooting: Diagnosis

**Date:** 2026-03-19  
**Model:** `EM-tenure/` — EM algorithm for employment transitions with tenure
and nonemployment durations  
**Symptom:** Estimated θ₁ ≈ 0.465, θ₀ ≈ 0.418, π ≈ 0.266 on QLFS data.
Expected: θ₁ ≈ 0.97, θ₀ ≈ 0.03, π ≈ 0.03 (from simpler binary-sequence
models and prior literature). Observed raw quarterly transition rate ≈ 8–9%.

---

## Data Terminology

Throughout this document, three data states are distinguished:

| Term | Description |
|------|-------------|
| **Raw data** | `data/raw/df_qlfs_A.rds` — tenure in months (integers), timegap as labelled factor codes (0–8, 99) |
| **Post-ingest data** | After `scripts/ingest_data_3waves_SA.R` — all variables in years, timegap midpoints, zeros for wrong state |
| **Post-filter data** | After `EM-tenure/estimate_pipeline.R` ad-hoc filter — rows with (y=1 & g≤0) or (y=0 & d≤0) dropped |

The EM-tenure model operates exclusively on **post-ingest** data after the
pipeline filter is applied. All issues below refer to **post-ingest** data
unless explicitly stated.

---

## Post-Ingest Data Encoding

The ingest script produces this structure for the duration variables:

| Variable | When employed (y=1) | When nonemployed (y=0) |
|----------|---------------------|------------------------|
| `tenure` | Raw months / 12 (years). Negative raw values → `NA` → row dropped | **Forced to 0** by `if_else(y==0, 0, tenure)` |
| `timegap` | **Forced to 0** by `if_else(y==1, 0, timegap)` | Midpoint of survey category (months) / 12, or `(age-educ-6)` years for never-worked |

**Consequence:** Every observation carries a zero in the "wrong" duration
variable. This is intentional for storage (you can't have both tenure and
timegap simultaneously) but is catastrophic for the EM model as shown below.

---

## Issue 1 — Zero Durations for the "Wrong" State (CRITICAL)

### What happens

The ingest script sets `timegap = 0` for all employed observations and
`tenure = 0` for all nonemployed observations. This is the right duration for
the _observed_ state, but the EM model must also evaluate durations for
_latent_ histories that differ from the observed state (misclassification).

### How it breaks the EM

Under the emission model, when the latent history has `h_t = 0` (truly
nonemployed) but the observation is `s_t = 1` (observed as employed), the
E-step evaluates:

```
log_emg(g_t, lambda_g, sigma2_g)
```

where `g_t` is the observed tenure. For a truly nonemployed person, `g_t = 0`
(set by the ingest script). Since `log_emg(0, ...) = -Inf`, any latent history
requiring misclassification at that wave receives log-probability `-Inf` —
**zero posterior weight regardless of the parameter values**.

The symmetric argument applies to misclassified nonemployed (h=1, s=0):
`log_emg(d_t, ...) = -Inf` because `d_t = 0` for employed people.

### Consequence chain

```
zero duration for wrong state
  → log_emg(0) = -Inf for all misclassifying histories
  → posterior γ_{ih} = 0 for all h with misclassification at that wave
  → sufficient statistic M (misclassification count) ≈ 0
  → M-step: π̂ = M / (3·Σwᵢ) ≈ 0
```

But π is initialized to 0.03 and the binary-sequence information still
pushes against π = 0. The EM therefore reaches a compromise that is
inconsistent with both the duration data and the binary-sequence data,
ultimately converging to a local optimum far from the true parameters.

The model cannot distinguish "misclassification" from "measurement noise"
when the inappropriate duration is exactly zero — it collapses all such
observations to histories with no misclassification.

### Scope

- Affects **every** wave for **every** observation where misclassification
  is plausible
- Wave 1: employed person (s₁=1) with timegap₁=0 → histories (h₁=0,*)
  have `log_emg(0, λ_d, σ²_d) = -Inf`
- Waves 2–3: same logic in `.wave_emission_vec()` for the misclassified branch

### Remedy: Nearest-non-zero imputation with floor fallback

Replace zero durations for the wrong state with the **nearest non-zero
value** of the same variable from another wave of the same individual.
If no non-zero value exists at any wave, use a very small floor.

**Algorithm:** For each individual _i_ and each wave _t_ where the
wrong-state duration is zero:

1. **Look at the same variable at other waves.** For `timegap_t = 0`
   (because y_t = 1), look at `timegap_{t'}` for t' ≠ t where
   `timegap_{t'} > 0`. Choose the nearest wave (|t' − t| = 1 preferred
   over |t' − t| = 2).

2. **If a non-zero value exists**, use it as the imputed wrong-state
   duration. Rationale: if the person is truly nonemployed but observed
   as employed, their nonemployment duration is plausibly similar to
   the duration they report at a nearby wave where they were correctly
   classified as nonemployed.

3. **If all three waves have zero** for that variable (e.g., the person
   is observed as employed at all three waves), use a very small floor:

   ```
   floor = 0.25 / 12 ≈ 0.021 years ≈ 0.25 months ≈ 1 week
   ```

   This is deliberately small. A person employed at all three waves with
   zero timegap throughout has very low plausible nonemployment duration
   if misclassified. The small floor makes the EMG density finite (so
   misclassification is not structurally ruled out) but assigns low
   likelihood (so the EM does not over-attribute misclassification).

4. **Symmetric treatment for tenure:** For `tenure_t = 0` (because
   y_t = 0), look at `tenure_{t'}` for t' ≠ t where `tenure_{t'} > 0`.
   If none exist, use the same floor (0.021 years).

**Example:**

| Wave | y | timegap (raw) | tenure (raw) | timegap (imputed) | tenure (imputed) |
|------|---|--------------|-------------|-------------------|------------------|
| 1 | 0 | 0.375 | 0 | 0.375 | **2.5** (from wave 2) |
| 2 | 1 | 0 | 2.5 | **0.375** (from wave 1) | 2.5 |
| 3 | 0 | 0.625 | 0 | 0.625 | **2.5** (from wave 2) |

**Note on interaction with the discrete/censored timegap model (Issue 3):**
Under the discrete model, the imputation provides a **timegap category code**
(integer 1–7) rather than a continuous midpoint. The nearest-non-zero logic
is the same but operates on category codes. If all three waves have y = 1,
assign category 1 ([0, 0.25) years) as the floor — the shortest possible
nonemployment spell.

---

## Issue 2 — Never-Worked Individuals (~50% of Sample) (CRITICAL)

### What happens in the raw data

The QLFS survey asks nonemployed respondents how long they have been
searching/out of work. For individuals who have **never held a job**
(`neverworked = 1`), this question is not applicable (timegap code = 0,
"Not applicable"). The ingest script imputes their nonemployment duration as:

```r
timegap = (age - educ - 6) / 12  # years since leaving school, in years
```

This produces values ranging from ~5 years (young school leaver) to ~30+
years (older never-worker).

### Post-ingest distribution

- `neverworked1` mean ≈ 0.495 (approximately 50% of the sample after ingest)
- Imputed timegaps for never-worked: typically 5–30 years
- The EM model's EMG distribution `Exp(λ_d) * N(0, σ²_d)` expects:
  - Mean ≈ 1/λ_d + ε, where λ_d is determined by θ₀ ≈ 0.03
  - At θ₀ = 0.03: λ_d = -log(0.97)/3 ≈ 0.0101 per month ≈ 0.121 per year
  - Mean of EMG ≈ 1/0.121 ≈ 8.3 years
- This might seem fine, but the quarterly increment model then expects
  increments of ~0.25 years with variance ~2σ²_d

### Consequence chain

```
never-worked timegap = 10–30 years
  → EMG(λ_d, σ²_d) requires λ_d ≈ 0.05–0.1 (mean ≈ 1/λ ≈ 10–20 years)
  → θ₀ = 1 - exp(-λ_d × 3) ≈ 0.14–0.26  (to fit the long durations)
  
  BUT:
  
  binary sequences show θ₀ ≈ 0.03
  → tension: duration data wants high θ₀, binary data wants low θ₀
  → EM compromises at θ₀ ≈ 0.42 (observed result)
```

Additionally, the increment model for nonemployment continuations expects:
```
Δd_t = d_t - d_{t-1} - 0.25 ~ N(0, 2σ²_d)
```
But for never-worked individuals, `d_t ≈ d_{t-1}` (same person, one quarter
later, same "years since school"), so Δd_t ≈ 0 ≠ 0.25. The systematic
mismatch inflates σ²_d and distorts the M-step for θ₀.

### Scope

- Approximately 50% of the post-ingest sample
- Concentrated entirely in the nonemployed subsample (since they have y=0)
- The never-worked imputation is applied at **each wave independently**, so a
  person who is never-worked at all three waves contributes three anomalous
  duration observations

### Remedy

Filter out never-worked individuals before passing to EM:

```r
df_em <- df_qlfs |>
  filter(neverworked1 == 0 | is.na(neverworked1))
```

Note: `neverworked` is NA for employed individuals in the raw data (the
question is not asked). After the ingest script, these NAs remain. The filter
should preserve employed individuals (neverworked = NA) and exclude only
those with neverworked = 1. Alternatively:

```r
df_em <- df_qlfs |>
  filter(!(!is.na(neverworked1) & neverworked1 == 1))
```

The model was designed for workers who have labour market histories, not for
long-term labour market non-participants.

---

## Issue 3 — Categorical Timegap Data (Gaussian Violation) (HIGH)

### What happens

The raw QLFS timegap variable is categorical with 8 buckets. The ingest
script maps these to midpoint values in months, then divides by 12 to get
years. The post-ingest values are therefore drawn from the discrete set
{0, 0.125, 0.375, 0.625, 0.875, 2.0, 4.0, 7.5} years.

### Model assumption violated

The EM model assumes that nonemployment duration increments follow:
```
Δd_t = d_t - d_{t-1} - 0.25 ~ N(0, 2σ²_d)
```

For categorical midpoint data, `d_t` is one of 8 discrete values and
`d_{t-1}` is one of 8 discrete values. The increment `Δd_t` takes values
from a discrete grid (differences of multiples of 1.5/12, 4.5/12, ...).
This is not Gaussian.

### Specific distortions

1. **Point mass at zero**: Many nonemployed people stay in the same bucket
   across waves (e.g., category 1 → category 1, both mapped to 0.125 years).
   Their increment is `0.125 - 0.125 - 0.25 = -0.25`. This is a systematic
   non-zero residual that inflates σ²_d.

2. **Large discrete jumps**: If someone transitions from bucket 4 (10.5/12
   ≈ 0.875 years) to bucket 5 (24/12 = 2.0 years), the increment is
   `2.0 - 0.875 - 0.25 = 0.875 years`. This is a large outlier under
   N(0, 2σ²_d) and distorts σ²_d upward.

3. **σ²_d inflation → θ₀ distortion**: Inflated σ²_d affects the joint
   FOC for θ₀ (via the EMG gradient term in the Newton step), pulling
   θ₀ away from its true value.

### Scope

Affects all post-ingest timegap values for the approximately 50% of
observations that are nonemployed. The distortion is systematic but
smaller in magnitude than Issues 1 and 2.

### Remedy: Discrete interval-censored timegap model

Replace **all** Gaussian and EMG emission densities on the nonemployment
(timegap) side with an **interval-censored exponential model** that works
directly with the survey's category structure. The tenure (employment)
emission model is unchanged.

#### Underlying model

The CTMC-linked exponential is preserved: true nonemployment duration
D ~ Exp(λ_d), where λ_d = −log(1 − θ₀)/3. What changes is how the
**observed** timegap data enters the likelihood: instead of evaluating a
continuous density at a midpoint, we compute the **probability that the
true duration falls in the observed category interval**.

#### Category intervals

Let K = 7 categories with boundaries [a_k, b_k) in years:

| k | Interval (months) | Interval (years) | a_k | b_k |
|---|-------------------|-------------------|-----|-----|
| 1 | [0, 3) | [0, 0.25) | 0 | 0.25 |
| 2 | [3, 6) | [0.25, 0.5) | 0.25 | 0.5 |
| 3 | [6, 9) | [0.5, 0.75) | 0.5 | 0.75 |
| 4 | [9, 12) | [0.75, 1.0) | 0.75 | 1.0 |
| 5 | [12, 36) | [1.0, 3.0) | 1.0 | 3.0 |
| 6 | [36, 60) | [3.0, 5.0) | 3.0 | 5.0 |
| 7 | [60, ∞) | [5.0, ∞) | 5.0 | ∞ |

#### The six emission cases involving timegap

A complete audit of the TeX document, E-step code, M-step code, emission
functions, and simulator confirms exactly **six emission cases** where
λ_d and/or σ²_d appear. All six are replaced:

**Case 1: Wave 1 matched nonemployment** (t=1, h₁=0, s₁=0)

- Current: ℓ₁ = f_EMG(d₁; λ_d, σ²_d)
- Proposed: ℓ₁ = P(D ∈ [a_c, b_c) | λ_d) = exp(−λ_d · a_c) − exp(−λ_d · b_c)
- Justification: memoryless property of exponential → residual life is
  Exp(λ_d). The interval probability is the natural likelihood for
  interval-censored data.

**Case 2: Wave 1 misclassified as nonemployed** (t=1, h₁=1, s₁=0)

- Current: ℓ₁ = f_EMG(d₁; λ_d, σ²_d)  [returns −∞ when d₁ = 0]
- Proposed: ℓ₁ = P(D ∈ [a_c, b_c) | λ_d)  [same as Case 1]
- Note: combined with nearest-non-zero imputation (Issue 1), the
  employed person now has an imputed timegap category, so the −∞
  problem is eliminated.

**Case 3: Nonemployment continuation, previous correctly observed**
(t ≥ 2, h_{t-1}=0, h_t=0, s_t=0, s_{t-1}=0)

- Current: ℓ_t = ϕ(Δd_t / √(2σ_d)) / (√2 σ_d)  — Gaussian increment
- Proposed: ℓ_t = P(c_t | c_{t-1}, λ_d)  — conditional category transition
  probability

The transition probability is:

```
P(c_t | c_{t-1}) = P(D_{t-1} + 0.25 ∈ [a_{c_t}, b_{c_t}) | D_{t-1} ∈ [a_{c_{t-1}}, b_{c_{t-1}}))
                 = [exp(−λ_d · L) − exp(−λ_d · U)] / [exp(−λ_d · a_j) − exp(−λ_d · b_j)]

where L = max(a_j, a_k − 0.25), U = min(b_j, b_k − 0.25), j = c_{t-1}, k = c_t
and b_j = ∞ for rightmost category.
```

The 7×7 transition matrix is **sparse** — at most 2 non-zero entries per
row. Categories 1–4 have **deterministic** transitions (the narrow 3-month
intervals shift exactly one category forward with a 0.25-year increment):

| From | After +0.25 yr | Destination |
|------|---------------|-------------|
| 1: [0, 0.25) | [0.25, 0.5) | 2 (deterministic) |
| 2: [0.25, 0.5) | [0.5, 0.75) | 3 (deterministic) |
| 3: [0.5, 0.75) | [0.75, 1.0) | 4 (deterministic) |
| 4: [0.75, 1.0) | [1.0, 1.25) | 5 (deterministic) |
| 5: [1.0, 3.0) | [1.25, 3.25) | 5 or 6 (probabilistic) |
| 6: [3.0, 5.0) | [3.25, 5.25) | 6 or 7 (probabilistic) |
| 7: [5.0, ∞) | [5.25, ∞) | 7 (deterministic) |

**Case 4: Nonemployment continuation, previous misclassified**
(t ≥ 2, h_{t-1}=0, h_t=0, s_t=0, s_{t-1}=1)

- Current: ℓ_t = f_EMG(d_t; λ_d, σ²_d)
- Proposed: ℓ_t = P(D ∈ [a_c, b_c) | λ_d)  — marginal interval probability
  (same as Cases 1/2, since previous category is unobserved)

**Case 5: Within-panel nonemployment start**
(t ≥ 2, h_{t-1}=1, h_t=0, s_t=0)

- Current: ℓ_t = ϕ((d_t − 0.25) / σ_d) / σ_d  — Gaussian level
- Proposed: ℓ_t = 1{c_t = 1}  — deterministic
- Justification: a new spell of at most 0.25 years must be in category 1
  ([0, 0.25) years = [0, 3) months). Any other category is structurally
  impossible.

**Case 6: Misclassified as nonemployed, t ≥ 2** (h_t=1, s_t=0)

- Current: ℓ_t = f_EMG(d_t; λ_d, σ²_d)
- Proposed: ℓ_t = P(D ∈ [a_c, b_c) | λ_d)  — marginal interval probability

#### Parameter elimination

σ²_d **is eliminated from the model**. The nonemployment measurement
variance no longer exists as a parameter because:

- The categorisation *is* the survey's observation mechanism — there is no
  continuous measurement corrupted by Gaussian noise.
- Continuation increments now contribute category transition probabilities
  that depend only on λ_d.
- Within-panel starts have deterministic emission (no parameters).

The model drops from 6 to 5 parameters:
(α, θ₀, θ₁, π, σ²_g).

#### M-step implications

- **θ₀ FOC**: The EMG gradient term E_d is replaced by an
  interval-probability gradient:
  ```
  ∂/∂θ₀ log P(D ∈ [a_k, b_k)) = [1/(3(1−θ₀))] · (−a_k exp(−λa_k) + b_k exp(−λb_k)) / (exp(−λa_k) − exp(−λb_k))
  ```
  The Brent solver structure remains the same — only the gradient function
  changes.
- **σ²_d**: Removed. Sufficient statistics S_d, N_d, S_d^start, N_d^start
  are no longer computed.
- **θ₁**: Unchanged (still uses EMG gradient on tenure side).
- **σ²_g**: Unchanged.

#### Theoretical advantages

1. **Respects the DGP**: the survey produces categories, not continuous
   measurements. The interval-censored model is the correct likelihood.
2. **Eliminates σ²_d**: removes a parameter that was not identified by the
   actual observation mechanism.
3. **Eliminates the Gaussian increment artefact**: same-category stays no
   longer produce a systematic −0.25 residual.
4. **More informative starts**: the deterministic c_t = 1 emission is exact.
5. **Sparse transition matrix**: at most 9 non-zero entries in the 7×7 matrix.
6. **CTMC link preserved**: Exp(λ_d) with λ_d = −log(1−θ₀)/3 is unchanged.

---

## Issue 4 — Wave 1 Emission for Misclassified Employed at Wave 1 (HIGH → RESOLVED by Issues 1+3)

### What happens

For wave 1, the E-step evaluates EMG for all duration observations
regardless of whether they are matched or misclassified (since all wave-1
durations are "left-censored"). For a history `h₁ = 1` (truly employed at
wave 1) with `s₁ = 0` (observed as nonemployed), the code evaluates:

```r
log_emg(d1[mask_non], lambda_d, sigma2_d)
```

where `d1 = timegap1`. For an actually employed person (`y₁ = 1`), the
ingest script sets `timegap1 = 0`. So:

```r
log_emg(0, lambda_d, sigma2_d) = -Inf
```

This makes **history (1,*,*) with s₁ = 0 impossible**, meaning
misclassification at wave 1 for an employed person is structurally ruled out
by the data encoding — not by the model.

### Impact on π estimation

If misclassification at wave 1 is impossible for employed people (due to
zero timegap), then the only way π enters the likelihood is through
misclassification at waves 2 and 3. This halves the effective information
about π and biases estimates.

More critically: the restriction is asymmetric. Employed people cannot be
misclassified as nonemployed at wave 1, but nonemployed people could (in
principle) be misclassified as employed — except that Issue 1 applies to
them too at waves 2 and 3. The net effect is that the EM sees essentially
no misclassification signal.

### Resolution

This issue is a special case of Issue 1 and is **fully resolved** by the
combination of the Issue 1 and Issue 3 remedies:

- **Issue 1 fix** (nearest-non-zero imputation): the employed person at
  wave 1 now has an imputed timegap category from a nearby wave, so the
  duration is no longer zero.
- **Issue 3 fix** (discrete/censored model): the emission density is now
  P(D ∈ [a_c, b_c) | λ_d) — an interval probability that is always
  positive for any category c ∈ {1,...,7}. Even if the imputed category
  is c = 1 (the floor), the probability 1 − exp(−λ_d · 0.25) > 0.

No separate fix is needed for Issue 4.

---

## Issue 5 — Potential Double-Counting in `.wave_emission_vec()` (MEDIUM, LIKELY BENIGN)

### What to check

In `estep.R`, the function `.wave_emission_vec()` handles each history `j`
with two blocks: "correctly classified at wave t" and "misclassified at
wave t". For the correctly classified block, emissions are added for
`h_prev==1, h_curr==1` (employment continuation) etc. For the misclassified
block, emissions are added when `s_t != h_t`.

The concern is: for the case `h_t = 1, s_t = 1` (matched employment), does
the code accidentally also enter the misclassified branch and add a
misclassification emission?

**Analysis:** The misclassified-as-nonemployed branch is:
```r
} else {  # hc == 1
  mask <- (s_t == 0)
  ...
}
```
When `s_t == 1`, `mask = (s_t == 0)` is `FALSE`, so nothing is added.
Similarly for the misclassified-as-employed branch.

**Conclusion:** No double-counting. The code correctly handles the cases.
However, note that the "correctly classified" block conditions on `h_curr`
and `h_prev` but the "misclassified" block conditions only on `h_curr`. This
means that for `h_prev == 0, h_curr == 1` (employment start), the correctly
classified block adds `log_emission_start_g(g_t)` when `s_t == 1`, but the
misclassified-as-nonemployed block adds `log_emg(d_t, ...)` when `s_t == 0`.
This is correct per the TeX specification.

---

## Issue 6 — Ad-Hoc Filter Masks the Root Problem (LOW)

### What happens

The `estimate_pipeline.R` ad-hoc filter drops rows where `(y==1 & g<=0)` or
`(y==0 & d<=0)`. From the session log: this drops 7,539 rows (1.9%).

### What it does NOT fix

The filter removes rows where the observed duration is zero for the **observed
state**. But the critical issue (Issues 1 and 4) is that the duration for the
**wrong** state is zero. For example:

- A person with `y₁=1, tenure₁=0.5, timegap₁=0` passes the filter
  (tenure > 0 ✓), but `timegap₁=0` still causes `-Inf` for any misclassifying
  history at wave 1.
- A person with `y₁=0, timegap₁=0.125, tenure₁=0` passes the filter
  (timegap > 0 ✓), but `tenure₁=0` still causes `-Inf` for histories with
  `h₁=1`.

The filter is necessary but not sufficient. It prevents the E-step validity
check from throwing an error, but the zeros for the wrong state remain.

---

## Consequence Summary

| Issue | Primary effect | Estimated contribution | Remedy status |
|-------|---------------|----------------------|---------------|
| 1: Zero wrong-state durations | Kills misclassification signal; π → distorted | **Critical** | → Nearest-non-zero imputation |
| 2: Never-worked ~50% | Pulls θ₀ up dramatically to fit long durations | **Critical** | → Exclude never-worked |
| 3: Categorical timegaps | Inflates σ²_d; distorts θ₀ via joint FOC | High | → Discrete/censored model |
| 4: Wave 1 asymmetry | Asymmetric information about π | High | → Resolved by Issues 1+3 |
| 5: Double-counting | None (no bug found) | None | N/A |
| 6: Ad-hoc filter | Insufficient | Low | Superseded by Issues 1+3 |

The dominant issues are 1 and 2. Issue 1 structurally prevents the EM from
learning π, and Issue 2 forces θ₀ upward to fit the never-worked duration
distribution. Together they explain why the model converges to
θ₁ ≈ θ₀ ≈ 0.42 — the EM collapses to a near-uniform mixing model where
the distinction between employed and nonemployed is blurred.

---

## Proposed Remedies (Priority Order)

### Remedy R1: Nearest-non-zero imputation (addresses Issues 1 & 4)

Replace zero wrong-state durations with the nearest non-zero value of the
same variable from another wave of the same individual. If no non-zero
value exists at any wave, use a floor of 0.25/12 ≈ 0.021 years.

See Issue 1 Remedy section above for the full algorithm and example.

**Theoretical justification:** A misclassified observation has an
"inappropriate" duration that should reflect a plausible spell length.
The nearest-non-zero approach exploits the within-individual panel structure.
The floor (≈ 1 week) is deliberately small: a person employed at all three
waves has very low plausible nonemployment duration if misclassified.

### Remedy R2: Exclude never-worked individuals (addresses Issue 2)

Filter to individuals who have ever worked, or at minimum who have a valid
categorical timegap (not imputed from age-educ-6):

```r
df_em <- df_qlfs |>
  filter(is.na(neverworked1) | neverworked1 == 0) |>
  filter(is.na(neverworked2) | neverworked2 == 0) |>
  filter(is.na(neverworked3) | neverworked3 == 0)
```

Note: `neverworked` is `NA` for employed individuals (the question is not
asked), so filtering `neverworked != 1` preserves all employed people.

**Theoretical justification:** The model assumes a Markov chain between
employment and nonemployment. Never-worked individuals have never entered the
Markov chain; their "nonemployment duration" is not the same quantity as the
time since last job for a displaced worker. The model is misspecified for them.

### Remedy R3: Discrete/censored timegap model (addresses Issue 3)

Replace all 6 Gaussian/EMG emission densities on the nonemployment side with
an interval-censored exponential model. This eliminates σ²_d from the model
(6 → 5 parameters) and resolves the Gaussian violation.

See Issue 3 Remedy section above for the full specification of all 6 emission
cases, the transition matrix structure, the M-step implications, and the
theoretical advantages.

**Theoretical justification:** The QLFS timegap is a categorical variable
with known interval boundaries. The interval-censored Exp(λ_d) model is the
correct likelihood for this data type. It preserves the CTMC link
(λ_d = −log(1−θ₀)/3) while respecting the discrete observation mechanism.

### Previous Remedy R3 (superseded): Restrict timegap range

The previous approach of restricting to timegap < 1 year was a partial
workaround. It is superseded by the discrete/censored model (new R3), which
handles all categories correctly without sample restriction. It may still be
useful as a **diagnostic subgroup analysis** to verify that the model performs
well on the narrowest categories (1–4) where transitions are deterministic.

---

## Verification Plan

The companion script `troubleshoot/diagnose_em_tenure.R` implements:

1. **Synthetic validation** (§1) — confirms code is correct on clean data
2. **Real data diagnostics** (§2) — quantifies the above issues
3. **Progressive stress tests** (§3) — isolates which pathology drives the
   distortion by introducing them one at a time to synthetic data
4. **Targeted re-estimation** (§4) — applies R1+R2 to real data and compares
   to simpler model baseline (θ ≈ 2–3%, π ≈ 3%)
5. **Discrete/censored model diagnostics** (§5) — pre-implementation
   validation of the proposed interval-censored model:
   - Category distribution of timegap in QLFS data (after excluding
     never-worked)
   - Transition matrix: observed category-to-category transitions vs
     theoretical predictions from the sparse Exp(λ_d) model
   - Comparison of interval probabilities P(D ∈ [a_k, b_k)) vs EMG
     densities evaluated at midpoints, showing how midpoint-based
     densities distort the likelihood
   - Within-panel start category verification: what fraction of
     new nonemployment spells fall in category 1 (expect ~100%)
   - Nearest-non-zero imputation: distribution of imputed values and
     how they compare to the original zeros
