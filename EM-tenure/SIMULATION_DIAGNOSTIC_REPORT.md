# Simulation Diagnostic Report

## Purpose

The EM-tenure model produces implausible estimates on real QLFS data (θ₁≈0.50, θ₀≈0.44) despite passing synthetic validation. To understand what the model is fitting, we simulated data from three parameter sets and compared each against the real (never-worked-excluded) data across multiple dimensions.

## Parameter sets

| Parameter | Sim-EM | Sim-Expected | Sim-Adjusted |
|-----------|--------|--------------|--------------|
| α         | 0.775  | 0.47         | **0.693**    |
| θ₁        | 0.50   | 0.95         | **0.961**    |
| θ₀        | 0.44   | 0.05         | **0.120**    |
| π         | 0.14   | 0.03         | 0.03         |
| σ²_g      | 0.12   | 0.01         | 0.01         |

- **Sim-EM**: current EM-tenure estimates (implausible).
- **Sim-Expected**: binary-sequence MLE on the *full* sample (including never-worked).
- **Sim-Adjusted**: recalibrated for the *filtered* population using empirical moments from `df_real`. This is the fair benchmark.

### Population mismatch correction

The binary MLE was estimated on the full sample where ~55% of nonemployed had *never worked*. The EM-tenure pipeline excludes these, shifting:

- **α**: 0.47 → 0.69 (employed count unchanged, nonemployed halved)
- **θ₀**: 0.05 → 0.12 (never-worked were "stuck" in state 0, deflating the exit rate)
- **θ₁**: ~0.95 → ~0.96 (minimal — never-worked don't affect retention)

## Key findings

### 1. Sim-Adjusted matches the binary structure well

| Dimension                      | Sim-Adjusted | Sim-EM | Sim-Expected |
|--------------------------------|-------------|--------|--------------|
| Employment rate (MAD)          | **0.024**   | 0.159  | 0.239        |
| Sequence pattern JSD           | **0.004**   | 0.206  | 0.036        |
| θ₁ naive (ppt diff from Real)  | **1.4**    | 44.0   | 4.4          |
| θ₀ naive (ppt diff from Real)  | **4.7**    | 32.5   | 4.6          |

The Markov transition model is approximately correct for this population.

### 2. No parameter set matches tenure

| Statistic    | Real  | Sim-Adjusted | Sim-EM | Sim-Expected |
|-------------|-------|-------------|--------|--------------|
| Mean (yrs)  | 6.53  | 1.12        | 3.69   | 1.18         |
| Median (yrs)| 4.00  | 0.86        | 2.25   | 0.90         |
| Max (yrs)   | 76.8  | 11.5        | 60.3   | 12.5         |

With θ₁=0.96, the CTMC link forces λ_g ≈ 1.07, giving E[tenure] ≈ 0.93 years. The exponential cannot produce the observed heavy-tailed distribution.

### 3. Tenure increments reveal massive measurement error

| Statistic   | Real | Sim-Adjusted | Sim-EM |
|------------|------|-------------|--------|
| SD          | 3.38 | 0.16        | 2.80   |
| Implied σ_g | 2.39 | 0.11        | 1.98   |

Real increment SD is 3.38 years — orders of magnitude larger than σ²_g=0.01 predicts. The EM inflates σ²_g to 0.12 to compensate, distorting everything else.

### 4. Timegap categories are non-monotone

| Cat | Real  | Sim-Adjusted | Sim-EM |
|-----|-------|-------------|--------|
| 1   | 12.1% | 6.1%        | 42.3%  |
| 5   | 18.3% | 7.9%        | 13.1%  |
| 7   | 40.2% | 74.3%       | 17.9%  |

Real data has a mode at category 7 (60+ months) — not the monotonically decreasing shape that Exp(λ_d) predicts. The exponential fails on the nonemployment side too.

## Diagnosis

| Dimension               | Sim-Adjusted | Sim-EM    |
|--------------------------|-------------|-----------|
| Employment rate          | ✅ 0.024     | ❌ 0.159   |
| Sequence patterns        | ✅ 0.004     | ❌ 0.206   |
| θ₁ naive                 | ✅ 1.4 ppt   | ❌ 44.0 ppt|
| Timegap JSD              | ✅ 0.070     | ❌ 0.084   |
| **Tenure mean**          | ❌ 5.41 yrs  | ✅ 2.84 yrs|
| **Tenure increment SD**  | ❌ 3.22      | ✅ 0.58    |

**Sim-Adjusted wins on transitions, Sim-EM wins on durations.** The EM compromises between two irreconcilable objectives:

1. **Binary sequence likelihood** pulls θ₁ → 0.96 (correct)
2. **Duration likelihood** pulls θ₁ → 0.50 (to lengthen expected tenure via the CTMC link)

The CTMC link forces these to share parameters, so the EM splits the difference.

## Conclusion

- The **Markov transition model is correct** for this population (α≈0.69, θ₁≈0.96, θ₀≈0.12, π≈0.03).
- The **exponential duration assumption fails** on both sides: tenure is heavy-tailed, timegap categories are non-monotone.
- The **CTMC link** is the binding constraint — forcing θ to explain both transitions and durations causes the EM to converge to an interior compromise wrong for both.
- The **population mismatch** (never-worked exclusion) was an important correction — without it, Sim-Expected appeared wrong on *everything*, masking the clean separation between binary-structure fit and duration-model failure.

## Next steps

1. Drop the CTMC link — estimate λ_g, λ_d as free parameters
2. Replace exponential with a heavier-tailed distribution for tenure
3. Re-examine σ²_g — real increment SD of 3.38 years suggests the continuation model may itself be misspecified
4. Re-run binary MLE on the filtered sample for exact adjusted θ values
