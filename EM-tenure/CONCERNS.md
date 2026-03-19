# CONCERNS.md — Methodological concerns and limitations

## 1. Observation-level loop

The E-step loops over individuals in R (`for (i in seq_len(N))`). For large datasets (N > 50k), this will be slow. Potential mitigations:
- Vectorise the E-step using matrix operations (main bottleneck is that each individual has different `s`, `g`, `d` vectors, making full vectorisation non-trivial but possible)
- Port the inner loop to C++ via Rcpp

## 2. EMG numerical stability

The EMG log-density uses the `erfc` function, which can underflow for large arguments. The current implementation uses the `log(erfc(...))` formulation with `pnorm()` under the hood, which is reasonably stable but could be improved with a dedicated `log_erfc` implementation for extreme tails.

## 3. Local optima

The EM algorithm is guaranteed to converge to a local maximum. The observed-data likelihood is a mixture over 8 latent histories, which may have multiple modes. Recommendations:
- Run from multiple random initialisations
- Use the synthetic data generator to validate that the true parameters are recovered

## 4. Identification with small pi

When misclassification is very small (π < 0.01), the model may struggle to distinguish π from measurement noise σ². The sufficient statistics for π come from all waves, while σ² comes only from correctly-classified consecutive observations. With very small π, nearly all responsibilities concentrate on the "correct" history, providing little information to update π.

## 5. Duration positivity

The abs() transformation in the simulator ensures durations are non-negative, but the Normal increment model can produce negative increments. In the real data, negative increments (g_t < g_{t-1}) are informative about measurement error. The model handles this correctly through the Normal density, which assigns positive probability to negative increments.

## 6. Pooled variance formula

The M-step pools continuation increments (variance 2σ²) with within-panel starts (variance σ²) to estimate σ². This is correct under the model but assumes both sources are equally reliable. If starts are noisier (e.g., due to rounding at spell boundaries), the pooled estimate may be biased.

## 7. Stationarity assumption

When `stationary = TRUE`, α is constrained to the ergodic distribution. This is appropriate for long-running labour markets but may be violated during recessions or structural breaks. The non-stationary version (default) allows α to differ from the ergodic value.

## 8. Symmetric misclassification

The model assumes P(s=1|h=0) = P(s=0|h=1) = π. If misclassification is asymmetric, the model will estimate an "average" π. Extending to asymmetric misclassification (π₁₀, π₀₁) is straightforward but adds one parameter and may be poorly identified with only 3 waves.

## 9. Three-wave limitation

With only T=3 waves, the model has 8 latent histories but 2³ = 8 observed state patterns. The continuous duration data provides additional identification, but the transition probability estimates rely on only 2 transitions per individual. Standard errors (not yet implemented) will be large.

## 10. Bootstrap SEs (deferred)

Standard errors via the observed information matrix are complex due to the mixture structure. The recommended approach is non-parametric bootstrap. This is deferred to a future implementation.

## 11. Weight handling

Survey weights enter both the log-likelihood and the sufficient statistics. The M-step updates are weighted analogues of the unweighted MLEs. If weights are highly variable, the effective sample size may be much smaller than N, inflating standard errors.

## 12. Missing duration data

The current implementation requires all 6 duration columns to be non-NA. The ingest script drops rows with any NA tenure/timegap. This is conservative; a more sophisticated approach would handle partial missingness in the E-step.

## 13. Log-likelihood decrease under the linked specification (GEM, not full EM)

> **Status: RESOLVED (2026-03-15).**  The M-step now uses Newton-Raphson to
> solve the joint FOC for θ₁ and θ₀, which includes both the Markov prior
> gradient and the EMG emission gradient.  This makes the M-step a true
> maximisation of the full Q function, restoring the standard EM monotonicity
> guarantee.  The E-step accumulates `emg_g_x/w` and `emg_d_x/w` (weighted
> (x, w) pairs for all EMG observations) and passes them to the Newton solver
> in `mstep.R`.  The monotonicity threshold in `em_driver.R` has been
> tightened to -1e-8 (from the old GEM-tolerant -1e-4 × |LL|).

The historical analysis below is retained for reference.

### Symptom

When running the EM algorithm with the linked specification (λ derived deterministically from θ via the CTMC link), the observed-data log-likelihood can **decrease** on early iterations, sometimes substantially (empirically: drops of 1–5 nats on |LL| ≈ 35 at n = 10). This violates the standard EM monotonicity guarantee.

### Root cause

The standard EM guarantee relies on a key property: the M-step must maximise (or at least not decrease) the complete-data expected log-likelihood Q(Θ | Θ^old) = E_{H|Y,Θ^old}[log p(Y, H | Θ)] with respect to Θ. When this holds, the observed-data log-likelihood L(Θ) is guaranteed to be non-decreasing at each iteration.

The Q function decomposes additively into three blocks:

```
Q(Θ | Θ^old) = Q_prior(α, θ₁, θ₀) + Q_misclass(π) + Q_emission(σ²_g, σ²_d, λ_g, λ_d)
```

The **Markov prior block** Q_prior depends on (α, θ₁, θ₀) and is maximised in closed form by the standard sufficient-statistic formulas:

```
θ̂₁ = T₁₁ / D₁,   θ̂₀ = T₀₁ / D₀,   α̂ = C₁ / (C₁ + C₀)
```

The **emission block** Q_emission depends on (σ²_g, σ²_d, λ_g, λ_d). For σ², the closed-form pooled MLE (Eq 21) correctly maximises Q_emission with respect to σ². For λ, the EMG terms in Q_emission contribute a non-trivial dependence.

Under the **linked specification**, λ_g = −log(1−θ₁)/3 and λ_d = −log(1−θ₀)/3. This means that updating θ in the Markov block **simultaneously changes λ** in the emission block. However, the closed-form M-step for θ only maximises Q_prior — it does not account for the Q_emission terms that depend on λ(θ). When θ changes substantially (as on early iterations), the induced change in λ can decrease Q_emission by more than Q_prior increases, causing Q and hence L to decrease.

### Formally

Let Θ^new denote the M-step update. In a standard (non-linked) EM:

```
Q(Θ^new | Θ^old) ≥ Q(Θ^old | Θ^old)   (M-step maximises Q)
```

Under the linked specification, the M-step maximises only Q_prior + Q_misclass + Q_σ but not Q_λ:

```
Q_prior(θ^new) ≥ Q_prior(θ^old)           ✓ (Markov block maximised)
Q_emission(σ^new, λ(θ^new)) ≶ Q_emission(σ^old, λ(θ^old))   ✗ (not guaranteed)
```

Therefore:

```
Q(Θ^new | Θ^old) ≶ Q(Θ^old | Θ^old)     (full Q may decrease)
```

This makes the algorithm a **Generalised EM (GEM)** rather than a full EM. GEM algorithms converge to a stationary point but do not guarantee monotone ascent of the observed-data log-likelihood at every iteration.

### TeX claim

The TeX document states: *"Monotone ascent is guaranteed because every M-step update maximizes Q in closed form."* This is **incorrect** for the linked specification. It would be correct only if:

1. λ were treated as a free parameter (unlinked), with its own M-step maximising Q_emission with respect to λ independently of θ; or
2. The M-step for θ maximised the full Q(Θ | Θ^old) including the EMG terms — which has no closed-form solution because Q_emission involves λ(θ) inside both the exponential and erfc terms of the EMG density.

### Empirical behaviour

In practice, the LL decrease is largest on the first 1–3 iterations when the initial θ values are far from their optimum and the induced λ change is large. As θ stabilises, the LL decrease becomes negligible and the algorithm converges normally. In all tested cases (synthetic data, n = 10 to n = 500), the algorithm converges to a reasonable solution despite early LL drops.

### Current mitigation

The monotonicity check in `em_driver.R` uses a threshold of −1e-4 × |LL|. Drops within this band are silently accepted. Drops exceeding this threshold trigger a warning. This prevents false alarms from the GEM behaviour while still catching genuine bugs (incorrect sufficient statistics, coding errors).

