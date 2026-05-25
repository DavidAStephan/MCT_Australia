# Gibbs port plan (Phase 4 / option B from plan.md)

Port the NY Fed MCT MATLAB sampler (`MCT-NYFED/`) to R. Replaces the
Stan-based estimator with a hand-rolled Gibbs sampler that uses
Kim-Shephard mixture-of-normals for stochastic volatility and Carter-Kohn
forward-filter backward-sample for the latent paths. The lit-backed
expected runtime is **10-25 min per single chain at T=435 N=11**, vs
Variant A's current 93 min.

This is a multi-day implementation. Do not start unprompted — see
"decision triggers" below.

## Why Gibbs beats Stan here

Per `plan.md` Phase 3 session log: Stan's autodiff overhead on the
matrix-shaped state of our MCT model is the binding constraint. The KF
marginalization (Akf/Akfs) is mathematically correct and 4× faster than
the naïve dense version, but it still loses to Ac by ~8× because each
filter step has O(D²) matrix ops that Stan's reverse-mode autodiff has
to tape-and-backprop. Conjugate Gibbs updates skip autodiff entirely —
each parameter block is sampled from its full conditional in O(1) using
direct linear algebra.

NY Fed's specific tricks:
- **Kim-Shephard mixture** (1998): the log-χ²₁ distribution arising in
  log-volatility updates is approximated by a 10-component mixture of
  normals. Indicator variables select the component, making the SV
  update a conditionally-Gaussian Carter-Kohn step.
- **Multi-move sampling**: latent factor and sector states sampled in
  one block via forward-backward, breaking the c-s funnel that Stan's
  HMC struggles with.

## Source material

Cached at `MCT-NYFED/` (read-only, ~2,400 lines MATLAB). Key files:

| File | Lines | Purpose |
|---|---|---|
| `functions/estimate_MCT.m` | 487 | Main Gibbs loop; orchestrates everything |
| `functions/Kalman_filter.m` | 175 | Forward filter for the latent factor |
| `functions/simulation_smoother.m` | 70 | FFBS — draws latent paths conditional on θ |
| `functions/fast_smoother.m` | 232 | Efficient smoother (likely de Jong-Penzer) |
| `functions/update_vol.m` | 158 | Kim-Shephard mixture for SV updates |
| `functions/update_tvcoef.m` | 101 | Time-varying loadings (RW) Gibbs draws |
| `functions/update_theta.m` | 34 | Hyperparameter conjugate updates |
| `functions/update_gam.m` | 26 | Common factor parameters |
| `functions/update_ps.m` | 23 | Outlier model component probabilities |
| `functions/update_scl.m` | 40 | Outlier model scales |
| `functions/slice_sampling.m` | 68 | Slice sampler (non-conjugate steps) |
| `functions/decompose_update.m` | 236 | Revisions attribution |
| `functions/simulate_SSM.m` | 138 | Simulator (DGP) |
| `automated_PCE.m` | 653 | Top-level script (data prep + sampler call + outputs) |

Reference papers (cited in NY Fed README):
- Stock & Watson (2016), *Core Inflation and Trend Inflation*, REStat 98(4)
- Kim, Shephard, Chib (1998), *Stochastic volatility: likelihood inference
  and comparison with ARCH models*, REStud 65(3)
- Carter & Kohn (1994), *On Gibbs sampling for state space models*, Biometrika 81(3)

## Architecture for the R port

**Top-level interface**: mirror `fit_mct()`'s shape so the rest of the
project (postprocess, dashboard, compare_*_plot) doesn't change.

```r
fit_mct_gibbs(
  stan_data,                      # same shape as fit_mct
  variant = c("A", "B"),          # AR(1) common or RW common
  n_iter = 5000,                  # post-burn; Gibbs needs more iters than HMC
  n_burn = 5000,
  n_chains = 4,                   # ran in parallel via mclapply
  seed = 20260524,
  save_dir = "outputs/draws"
) -> mct_gibbs_fit
```

**Output object**: an S3 class `mct_gibbs_fit` that exposes the
methods downstream code already uses:

| Method | Returns | Used by |
|---|---|---|
| `$draws(variables, format)` | posterior::draws_* | All postprocess fns |
| `$diagnostic_summary()` | list(num_divergent, num_max_treedepth, ebfmi) | fit_mct watchdog |
| `$save_object(file)` | side effect: writes RDS | scripts/run_real_fits.R |

For `num_divergent` and `num_max_treedepth`, return 0 (Gibbs has neither).
For `ebfmi`, compute as `var(energy)/(2 mean(energy_diff^2))` analogue or
just return 1.0 (Gibbs doesn't use HMC energy; the diagnostic is a no-op).

Downstream targets (`trend_path`, `common_share_path`, etc.) inspect
`$draws(variables = "trend")` etc. As long as the Gibbs fit exposes
draws of the SAME named generated quantities (`trend[T]`,
`common_share[T]`, `sector_trend[T,N]`, `log_lik[n_obs]`), downstream
code is unchanged.

## Port order (low risk first; each step is independently validatable)

### Step 1 — Simulator (`simulate_SSM.m` → `R/gibbs/sim_ssm.R`) ✅ DONE 2026-05-26

Ported as [`R/gibbs/sim_ssm.R`](R/gibbs/sim_ssm.R). Tests at
[`tests/testthat/gibbs/test-sim-and-kalman.R`](tests/testthat/gibbs/test-sim-and-kalman.R)
verify output shape + that simulated AR(1) data hits the right
stationary variance (~10% MC tolerance).

### Step 2 — Kalman filter (`Kalman_filter.m` → `R/gibbs/kalman.R`) ✅ DONE 2026-05-26

Ported as [`R/gibbs/kalman.R`](R/gibbs/kalman.R). Tests verify
log-likelihood matches the analytic joint Gaussian density on a scalar
AR(1) SSM, both with and without missing observations.

**Convention difference from MATLAB**: the MATLAB version uses a fixed
`-(T·N/2)·log(2π)` constant regardless of missing data. The R port
accumulates the constant per-step based on actual non-missing counts,
so the returned value is the joint log-density of observed data only
(matches `dmvnorm`). If we want to cross-check against MATLAB output
later, add back `0.5 · n_missing · log(2π)` to the R port's value.

### Step 3 — Simulation smoother (`simulation_smoother.m` → `R/gibbs/sim_smoother.R`) ✅ DONE 2026-05-26

Ported via two files:

- [`R/gibbs/fast_smoother.R`](R/gibbs/fast_smoother.R) — disturbance +
  state smoother (port of fast_smoother.m §1+§2, ~95 of the original
  232 lines). MSE recursions (§3) intentionally **NOT** ported; not
  needed by the Gibbs loop. Add later if/when we want uncertainty
  bands on smoothed states.
- [`R/gibbs/sim_smoother.R`](R/gibbs/sim_smoother.R) — Durbin-Koopman
  simulation smoother (port of simulation_smoother.m, 70 lines).

Tests at
[`tests/testthat/gibbs/test-smoother.R`](tests/testthat/gibbs/test-smoother.R)
(9 new, 16 total at this point):
- fast_smoother outputs correct shapes
- simulation_smoother runs and returns valid draws
- **CRITICAL**: average of 500 simulation_smoother draws converges to
  fast_smoother's point estimate within MC tolerance. This validates
  that sim_smoother actually draws from the correct posterior.
- Both handle missing observations without crashing.

### Step 4 — Kim-Shephard SV update ✅ DONE 2026-05-26

Ported as [`R/gibbs/update_vol.R`](R/gibbs/update_vol.R) (158 lines).
Mixture-component constants (probs/means/vars × 10) copied verbatim
from MATLAB — calibrated to log(χ²₁), do not edit.

Tests at
[`tests/testthat/gibbs/test-update-vol.R`](tests/testthat/gibbs/test-update-vol.R)
(4 new, 29 total at this point):
- Runs on simulated SV data
- **CRITICAL**: 500-burn + 2000-draw Gibbs chain recovers truth log-σ
  path with no systematic bias (|log bias| < 0.20) AND positive
  correlation with truth (>0.5).
- Handles missing data.
- Mixture constants sum to 1 and stdvs = √vars (defensive).

### Step 5 — Hyperparameter updates ✅ PARTIAL DONE 2026-05-26

`update_gam.m` → [`R/gibbs/update_gam.R`](R/gibbs/update_gam.R) (26 →
33 lines). Inverse-Gamma conjugate update for SV step sizes.
Vectorised over N independent gammas (Tx N input ⇒ length-N output).

`update_theta.m` (MA(q) coef update) **NOT** ported — Variant A
has no MA errors. Port when/if we add MA(3) noise to match NY Fed.

Tests at
[`tests/testthat/gibbs/test-update-gam.R`](tests/testthat/gibbs/test-update-gam.R)
(4 new, 33 total at this point):
- Runs on vector and matrix input
- **CRITICAL**: posterior concentrates around truth as T → 10K
- Reduces to prior when T = 0 (defensive)

### Step 6 — Time-varying loadings ✅ DONE 2026-05-26

Ported as [`R/gibbs/update_tvcoef.R`](R/gibbs/update_tvcoef.R) (101 →
~120 lines). Carter-Kohn FFBS on the (N*K)-dim vec(alpha_t) RW state.

Tests at
[`tests/testthat/gibbs/test-update-tvcoef.R`](tests/testthat/gibbs/test-update-tvcoef.R)
(3 new, 42 total at this point):
- Runs on small TVP regression
- **CRITICAL**: FFBS draws recover true alpha path with correlation >0.8
- Handles missing observations

Note: our Variant Ac uses CONSTANT loadings, so this isn't strictly
needed for v1 of the Gibbs port. Port done to keep the option open
for adding TVP loadings later.

Original step text (preserved for reference; the rest below is unchanged):

Forward filter + backward sample (Carter-Kohn). ~1 day.

Validation: posterior of latent state from many calls (e.g. 1000) to
sim_smoother on fixed θ should match the marginal posterior derived from
running Stan with `c` and `s` sampled (Variant Ac). Compare medians +
intervals; both should be within MC error.

### Step 4 — Kim-Shephard SV (`update_vol.m` → `R/update_vol.R`)

**The subtle one.** 10-component Gaussian mixture approximation to
log(χ²₁), index variables sampled from full conditional, then the
log-volatility path drawn via Carter-Kohn conditional on the mixture
indices. ~2 days; bugs here are insidious.

Validation: run the SV update on simulated log-volatility data with
known truth. Repeated draws should have correct stationary mean and
variance (i.e. unbiased recovery).

Reference: Kim, Shephard, Chib 1998 §3. The mixture-component table
(weights w_i, means m_i, variances v²_i, i=1..10) is just hardcoded
constants — copy directly from MATLAB.

### Step 5 — Hyperparameter updates (`update_theta.m`, `update_gam.m` → `R/update_hyper.R`)

Conjugate Gibbs draws (Normal-Inverse-Gamma for variances, Normal for
means, Truncated Normal for AR(1) ρ in [-1,1]). ~0.5 days.

Validation: run with known full-conditional means and verify draws have
correct distribution via QQ-plot or Kolmogorov-Smirnov.

### Step 6 — Time-varying loadings (`update_tvcoef.m` → `R/update_lambda.R`)

RW state update for each non-reference sector loading. Same Carter-Kohn
pattern as the latent factor. ~1 day.

Validation: drive with known true lambda paths and check recovery.

### Step 7 — Main Gibbs loop (`estimate_MCT.m` → `R/fit_mct_gibbs.R`)

Orchestrate the per-iteration sweep:
1. Update latent factor c (sim_smoother)
2. Update sector trends s_i (sim_smoother per sector)
3. Update log-volatilities h_c, h_s, h_e (Kim-Shephard SV)
4. Update loadings lambda_i (RW update)
5. Update hyperparameters (conjugate)
6. (If outlier model) Update mixture probabilities + scales

Per-iter cost: dominated by sim_smoother calls. T=435 means ~T scalar
ops per smoother pass × number of smoother calls per iter (~N+3 = 14).
Total: ~6K ops/iter, negligible without autodiff. Many iters (~5K
post-burn) gives ~30M ops, runs in seconds-to-minutes.

~0.5 days.

### Step 8 — Output wrapper (`R/mct_gibbs_fit.R`)

S3 class with `$draws()`, `$diagnostic_summary()`, `$save_object()`.
Also generate the `trend[T]`, `common_share[T]`, `log_lik[n_obs]` as a
post-processing step on the saved draws, matching the Stan generated
quantities exactly. ~0.5 days.

### Step 9 — Adaptation for Variant A (AR(1)) + mixed-frequency obs

NY Fed has Variant B (RW common) and monthly obs only. Two
project-specific tweaks:

- **AR(1) common factor**: replace the random walk c_t = c_{t-1} + ν
  with c_t = ρ c_{t-1} + ν. Affects sim_smoother (state transition
  matrix) and hyperparameter update for ρ (conjugate truncated normal
  given Gaussian errors).
- **Mixed-frequency obs**: quarterly observations are an average of 3
  consecutive monthly latents. Inside the Kalman filter, treat this
  as a missing-data pattern with an aggregated observation operator
  (same logic as `stan/mct_aus_Akf.stan` — see lines on Z^Q for
  quarterly obs).

~1-2 days. The tricky bit: making sure the SV mixture is correctly
parameterised for both monthly and quarterly noise.

### Step 10 — End-to-end validation

- Sim recovery: fit ported Gibbs on simulated data, confirm posterior
  medians fall in 95% intervals of the truth at most t.
- Real-data agreement: fit on real ABS data, compare hyperparameter
  posteriors + trend / common_share to Variant Ac's Stan fit. Should
  match within MC noise (the model is the same; the sampler differs).
- Speed test: full T=435 N=11 fit timed end-to-end. Target ≤ 25 min
  for single chain on M1.

## Adaptation decisions (Australian Variant A specifics)

- **N=11 sectors** vs NY Fed's 17. Data shape only.
- **2000-2019 demean window** vs NY Fed's full-sample demean. Already
  done in `R/data_prep.R`; reuse.
- **Cmdstanr `$save_object()`**: emulate via `saveRDS()`.
- **`format = "file"` target pattern**: keep — `fit_mct_gibbs` returns
  the saved path. Downstream targets `readRDS(fit_*)`.

## Critical differences when porting MATLAB → R

- **1-indexed** (both). No translation needed.
- **Column-major** (both). Matrix mults map cleanly.
- **MATLAB `cell` arrays** → R `list`. Add unit-test cases that
  exercise both heterogeneous and homogeneous shapes.
- **MATLAB `bsxfun`** → R recycling or `outer()`. Watch for shape
  mismatches; MATLAB is more permissive.
- **MATLAB `randn`** → R `rnorm`. Stream order differs; reproducibility
  requires manual seeding of every distribution call.
- **MATLAB `mvnrnd(mu, Sigma)`** → R `MASS::mvrnorm()`. Same arg order,
  different underlying decomposition; expect tiny numerical drift.
- **MATLAB built-in `chol`** → R `chol()`. Both upper triangular by
  default. ✓

## Effort & risk

| Phase | Days | Risk |
|---|---|---|
| 1. Simulator | 1 | Low |
| 2. KF | 1 | Low (can cross-check Akf) |
| 3. Sim smoother | 1 | Med (numerics) |
| 4. Kim-Shephard SV | 2 | **High** (fiddly mixture indexing) |
| 5. Hyperparam updates | 0.5 | Low |
| 6. Loadings update | 1 | Med |
| 7. Main loop | 0.5 | Low |
| 8. Output wrapper | 0.5 | Low |
| 9. AR(1) + mixed-freq adaptation | 1.5 | **High** (changes MATLAB structure) |
| 10. Validation & debug | 2-3 | Med |
| **Total** | **~11-12 days** | (1.5-2 weeks focused) |

**Highest-risk step**: Kim-Shephard mixture (#4). Off-by-one in the
mixture-component table or a wrong covariance computation in the
augmented Gaussian step yields silently biased SV. Test with synthetic
SV data and check posterior coverage. Reference implementation in
`stochvol` R package (Kastner & Frühwirth-Schnatter 2014) is a good
sanity check — fit a simple SV model with both `stochvol::svsample()`
and the ported Kim-Shephard and compare.

## Output format compatibility checklist

Downstream code that consumes a fit object:

- [ ] `R/postprocess.R::extract_path_matrix(fit, param)` — expects
      `fit$draws(variables=param, format="draws_matrix")` returning a
      draws-by-T matrix.
- [ ] `R/postprocess.R::summarise_path(draws, dates)` — pure on matrix
      input, no fit-object dependence.
- [ ] `R/postprocess.R::sector_contributions(fit, weights, dates, groups)`
      — expects `fit$draws("sector_trend")` as draws × T × N.
- [ ] `R/postprocess.R::posterior_revisions(current_fit, prior_fit, ...)`
      — uses `extract_path_matrix` on both.
- [ ] `R/compare_variants.R::loo_for_fit(fit)` — expects
      `fit$draws("log_lik")`. **Action**: confirm log_lik shape matches
      Stan's `vector[n_obs_m + n_obs_q] log_lik`.
- [ ] `R/fit_stan.R::fit_mct()` watchdog reads
      `fit$diagnostic_summary()$ebfmi` — return numeric vector of
      length n_chains; values of 1.0 are a fine no-op.

## Decision triggers (when to actually do this work)

Only start the port when one of these is true:

1. **Production runtime becomes painful**: e.g., we want intra-month
   re-fits, multi-vintage backfills, or sensitivity sweeps.
2. **We extend to N=30-50 sub-indices** (closer to NY Fed's 17 PCE
   sub-indices in spirit). Per the scaling argument in `plan.md`, the
   Stan-side cost goes as N²; the Gibbs cost stays linear. Crossover
   somewhere around N≈20.
3. **The KF approach is needed for a different reason** (e.g., we want
   exact marginal likelihood for model comparison via Bayes factors,
   which Stan's `bridgesampling` can't deliver reliably here).
4. **There's an RA/grad student available** for whom this is a
   self-contained two-week project.

Until one of these triggers, **Ac on the monthly cron is the
production answer** (see `plan.md`).

## Open questions to resolve before starting

- Does the NY Fed code include an AR(1)-common variant, or is it
  strictly RW-common? Inspect `update_gam.m` to find out — if not, add
  the AR(1) update from scratch (conjugate truncated normal given
  Gaussian errors).
- License compatibility: NY Fed Federal Reserve Bank code is typically
  in the public domain. Confirm in `MCT-NYFED/LICENSE`.
- Should Gibbs be a permanent replacement for Stan, or a "speed mode"
  toggled per fit? Suggest the latter for safety — keep Ac as the
  reference until Gibbs has months of production use.

## References

- `plan.md` — context, Phase 1-3 results, why we ended up here.
- `MCT-NYFED/README.md` — original NY Fed README with model spec.
- `stan/mct_aus_Akf.stan` and `mct_aus_Akfs.stan` — already-validated
  KF math, useful as cross-checks for the R port's filter step.
- `tests/testthat/test-stan-simulated.R` — existing simulator + recovery
  tests that the ported sampler should also pass.
