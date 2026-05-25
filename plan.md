# Sampling speedup plan

Goal: cut Stan fit time for Variant A from ~90 min (M1) / ~3 h (x86 runner)
to <30 min, ideally <10 min. The current ~19,650-parameter HMC is far from
the state of the art for this model class.

## Diagnosis

Variant A at T=435, N=11 samples **~19,650 parameters**, of which
**~19,575 (>99.6%) are standard-normal innovations for non-centred
random walks** (`z_*` matrices in
[stan/mct_aus_A.stan:56-61](stan/mct_aus_A.stan#L56-L61)).
The model's "real" degrees of freedom — hyperparameters + log-vol paths
+ loadings — are about **5,650**. HMC is exploring 14k extra dimensions
that exist only because we've written the latent states as parameters
rather than integrating them out.

Empirical evidence the geometry scales badly:
- T=120 (sim, ~5,500 params): ~5 min/fit
- T=435 (real, ~19,650 params): ~93 min/fit
- 3.6× more params → **19× more time** (superlinear)

No within-chain parallelism (`stan_threads=FALSE`, no `reduce_sum`).
No profiling has been done.

## Literature findings

| Source | Takeaway |
|---|---|
| [Sørensen 2026, arXiv:2603.04003](https://arxiv.org/pdf/2603.04003) | NUTS-Kalman on a structurally similar DSEM-with-latent-AR(1) model is **8-19× faster** than brute-force NUTS-Joint (which is what we do). N=100, T=200: **534 min vs 6339 min** to ESS=10K. Open-source Stan reference code at [osf.io/p754m](https://osf.io/p754m). |
| [Arnold's ssmodels-in-stan](https://github.com/jrnold/ssmodels-in-stan) | Canonical hand-rolled KF library for Stan. Filter + smoother + backwards sampler. Use as crib for our implementation. |
| [Kokkala blog](http://www.juhokokkala.fi/blog/posts/kalman-filter-style-recursion-to-marginalize-state-variables-to-speed-up-stan-inference/) | Concrete walk-through of marginalizing latent state via KF recursion inside a Stan `transformed parameters` block. |
| [Stan issue #2854](https://github.com/stan-dev/stan/issues/2854) | `gaussian_dlm_obs` does NOT support time-varying W. Our SV makes W_t time-varying ⇒ we must hand-roll. |
| [Davis et al. 2025, CPT:PSP](https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1002/psp4.13238) | `reduce_sum` within-chain parallelism gives **2-4× per chain** on likelihoods like ours. Compounds with marginalization. |
| [Stan forum 2024 — SV factor model](https://discourse.mc-stan.org/t/stochastic-volatility-for-latent-factor-models/34140) | Bob Carpenter on Stan-side SV+factor models: tighten priors, normalise parameter scales via `<multiplier>`, vectorise loops, replace `diag_matrix()` with elementwise ops. Several easy wins likely apply to our code. |

## Options, ranked by speedup × effort

| # | Approach | Expected speedup | Effort | Compounds with |
|---|---|---|---|---|
| 1 | **KF marginalization** of (c, s) conditional on (h, λ) | **8-15×** | 2-3 days | #2, #3 |
| 2 | Drop time-varying λ (constant loadings) | 1.5-2× | 1-2 h | all |
| 3 | `reduce_sum` within-chain parallelism | 1.5-2× per chain | 2-4 h | #1 |
| 4 | Profile current model (Stan `profile()`) | informs above | 1-2 h | — |
| 5 | Tighter mass-matrix adaptation | 1.2-1.5× | 1 h | all |

Stacked best-case: **#1 + #2 + #3 + #5 ≈ 20-50× speedup** → monthly fit
in 2-5 min on M1, well within the runner's free tier (which would let us
restore the full A+B+C cron build).

## Roadmap

### Phase 1 — cheap wins, no architectural risk (this session)

- [ ] **#2 Constant-λ Stan variant.** Create `stan/mct_aus_A_constlambda.stan`
  as a copy of A with `lambda` as a time-invariant `N-1` vector instead
  of an RW. Add a `variant = "A_const"` branch in `fit_mct()`. Rationale:
  `sigma_lambda ~ N(0, 0.01)` already says loadings barely move; removing
  4,350 z_lambda params is nearly free.
- [ ] **#5 Adapt config tune.** Try `init_buffer = 200, term_buffer = 200,
  window = 50` to give the mass-matrix more time to converge. Same total
  warmup (1500); just better-allocated.
- [ ] **#4 Profile.** Add `profile("rw_construction") {...}`,
  `profile("monthly_lik") {...}`, `profile("quarterly_lik") {...}` to
  `mct_aus_A.stan`. Run one short fit (sim T=120) and read
  `fit$profiles()`. Confirms whether time is in TP-block matrix
  construction or likelihood.
- [ ] **Validation.** After each change: fit on sim data at T=120, confirm
  posterior recovers truth. Then fit on real data at T=435, confirm trend
  + common_share medians match baseline within MC error and LOO ranking
  preserved.

Expected end-of-phase: 2-4× speedup (90 min → 25-45 min) with confidence
nothing has broken.

### Phase 2 — within-chain parallelism (next session)

- [ ] **#3 `reduce_sum` over observations.** Partition `(y_m, y_q)`
  contributions into a `partial_sum_lpdf` and call
  `reduce_sum(partial_sum_lpdf, ...)`. Build with
  `cpp_options = list(stan_threads = TRUE)`, sample with
  `threads_per_chain = 2` (so 4 chains × 2 threads = 8 cores on M1).
- Validation as above.

Expected end-of-phase: +1.5-2× on top of Phase 1, i.e. cumulative 4-8×
(90 min → 11-22 min).

### Phase 3 — Kalman-filter marginalization (multi-day commit)

This is the headline structural change. Decide GO/NO-GO after Phases 1 & 2
based on whether the cheap stack already brought us to acceptable runtime.

Algorithm (per Sørensen §3.2, adapted):

1. **Sampled parameters**: hyperparameters (~50) + log-vol RW innovations
   (`z_hc`, `z_hs`, `z_he` — ~1,300 at T=435) + `lambda` (constant,
   N-1=10 if Phase 1 done; else ~4,350 with RW). Initial condition
   means for c, s.
2. **Marginalised** in the KF: latent paths `c` (AR(1) with SV) and
   `s_{i,t}` (RW with SV) — currently ~14,000 z's, all gone.
3. **State vector for the filter**: augment with 2 lags to handle
   quarterly accumulator: `[c_t, c_{t-1}, c_{t-2}, s_{1,t..t-2},
   ..., s_{N,t..t-2}]` → dim `3(N+1) = 36`.
4. **Time-varying transition**: `T_t` is mostly the identity shift; the
   AR(1) row has `[ρ, 0, 0, ...]`. Process noise `W_t = diag(exp(h_c[t]),
   0, 0, exp(h_s[1,t]), 0, 0, exp(h_s[2,t]), ...)`.
5. **Observation matrix `Z_t`** picks rows for whichever sectors are
   observed at time t. For monthly: rows are `[λ_i, 0, 0, 0_{...}, 1,
   0, 0, 0_{...}]` (current c and current s_i). For quarterly: rows are
   `[λ_i/3, λ_i/3, λ_i/3, 0, ..., 1/3, 1/3, 1/3, 0, ...]` (average over 3
   lags of c and s_i).
6. **Filter recursion** as Sørensen Eq 2.6 (predict + update). Likelihood
   contribution per t as Eq 2.7.
7. **Smoother** in `generated quantities` for posterior recovery of c,
   s, trend, common_share, log_lik.

Implementation steps:
- [ ] Write `stan/mct_aus_A_kalman.stan` as new variant. Crib state-space
  primitives from Arnold's `ssmodels-in-stan`.
- [ ] Validate against current A on sim data (T=120, N=5): both should
  produce same posterior medians + intervals on truth.
- [ ] Validate on real data (T=435): trend, common_share, log_lik should
  match A baseline within MC noise. LOO should be within ±2 SE of A's
  baseline elpd.
- [ ] If it works: deprecate the old A model; update cron back to full
  A+B+C (since C's 108 min becomes ~10 min).

Risks:
- KF numerical stability on T=435 with very small SV-derived process
  noise. Mitigate with Joseph form for covariance update.
- Quarterly accumulator + lag augmentation is error-prone. Validate
  the obs equation arithmetically on a 6-month sim before scaling.
- 36-dim state covariance updates are O(36³) per step → 47K ops/step ×
  435 steps = 20M ops per filter pass. Should be ~10× cheaper than the
  current TP-block matrix reconstructions. Worth profiling at the end.

### Phase 4 — fallback (only if Phases 1-3 still too slow)

The deferred NY Fed Gibbs port via Kim-Shephard mixture is the
gold-standard fast SV approach (10-25 min single chain on T=435 per the
handoff). Estimated 1-2 weeks of R porting from the MATLAB reference at
`MCT-NYFED/functions/`. Park unless the Stan stack genuinely can't get
under, say, 15 min/fit after all of Phases 1-3.

## Session log — 2026-05-25

**Files created this session:**
- `stan/mct_aus_Ac.stan` — Phase 1 #2 (constant lambda). Compiles in ~10s.
- `stan/mct_aus_Act.stan` — Phase 2 #3 (Ac + `reduce_sum`). Compiles in ~24s
  (slower due to threading machinery).
- `stan/mct_aus_Akf.stan` — Phase 3 #1 **DRAFT**. Compiles in ~24s.
  Smoother in `generated quantities` is a TODO placeholder — filter-only
  works for measuring sampling speed, but trend/common_share extraction
  needs RTS smoother written before this can replace Ac downstream.
- `scripts/benchmark_Ac.R` — head-to-head Ac vs cached A.
- `reference/sorensen/dsem_kalman_non_centered.stan` + `compile_models.R`
  — Sørensen 2026 reference code for the NUTS-Kalman pattern.

**Mid-warmup Ac measurement (T=435 N=11 real data):**
- 600 iters/chain in ~15 min → **~40 iter/min/chain** (vs A's ~32)
- Implied total: ~75 min vs A's 93 min → **~20% speedup from #2 alone**
- Less than the 1.5-2× hoped; consistent with removing 22% of params
  saving ~22% of gradient work (linear, not multiplicative).
- Final wall-time + validation pending (benchmark still running).

**Profile() blocks are expensive — banned from production runs.**
First Ac benchmark with `profile("...") { ... }` blocks ran 5× slower
than baseline A. Stan's profile runs per leapfrog step, not per iter;
with max_treedepth=12 the overhead compounds. Re-add only for one-off
diagnostic runs.

**Akf math notes (to validate later):**
- Augmented state dim D = 3(N+1) = 36 at N=11.
- Indexing: x[1..3] = (c_t, c_{t-1}, c_{t-2});
  x[3i+1..3i+3] = (s_{i,t}, s_{i,t-1}, s_{i,t-2}) for i=1..N.
- Transition matrix sparse (~24 non-zeros at N=11); current draft uses
  dense `quad_form_sym(P, Tm')`. Optimisation hook: replace with
  cell-by-cell updates to save O(D^3) → O(D) per step.
- Process noise diagonal: exp(h_c[t]) at [1,1]; exp(h_s[i,t]) at
  [3i+1, 3i+1].
- Monthly obs equation: Z row has lambda[i] at col 1, 1 at col 3i+1.
- Quarterly obs equation: Z row has lambda[i]/3 at cols 1..3,
  1/3 at cols 3i+1..3i+3.
- Initial conditions: c_1 ~ N(0, exp(h_c[1])/(1-rho^2)); lags c_0, c_{-1}
  diffuse (var 100); s_{i,1} ~ N(s_init[i], 0.25), lags diffuse.
- Update step: non-Joseph form (cheaper but less stable). Switch to
  Joseph if R/F instability shows up in validation.
- **PRE-CONDITION**: input obs arrays must be sorted by t ascending.
  Currently NOT enforced in build_stan_data(). Add a sort there before
  calling fit_mct(variant="Akf").
- Smoother TODO: Rauch-Tung-Striebel backward pass to recover smoothed
  (c, s) — needed for trend/common_share/log_lik in generated quantities.

## Session log — 2026-05-26 Phase 3 results (T=120 sim)

**Phase 3 KF marginalization, dense + sparse implementations both validated:**

| Variant | Wall (min) | vs Ac | vs Akf | Divergences | ebfmi |
|---------|-----------|-------|--------|-------------|-------|
| Ac      | 5.28      | 1×    | 33×    | 20          | 0.79  |
| Akf     | 174.83    | 0.03× | 1×     | **0**       | 0.92  |
| Akfs    | 40.95     | 0.13× | **4.27×** | **0**    | 0.945 |

**Math: validated.** Akf and Akfs hyperparameter posteriors match Ac to 3 decimal places across rho/mu_hc/sigma_hc/lambda_const. Zero divergences on both KF variants vs Ac's 20.

**Sparse implementation (mct_aus_Akfs.stan) wins ~4× over dense Akf** by exploiting that the D=36 transition matrix has only 24 non-zeros (1 per row). Computing `T·P·T'` cell-by-cell drops the predict step from O(D³)=47K ops to O(D²)=1296 ops. Obs update similarly: Z·P·Z', P·Z', and (I−K·Z)·P all done with explicit lookups, never materialising Z.

**Scaling argument: Akfs is the wrong architecture for production scale.**
- Per-step cost O(T·D²) where D=3(N+1).
- T=120 N=5: 41 min.
- T=435 N=11: ~41 × 3.6(T-ratio) × 4(D-ratio²) ≈ **570 min/fit**.
- Versus Ac at production: ~93 min.
- The autodiff overhead on matrix-shaped state dominates the savings from removing latent-state innovations from the sampler.

**Implication.** Sørensen's 8-19× speedup happens because his scalar latent KF has 5 ops/step. Our multi-sector setup forces ~D² ops/step which Stan's autodiff can't compete with the direct sampling of Ac. The KF strategy only beats sampling when the state is small (≤3) or the data is much wider (multiple subjects/series).

**Decision point.** Three remaining paths:

- **A. Further Akfs optimization** (block-P, cumulator state, alloc reduction): 2-4 days for maybe 2-3× more. Still loses to Ac at production scale on this hardware. Marginal value.
- **B. Phase 4 Gibbs port** from NY Fed MATLAB to R: ~1 week of careful porting. 10-25 min/fit per NY Fed timings — the only credible path to a transformative speedup. Loses Stan ecosystem.
- **C. Accept Ac (or A) as production**: 0 effort. 93 min/fit on cron, monthly automation already works.

Current recommendation: **C, with the Phase 3 work documented** as a tried-and-found-not-helpful branch. If runtime ever becomes critical (e.g., we want intra-month re-fits, real-time updates, or multi-vintage backfills), revisit B.

## Session log — 2026-05-25 head-to-head results (T=120 sim)

Clean apples-to-apples comparison: same simulated data, same 4 chains ×
750+750 iters, same seed, same adapt_delta=0.97.

| Variant      | Wall (min) | vs A    | Divergences | ebfmi min |
|--------------|------------|---------|-------------|-----------|
| A (baseline) | 6.02       | 1.00×   | 4           | 0.71      |
| Ac (const λ) | 5.10       | 1.18×   | 20          | 0.79      |
| Act (2 thr)  | 30.76      | 0.20×   | 12          | 0.69      |

**Phase 1 #2 verdict (const λ).** Real but marginal: 1.18× speedup with
5× more divergences. Lambda RW is doing slightly more identification work
than the prior `sigma_lambda ~ N(0, 0.01)` suggested.

**Phase 2 #3 verdict (reduce_sum).** Net 5× *regression* on this M1
with 4 chains × 2 threads = 8 threads / 8 cores. Pure OS oversubscription.
Could be revisited with 2-chain × 4-thread or 1-chain × 8-thread configs,
but those sacrifice the 4-chain mixing diagnostics the project wants.
**Skip on this hardware unless we drop to 2 chains.**

**Conclusion.** The cheap-stack ceiling is ~18%. The only credible path
to a large speedup is **Phase 3 KF marginalization**. Move there next.

## Validation discipline (apply to every phase)

Before declaring any speedup variant "done":

1. Sim recovery: fit on simulated data with known truth at T=120, N=5.
   Posterior medians fall in 95% intervals of the truth for the headline
   `trend` series at most t.
2. Real-data match: fit on T=435 real ABS data. Compare to baseline:
   - Latest-month trend median: must be within ±0.05 (current: 0.65).
   - Common-share latest median: within ±0.02 (current: 0.074).
   - LOO elpd: within ±2 SE of baseline.
3. Diagnostics: no divergences in >0.5% of post-warmup transitions;
   ebfmi ≥ 0.3 all chains; R̂ ≤ 1.01 on all monitored params.

If any validation fails, revert and diagnose before moving on.

## Code sketch (Phase 2 — `reduce_sum` over observations)

Once Phase 1 lands, drop-in for `mct_aus_Ac.stan`. Two parallel partial
sums (monthly and quarterly likelihoods are conditionally independent
given the latent paths). Build with `cpp_options = list(stan_threads = TRUE)`,
sample with `threads_per_chain = 2`.

```stan
functions {
  // Sum monthly log-lik over a slice of indices [start..end].
  // All time-varying paths are accessed by reference via the slice's
  // index arguments — Stan handles the slicing automatically.
  real partial_sum_monthly_lpdf(
      array[] int slice_idx, int start, int end,
      vector y_m, array[] int t_m, array[] int i_m,
      vector lambda, vector c, matrix s_trend, matrix h_e) {
    real lp = 0;
    for (k in start:end) {
      int t = t_m[k]; int i = i_m[k];
      real mu = lambda[i] * c[t] + s_trend[t, i];
      lp += normal_lpdf(y_m[k] | mu, exp(h_e[t, i] / 2));
    }
    return lp;
  }
  // Same shape for quarterly — three-month accumulator inside the loop.
  real partial_sum_quarterly_lpdf(...) { ... }
}
model {
  // ... priors and TP-block as before ...
  int grainsize = 1;  // let Stan tune
  array[n_obs_m] int idx_m;
  for (k in 1:n_obs_m) idx_m[k] = k;
  target += reduce_sum(partial_sum_monthly_lpdf, idx_m, grainsize,
                       y_m, t_m, i_m, lambda, c, s_trend, h_e);
  target += reduce_sum(partial_sum_quarterly_lpdf, idx_q, grainsize,
                       y_q, t_q, i_q, lambda, c, s_trend, h_e);
}
```

Expected gain at 2 threads/chain on M1: ~1.5-1.8× per chain (Amdahl-limited
because TP-block reconstruction of the SV paths stays serial). Sublinear
gain is fine — it stacks with #2 and #1.

## Sørensen reference code (cached at `reference/sorensen/`)

Downloaded `dsem_kalman_non_centered.stan` and `compile_models.R` from
Sørensen's `eg2-three-indicator-ar1` OSF folder. 111 lines, scalar
latent AR(1) factor with 3 observed indicators — structurally the
nearest published analogue to our setup.

**Key patterns to crib for our Phase 3 KF:**

1. **Pre-compute observation-side summaries outside the time loop.**
   He builds `pre_q[t]` and `pre_r[t]` (linear and quadratic functions of
   `y[t] - nu` weighted by inverse obs variance) once before the KF loop.
   Inside the loop the per-step KF update is then just **5 scalar ops**:

   ```stan
   real qt = pre_q[t] - a * c;
   real rt = pre_r[t] - 2.0 * a * pre_q[t] + a * a * c;
   real F_det_term = 1.0 + P * c;
   real P_over_F = P / F_det_term;
   sum_ll_dynamic += log(F_det_term) + rt - P_over_F * qt * qt;
   real a_post = a + P_over_F * qt;
   a = phi * a_post;
   P = phi_sq * P_over_F + sigma_proc_sq;
   ```

2. **Jungbacker-Koopman compression** of the U-dim observation onto a
   V_1-dim compressed space. He has V_1=1 (scalar latent) so the
   3-indicator obs collapses to scalars `qt`, `rt`. **Doesn't trivially
   apply to us** because we also have sector-specific trends `s_{i,t}`
   that need to be part of the state.

3. **Return `-0.5 * (T*N_obs*log(2pi) + T*log_det_Sigma + sum_ll_dynamic)`**
   — collects all the constants so the per-step loop only accumulates
   the data-dependent quadratic form.

4. **Non-centered parameterisation only for the BETWEEN-LEVEL random
   effects** (subject-specific intercepts etc.). The within-level latent
   state is fully marginalised — never appears as a parameter.

**What's different in our problem:**

- SV: `sigma_obs` → `exp(h_e[i,t]/2)`, `sigma_proc` → `exp(h_c[t]/2)`.
  Time-varying. The pre-loop precomputation has to move inside the loop
  for the parts that depend on h. Mitigation: pre-compute the
  `exp(h/2)` outside the KF (since h is fixed during a leapfrog step).
- Sector-specific trends `s_{i,t}`: extra latent random walks. Must
  go in the state vector (or be marginalised alongside `c`).
- Mixed frequency: quarterly obs is an average of 3 latent monthly
  values. Either augment state with 2 lags, OR process quarterly
  obs as a special multi-step likelihood contribution.

**Recommended state for our KF:** dim N+1 = 12 (joint
`(c_t, s_{1,t}, ..., s_{N,t})`). Augment with 2 lags only when there's
a quarterly obs at time t (handled by computing the obs likelihood
contribution against a virtual 3-step-averaged observation operator,
not by carrying lags in the state perpetually). Per-step KF cost ≈ 12³
= 1,728 ops × T=435 = 750k ops per filter pass. ~3 orders of magnitude
cheaper than the current TP-block reconstruction cost.

## Code sketch (Phase 3 — KF in Stan)

```stan
functions {
  // Returns the log marginal likelihood after running the Kalman filter
  // over T steps with time-varying W_t. State dim = 3*(N+1) = 36.
  real kalman_lpdf(/* args: state init, h_c, h_s, h_e, lambda, ...,
                     y_m, t_m, i_m, y_q, t_q, i_q, T, N */) {
    int D = 3 * (N + 1);
    vector[D] a = a_init;
    matrix[D, D] P = P_init;
    real loglik = 0;
    for (t in 1:T) {
      // Predict step
      vector[D] a_pred = T_mat * a;
      matrix[D, D] W_t = build_W_t(h_c[t], h_s[t]);   // time-varying!
      matrix[D, D] P_pred = quad_form_sym(P, T_mat')
                            + W_t;
      // Find observations at time t (both monthly and quarterly)
      // For each obs: build Z_row, compute innovation, update
      // (use Joseph form for numerical stability)
      a = a_pred;
      P = P_pred;
      for (obs in obs_at_t) {
        row_vector[D] Z = build_Z(obs, lambda);
        real H = obs_noise_var(obs, h_e[t]);
        real v = obs_y[obs] - Z * a;
        real F = Z * P * Z' + H;
        loglik += -0.5 * (log(2*pi()) + log(F) + v*v/F);
        // Joseph form update
        vector[D] K = P * Z' / F;
        a += K * v;
        P = (identity_matrix(D) - K * Z) * P * (identity_matrix(D) - K*Z)'
            + K * H * K';
      }
    }
    return loglik;
  }
}
model {
  // Priors on hyperparams and h_* and lambda as before
  target += kalman_lpdf(...);  // marginal lik
}
generated quantities {
  // Run Kalman smoother backwards to recover posterior of c, s
  // then derive trend, common_share, log_lik
}
```

## Open questions to revisit later

- Time-varying λ: does dropping it bias common-share estimates? Phase 1
  validation against baseline will answer.
- Whether to use Stan's profile() output to drive priors or just sampler
  config — depends on where time actually goes.
- Whether to lift Phase 4 (Gibbs port) into a separate `mct_aus_gibbs/`
  submodule rather than fighting Stan further.

## References (full list)

- Sørensen (2026), *Efficient Bayesian Estimation of DSEMs via State
  Space Marginalization*, arXiv:2603.04003 — the lit anchor.
- Arnold's [ssmodels-in-stan](https://github.com/jrnold/ssmodels-in-stan)
  — Stan KF/smoother code to crib.
- Kokkala blog on KF marginalization in Stan — concrete walk-through.
- Davis et al. (2025), CPT:PSP — within-chain parallelism gains.
- Stan SUG §2.5 — Stochastic Volatility Models reference.
- Stan SUG §3.4 — Parallelization.
- Stan issue [#2854](https://github.com/stan-dev/stan/issues/2854) —
  why we can't use `gaussian_dlm_obs`.
- Stan forum thread on [SV in latent factor models](https://discourse.mc-stan.org/t/stochastic-volatility-for-latent-factor-models/34140)
  — Carpenter's debugging advice.
