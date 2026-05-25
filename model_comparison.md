# MCT model spec: ours vs NY Fed

Side-by-side of the Australian MCT model implemented in this repo
(headline variant: A; alternates B and C documented below) against the
NY Fed Multivariate Core Trend model implemented in
[`MCT-NYFED/`](MCT-NYFED/) (read-only reference).

Anchors:
- Our spec: [`CLAUDE.md`](CLAUDE.md), [`stan/mct_aus_A.stan`](stan/mct_aus_A.stan)
- NY Fed spec: [`MCT-NYFED/automated_PCE.m`](MCT-NYFED/automated_PCE.m),
  [`MCT-NYFED/functions/estimate_MCT.m`](MCT-NYFED/functions/estimate_MCT.m)
- Stock & Watson (2016), *Core Inflation and Trend Inflation*, REStat 98(4)

## At a glance

| Feature                            | Ours (Variant A)                    | NY Fed MCT                          |
|------------------------------------|-------------------------------------|-------------------------------------|
| Country / data source              | Australia, ABS CPI 6401.0           | US, BEA PCE                         |
| Sectors (N)                        | 11 ABS CPI groups                   | 17 PCE sub-indices                  |
| Time span (T)                      | 1990m1 → present (~435 months)      | 1960m1 → present (~750+ months)     |
| Observation frequency              | **Mixed**: quarterly + monthly (since Apr-2024) | Monthly only                |
| Common factor                      | **AR(1) transitory cycle** (Variant A) | **RW common trend + i.i.d. cycle** |
| Sector trends                      | RW with SV                          | RW with SV                          |
| Time-varying loadings              | Yes, RW (`σ_λ ~ N(0, 0.01)`)        | Yes, RW (`s²_λ = 0.25²/60/12`)      |
| Stochastic volatility              | On common, sector trends, meas. noise | On common trend, common cycle, sector trends, meas. noise |
| MA error terms                     | **No** (Gaussian iid)               | **Yes, MA(3)** by default           |
| Outlier / fat-tail treatment       | **No**                              | **Yes** (scale mixture, slice-sampled indicators) |
| Cross-sector dependence            | No                                  | Optional (`i_depend` flag)          |
| Identification                     | `λ[Housing, t] = 1 ∀t`              | Normalised post-estimation (sign + scale) |
| Estimator                          | **HMC via Stan / cmdstanr**         | **Gibbs (Kim-Shephard mixture for SV)** |
| Single-fit wall time (production)  | ~93 min (M1, 4 chains × 3000 iters) | ~10–25 min (single chain × 3000 iters) |

## Measurement equation

**Ours (Variant A)** — see [`stan/mct_aus_A.stan:148-153`](stan/mct_aus_A.stan#L148-L153):

$$
\pi_{i,t} = \lambda_{i,t}\, c_t + s_{i,t} + \exp(h^\varepsilon_{i,t}/2)\,\varepsilon_{i,t},\quad \varepsilon_{i,t}\sim N(0, 1).
$$

The common factor `c_t` is a stationary AR(1) transitory cycle, not a
trend. Trend inflation in Variant A is recovered as

$$
\tau_t = \sum_i w_i\, s_{i,t}.
$$

**NY Fed** — see comments in
[`MCT-NYFED/functions/estimate_MCT.m`](MCT-NYFED/functions/estimate_MCT.m):

$$
\pi_{i,t} = \alpha_{\tau,i,t}\,\tau^c_t + \tau_i^t + \alpha_{\varepsilon,i,t}\,\varepsilon^c_t + e_{i,t},\quad
e_{i,t} \sim \text{MA}(q),\; q = 3.
$$

with `e_{i,t}` having stochastic volatility AND an outlier scale mixture
(see "Outlier treatment" below). Trend inflation in their model is

$$
\tau_t = \big(\sum_i w_i\, \alpha_{\tau,i,t}\big)\cdot \tau^c_t + \sum_i w_i\, \tau_i^t.
$$

The structural difference: their measurement equation has BOTH a
common trend `τ^c_t` (with time-varying loadings `α_τ`) and a common
cycle `ε^c_t` (with separate loadings `α_ε`). We chose between them
(Variant A picks AR(1) cycle; Variant B picks RW trend; Variant C
combines both but lost the 3-way LOO comparison decisively — see
[`README.md`](README.md)).

## Common-factor structure

| Component                  | Ours (A) | Ours (B) | Ours (C) | NY Fed       |
|----------------------------|----------|----------|----------|--------------|
| RW common trend $\tau^c$   | ✗        | ✓        | ✓        | ✓            |
| AR(1) common cycle $c$     | ✓ (stat.)| ✗        | ✗        | ✗            |
| i.i.d. common cycle $\varepsilon^c$ | ✗ | ✗      | ✓        | ✓            |
| Time-varying loadings on common trend | n/a | ✗ (just `c`) | ✗ (time-invariant) | ✓ |
| Time-varying loadings on common cycle | ✓ | n/a   | ✗ (time-invariant) | ✓ |

**Key gap:** NY Fed identifies both a common trend AND a common cycle
*simultaneously*, with separate time-varying loadings on each. Our
Variant C tried this but the LOO comparison strongly preferred the
simpler Variant A (elpd_diff = −71.6, SE = 11.2 — a 6σ rejection of C
vs A). The empirical Australian story is overwhelmingly transitory
co-movement, not a drifting common trend.

## State equations and SV

Both models use random walks with stochastic volatility for the latent
states. The SV structure is:

**Ours** — [`stan/mct_aus_A.stan:72-92`](stan/mct_aus_A.stan#L72-L92):

$$
\begin{aligned}
\lambda_{i,t} &= \lambda_{i,t-1} + \sigma_\lambda\, \eta^\lambda_{i,t},\\
s_{i,t}      &= s_{i,t-1} + \exp(h^s_{i,t}/2)\,\eta^s_{i,t},\\
h^c_t        &= h^c_{t-1} + \sigma_{h_c}\,\eta^{h_c}_t,\\
h^s_{i,t}    &= h^s_{i,t-1} + \sigma_{h_s,i}\,\eta^{h_s}_{i,t},\\
h^\varepsilon_{i,t} &= h^\varepsilon_{i,t-1} + \sigma_{h_\varepsilon,i}\,\eta^{h_\varepsilon}_{i,t},
\end{aligned}
$$

with non-centred parameterisation throughout. Initial conditions: AR(1)
common factor uses the stationary distribution; sector trends and
loadings have weak priors near 0/1 respectively.

**NY Fed** — see comments at the top of
[`MCT-NYFED/functions/update_vol.m`](MCT-NYFED/functions/update_vol.m).
Same RW-with-SV structure, but on FIVE separate volatility paths:
`σ_dτ_c` (common trend innovation), `σ_dτ_i` (sector trend
innovation), `σ_ε_c` (common cycle), `σ_ε_i` (sector noise), and
implicit `σ_λ` (loading drift). The SV update uses the
[Omori et al. (2007)] 10-component mixture approximation to the log-χ²
distribution.

## Outlier / fat-tail treatment

**Ours**: **none**. Sectors with idiosyncratic shocks (fuel, fruit &
veg, travel) get absorbed entirely through the SV on
`h^\varepsilon_{i,t}`. The [`next_session.md`](next_session.md) v2
backlog flags this as a sensitivity to run with winsorised series.

**NY Fed**: explicit **scale mixture for outliers** on both common and
sector-specific transitory shocks. See
[`MCT-NYFED/functions/update_scl.m`](MCT-NYFED/functions/update_scl.m)
and [`MCT-NYFED/functions/update_ps.m`](MCT-NYFED/functions/update_ps.m).
Each observation gets an indicator `s_{ε,c,t} ∈ {1, 2, ..., 10}` whose
value scales the innovation variance. `s = 1` is "normal-noise mode";
`s > 1` flags the observation as an outlier. The indicator is sampled
via the slice sampler from
[`MCT-NYFED/functions/slice_sampling.m`](MCT-NYFED/functions/slice_sampling.m).

This was on our backlog ([`next_session.md`](next_session.md) "NY Fed
item 3") but was deferred — HMC can't mix the discrete indicators
efficiently, so this is one of the things that motivates the Gibbs port
plan in [`gibbs_port_plan.md`](gibbs_port_plan.md).

## MA error terms

**Ours**: Gaussian iid measurement noise per sector
(`exp(h^ε_{i,t}/2) ε_{i,t}` with `ε ~ N(0, 1)`).

**NY Fed**: **MA(3) by default** on the measurement error
(`settings.n_lags = repmat(3, n, 1)` in automated_PCE.m). The MA
parameters `θ` are updated via
[`MCT-NYFED/functions/update_theta.m`](MCT-NYFED/functions/update_theta.m).

This was [`next_session.md`](next_session.md) "NY Fed item 2", also
deferred for the same reason as outliers — HMC mixes the MA
parameterisation poorly. The Gibbs sampler handles it natively via
conditional Gaussian updates given the latent factor path.

## Cross-sector dependence

**Ours**: zero by construction. Each sector's idiosyncratic component
`s_{i,t} + \exp(h^\varepsilon_{i,t}/2)\,\varepsilon_{i,t}` is
independent across `i`.

**NY Fed**: optional. The `i_depend` setting at
[`MCT-NYFED/automated_PCE.m:108`](MCT-NYFED/automated_PCE.m#L108) is
set to all-zeros (`i_depend = zeros(n, 1)`) by default, so this defaults
to OFF in their production runs too. The infrastructure exists in
[`MCT-NYFED/functions/update_tvcoef.m`](MCT-NYFED/functions/update_tvcoef.m)
for sectors to depend on another sector's lagged values via a `β`
coefficient.

## Mixed-frequency observations

**Ours**: **load-bearing.** ABS published a complete monthly CPI only
from October 2025 (back-series to April 2024). The model handles this
via a single state-space spec on a monthly latent clock:

- Monthly obs (Apr-2024 onward): `y_{i,t} = π_{i,t}` directly.
- Quarterly obs (1990–Mar-2024, placed at the last month of each
  quarter): `y_{i,t}^Q = (π_{i,t-2} + π_{i,t-1} + π_{i,t})/3` with
  measurement-noise variance averaged accordingly.

See [`stan/mct_aus_A.stan:155-172`](stan/mct_aus_A.stan#L155-L172).

**NY Fed**: monthly only. They have a `is_timeag` flag in
[`MCT-NYFED/automated_PCE.m`](MCT-NYFED/automated_PCE.m) that could
handle a series whose observations are time-averages of 12 monthly
latents (e.g., an annual series), but in their default run this is set
to false for every series.

## Identification

| Constraint                  | Ours                                  | NY Fed                                  |
|-----------------------------|---------------------------------------|-----------------------------------------|
| Common-factor sign + scale  | `λ_{Housing, t} = 1 ∀t` (hard pin)    | Post-estimation normalisation by `σ_dτ_c[1]` |
| Initial common-factor state | Stationary AR(1) (Variant A); diffuse `N(0, 2²)` (Variant B) | Diffuse `Σ_1 = 0` for trend block |
| Initial sector trends       | `s_{i, 1} ~ N(0, 0.5²)`              | Initial means set from data scale       |
| Initial loadings            | `λ_{i, 1} ~ N(1, 0.3²)` for non-ref  | Initialised at constant `y_scale / 16`  |

The NY Fed approach (no pinned sector) requires post-estimation
normalisation to break sign and scale indeterminacies; we pre-empt
those by pinning Housing's loading to 1.

## Priors

Both models use weakly informative priors but with different shapes
because Gibbs needs conjugate forms.

**Ours** — [`stan/mct_aus_A.stan:115-134`](stan/mct_aus_A.stan#L115-L134):

| Param                | Prior            | Comment                                   |
|----------------------|------------------|-------------------------------------------|
| `μ_h*`               | `N(0, 1)`        | Log-vol RW level                          |
| `σ_h*`               | `N(0, 0.1)`      | Tightened from 0.2 to prevent mode-trap   |
| `σ_λ`                | `N(0, 0.01)`     | Loading drift kept small                  |
| `ρ` (AR(1))          | `N(0.5, 0.3)`    | Truncated to (−1, 1)                      |
| `s_init`             | `N(0, 0.5)`      | Initial sector trend                      |
| `λ_init` (non-ref)   | `N(1, 0.3)`      | Initial loading                           |

**NY Fed** — [`MCT-NYFED/automated_PCE.m:80-95`](MCT-NYFED/automated_PCE.m#L80-L95):

| Param            | Prior                                          | Comment                            |
|------------------|------------------------------------------------|------------------------------------|
| `s²_λ` (loading var) | Inverse-gamma, `ν=12, scale=0.25²/(60·12)` | Tighter than ours                  |
| `s²_γ` (vol-RW var)  | Inverse-gamma, `ν=60, scale=1/(60·12)`     | Strongly anchored                  |
| `α_τ`, `α_ε`         | Conjugate Normal, `prec_beta = 1`           |                                    |
| `θ` (MA)            | Conjugate Normal, `prec_MA = 0.1`           |                                    |
| `p_s` (outlier mix)  | Beta, `(p̄=1−1/(48), n_obs=120)`           | Strong prior that outliers are rare |

NY Fed's tighter priors are a Gibbs requirement: a weakly informative
prior with no conjugate form forces a Metropolis-within-Gibbs step that
mixes badly.

## Estimation method

**Ours**: HMC via Stan / cmdstanr. 4 chains × 1500 warmup + 1500
sampling, `adapt_delta = 0.97`, `max_treedepth = 12`. The non-centred
RW parameterisation exposes ~19,650 standard-normal innovations as
sampled parameters at T=435 N=11; HMC traverses this 19k-dim space
each iteration. Total wall time: ~93 min per fit on M1.

**NY Fed**: Gibbs sampler. 3000 burn + 3000 post-burn, thinned by 2.
Each iteration does:
1. Sample latent factor `τ^c` via Carter-Kohn forward-filter
   backward-sample (FFBS).
2. Sample sector trends `τ_i` via FFBS per sector.
3. Sample SV paths via the Omori et al. 2007 mixture approximation:
   draw the mixture indicators from their full conditional, then the
   log-volatility path via FFBS conditional on the indicators.
4. Sample time-varying loadings `α_τ`, `α_ε` via conjugate Normal updates.
5. Sample MA parameters `θ` via conjugate Normal updates.
6. Sample outlier scale mixture indicators via slice sampling.
7. Sample hyperparameters (variances, AR coefficients) via Inverse-Gamma /
   truncated Normal conjugate updates.

No autodiff, no leapfrog steps, just direct sampling from full
conditionals. Single chain runs in 10–25 min on a comparable machine
(per the README and the Sørensen 2026 paper's analogous timings).

## Net implications — what would we GAIN by porting NY Fed?

1. **~5–10× speedup per fit** (10–25 min vs 93 min) — see
   [`gibbs_port_plan.md`](gibbs_port_plan.md) for the path.
2. **Outlier handling** would let us properly model fuel / fruit & veg
   without winsorising — a methodologically cleaner story.
3. **MA(3) measurement errors** would absorb short-run serial
   correlation that's currently hitting the SV process.
4. **Combined trend + cycle** common factor specification (NY Fed's
   default) — though our Variant C tried this and it lost to A on LOO,
   the test may have been unfair because C didn't have NY Fed's MA(3)
   + outlier scaffolding to absorb noise into the right places.

## Net implications — what we GAIN by sticking with Stan

1. **Mixed-frequency obs**: NY Fed has no public infrastructure for
   the AU quarterly-then-monthly setup. We'd have to add it.
2. **Reproducibility**: Stan's R-hat / ess / divergence diagnostics
   are battle-tested. Gibbs samplers need bespoke convergence
   diagnostics.
3. **Flexibility**: trying a new prior or a new factor structure is
   one line of Stan code vs derivation of a new full conditional.
4. **Maintenance**: the project is already on Stan; introducing a
   second sampler doubles the surface area for bugs and tests.

## What would happen if we ran the NY Fed model on AU data unchanged?

It wouldn't run — it expects monthly observations every period.
Adapting it to AU requires:
1. Modifying the Kalman filter to accept time-averaged quarterly obs
   (the `is_timeag` flag handles 12-month averages, not 3-month).
2. Adjusting `n` from 17 to 11.
3. Adjusting the demeaning window from full-sample to 2000–2019.
4. Translating the MATLAB-only data pipeline to read ABS CSVs/XLSXs.

All addressed in [`gibbs_port_plan.md`](gibbs_port_plan.md) step 9.

## Bottom line

**Our model is a *strict-subset specification* of the NY Fed model**
(missing MA(3), outliers, cross-dependence, combined trend+cycle)
**with the addition of mixed-frequency observations**. The NY Fed
model is more flexible but requires a different sampler. We have the
mixed-frequency story right; the NY Fed has the noise model right.

A v2 effort would: (a) port the NY Fed Gibbs sampler to R, (b) add the
mixed-frequency observation handling, (c) keep our identification
choice (Housing pinned) for simpler downstream interpretation. End
state: a model that does NY Fed's noise modelling, our mixed-frequency
handling, and runs in 10–25 min per fit. Path described in
[`gibbs_port_plan.md`](gibbs_port_plan.md); decision triggers spelled
out there.
