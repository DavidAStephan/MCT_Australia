# Australian Multivariate Core Trend Inflation

A Bayesian mixed-frequency dynamic factor model for Australian CPI, in the
spirit of the NY Fed's [Multivariate Core Trend (MCT)](https://www.newyorkfed.org/research/policy/mct)
inflation model (Stock & Watson 2016). Three variants of the common factor
are estimated and compared:

- **Variant A** — AR(1) transitory common factor; trend = weighted sum of
  sector-specific RW drifts.
- **Variant B** — random-walk common trend (canonical Stock-Watson);
  trend includes the common drift plus sector-specific drifts.
- **Variant C** — combined: BOTH a random-walk common trend AND an i.i.d.
  common cycle, each with stochastic volatility, with time-invariant
  loadings. Nests A (set trend factor to zero) and B (set cycle factor to
  zero) as restrictions, letting the data partition the common component.

See [`CLAUDE.md`](CLAUDE.md) for the full project brief and model specification.

## Status

Three Stan variants fit on the full T=435 sample (1990-Jan → 2026-Mar);
mixed-frequency observations (monthly from Apr-2024).

### Headline finding (2026-05-23)

**Variant A is preferred by LOO** in a 3-way comparison:

```
  elpd_diff  se_diff
A    0.0      0.0    ← winner
B  -13.0      4.4    (A > B by ~3σ)
C  -71.6     11.2    (A > C by ~6σ)
```

Substantively: Australian sectoral CPI co-movement is **predominantly
transitory**, not a drifting common trend. Variant C identifies BOTH a
common cycle and a common drift (the cycle loadings are all clearly
non-zero across sectors, and the cycle's log-vol RW step `sigma_hec=0.17`
is ~3× larger than the trend's `sigma_hc=0.065`), but the combined-factor
structure over-parameterises and loses predictive accuracy versus the
AR(1)-cycle-only A specification.

This is the **opposite** of the NY Fed's finding for US PCE, where the
Stock-Watson random-walk common trend is the canonical decomposition. For
the Australian group-level CPI panel, the data prefers an AR(1) cycle
interpretation of the common factor.

The dashboard headlines Variant A; Variant B and C are shown in the
variant-comparison panel.

## Requirements

- R ≥ 4.3
- [CmdStan](https://mc-stan.org/docs/cmdstan-guide/cmdstan-installation.html) ≥ 2.34 (installed via `cmdstanr::install_cmdstan()`)
- [Quarto](https://quarto.org) ≥ 1.4 — needs sudo on macOS; install
  interactively via `brew install --cask quarto` or the .pkg installer
  from quarto.org. The `tar_quarto(dashboard, ...)` target needs this on
  PATH.

R package dependencies are managed with `renv`; see `DESCRIPTION` for the
list. On first checkout:

```r
install.packages("renv")
renv::restore()
cmdstanr::install_cmdstan()   # first time only, ~10 min
```

## Running the pipeline

```r
targets::tar_make()           # full pipeline
targets::tar_visnetwork()     # inspect the DAG
targets::tar_load(fit_A)      # pull a specific target into the session
```

Fitting both Stan variants on the full sample takes ~1–2 hours on a modern
laptop (M-series Mac, 8 cores). The `targets` cache means data-prep changes
do not invalidate the fits; data updates do.

The rendered Quarto dashboard lands in `dashboard/_site/`.

## Project layout

See [`CLAUDE.md`](CLAUDE.md) for the canonical layout. In brief:

- `R/` — data fetch, prep, fitting, post-processing, viz helpers
- `stan/` — the three Stan models (`mct_aus_A.stan`, `mct_aus_B.stan`,
  `mct_aus_C.stan`)
- `dashboard/` — Quarto dashboard sources
- `data/` — raw ABS pulls (gitignored) and processed series
- `outputs/` — posterior draws (gitignored), figures, tables
- `tests/testthat/` — unit tests, including simulation-based Stan tests
- `_targets.R` — pipeline definition

## Key decisions

- **CPI panel**: the 11 standard ABS groups — Food & non-alcoholic beverages;
  Alcohol & tobacco; Clothing & footwear; Housing; Furnishings, household
  equipment & services; Health; Transport; Communication; Recreation & culture;
  Education; Insurance & financial services.
- **Reference sector for identification**: Housing (largest weight),
  $\lambda_{\text{Housing},t} = 1 \forall t$.
- **Sample start**: 1990 Q1 (latent monthly clock from 1990-01).
- **Demeaning window**: 2000–2019 per series.
- **Volatile sub-index handling**: no preprocessing — rely on stochastic
  volatility in the measurement noise to absorb supply-shock outliers
  (fuel, fruit & veg, holiday travel). Revisit only if posterior diagnostics
  demand it.
- **Expenditure weights source**: ABS *Consumer Price Index — 2025 Weighting
  Pattern* (released 26 Nov 2025, anchored to December 2024 expenditure
  shares), Table 1 (weighted average of eight capital cities). Cached at
  `data/raw/cpi_weights_2025.xlsx`. Used as fixed $w_i$ for the whole
  1990–present sample — see CLAUDE.md for the rationale.
  Source URL: https://www.abs.gov.au/system/files/83cf3a8e846b6b3f91ae7ea99ef968ab/Consumer%20Price%20Index%20-%202025%20Weighting%20Pattern.xlsx
- **Priors** (tightened from the brief's defaults — see commit history for
  the diagnosis): at T=435 the brief's `sigma_h ~ N(0, 0.2)` allowed
  pathological per-chain mode trapping (R̂ ≈ 18, ess_bulk ≈ 4). Current
  priors:
  - `mu_h* ~ N(0, 1)`, `sigma_h* ~ N(0, 0.1)+` (RW step ≤ 2 log-var units
    over the sample), `sigma_lambda ~ N(0, 0.01)+`
  - `rho ~ N(0.5, 0.3)` truncated to (−1, 1) — Variant A only
  - `s_init ~ N(0, 0.5)`, `lambda_init ~ N(1, 0.3)`, `c[1]` stationary for A
    / `N(0, 4)` for B
- **Initialisation**: all chains start at the same deterministic point
  (`R/fit_stan.R::default_init()` — all z innovations at 0, sigmas at small
  positive values, rho at 0.5). Eliminates the per-chain mode trapping
  Stan's default `U(-2, 2)` causes at long T. Stan's warmup-time random
  jitter still drives per-chain exploration.
- **Sampler config**: `chains = 4`, `iter_warmup = 1500`,
  `iter_sampling = 1500`, `adapt_delta = 0.97`, `max_treedepth = 12`.
  (`adapt_delta` raised from the brief's 0.95.)
- **Variant C scope**: implements the NY Fed MCT's structural choice of
  having BOTH a common-trend and a common-cycle factor (item 1 of the
  NY Fed comparison) plus a tighter log-vol drift prior matching their
  scale (item 5). Two NY Fed features are **deferred**:
  - **MA(3) sector idiosyncratic errors** (item 2). Stan with MA(q) on
    sector noise via the forward-recursion idiom blew up identification
    on T=120 simulated data — NUTS could not reliably keep the MA roots
    inside the unit circle even with tightly shrunk priors. NY Fed's
    Gibbs sampler enforces invertibility post-draw via root reflection,
    which has no straightforward HMC analogue. Would require either a
    partial-autocorrelation reparameterisation (q=3 case is tractable
    but non-trivial) or moving off Stan to a custom sampler.
  - **Outlier scale mixture / Student-t robust likelihood** (item 3).
    Combined with the two-common-factor structure, the extra distributional
    flexibility led to chronic warmup multi-modality (50% of chains stuck
    at ebfmi < 0.1). NY Fed sidesteps this via slice-sampled Bernoulli
    outlier indicators in their Gibbs sampler; the HMC analogue is
    marginalisation, but the resulting posterior geometry was infeasible
    on the test cases used.
  Both deferred items are worth pursuing in a v2 — either via PAC
  reparameterisation + careful sampler tuning, or by porting the NY Fed's
  Gibbs sampler to R (see `MCT-NYFED/functions/`).
- **Variant C multi-modality**: even with all simplifications, ~1-in-4
  chains gets stuck during warmup in an alternate factor-allocation mode
  on small samples. The healthy chains converge to the right answer
  (95%+ posterior coverage in simulation); the stuck chain inflates R-hat.
  `fit_mct()` emits a warning identifying chains with `ebfmi < 0.3`;
  discard those before summarising. On the full T=435 real sample, more
  data should make this rarer — verify post-fit.
