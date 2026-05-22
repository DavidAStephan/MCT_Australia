# Australian Multivariate Core Trend Inflation

A Bayesian mixed-frequency dynamic factor model for Australian CPI, in the
spirit of the NY Fed's [Multivariate Core Trend (MCT)](https://www.newyorkfed.org/research/policy/mct)
inflation model (Stock & Watson 2016). Two variants of the common factor are
estimated and compared.

See [`CLAUDE.md`](CLAUDE.md) for the full project brief and model specification.

## Status

🚧 In progress — Step 1 (project skeleton).

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
- `stan/` — the two Stan models (`mct_aus_A.stan`, `mct_aus_B.stan`)
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
