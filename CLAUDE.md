# Project: Australian Multivariate Core Trend (MCT) Inflation

## Context

I'm building an Australian analogue of the NY Fed's Multivariate Core Trend (MCT)
inflation model (https://www.newyorkfed.org/research/policy/mct). The NY Fed model
is a Stock-Watson (2016) dynamic factor model with time-varying parameters and
stochastic volatility, fit to 17 monthly PCE sub-indices, decomposing each
sector's inflation into a common trend, sector-specific trend, common transitory
shock, and sector-specific shock. The output is (a) a measure of trend inflation,
and (b) a decomposition of where persistence is coming from (common vs sector-
specific, which sectors are dominant).

The Australian challenge is that the ABS only began publishing a complete monthly
CPI in October 2025 (back-series to April 2024). Quarterly CPI by group goes back
to the 1980s. The model therefore needs to handle **mixed-frequency** data
honestly: quarterly observations pre-April-2024, monthly observations after,
estimated jointly in a single state-space model.

Reference reading:
- Stock & Watson (2016), "Core Inflation and Trend Inflation", REStat 98(4)
- NY Fed MATLAB code: https://github.com/MCT-Inflation-NYFed/MCT-PCE
- ABS monthly CPI documentation: https://www.abs.gov.au/statistics/economy/price-indexes-and-inflation/consumer-price-index-australia/latest-release
- RBA Nov 2025 SMP Technical Note: https://www.rba.gov.au/publications/smp/2025/nov/technical-note-the-transition-to-a-complete-monthly-cpi.html

## Goal

A reproducible R pipeline that:

1. Fetches and cleans Australian CPI data (quarterly historical + monthly recent)
2. Estimates a Bayesian mixed-frequency MCT in Stan via `cmdstanr`
3. Estimates **two model variants** of the common factor and compares them
4. Produces a Quarto dashboard summarising trend inflation, common-vs-sector
   decomposition, and sector contributions
5. Re-runs automatically on the day of each new ABS Monthly CPI release

Everything orchestrated with `{targets}` so partial re-runs are cheap.

## Model specification

Let $\pi_{i,t}$ be the annualised, demeaned month-on-month inflation of CPI group
$i \in \{1,\dots,N\}$ at month $t$, on a latent monthly clock. There are $N=11$
ABS CPI groups in the baseline.

**Measurement equation (monthly latent):**

$$\pi_{i,t} = \lambda_{i,t} \, c_t + s_{i,t} + \exp(h^\varepsilon_{i,t}/2)\,\varepsilon_{i,t}$$

**Mixed-frequency observation:** when only quarterly data exists at month $t$
(where $t$ is the last month of a quarter), the observation equals the average
of the three latent monthly inflations $(\pi_{i,t-2} + \pi_{i,t-1} + \pi_{i,t})/3$.

**State equations (random walks):**

- $\lambda_{i,t} = \lambda_{i,t-1} + \sigma_\lambda \eta^\lambda_{i,t}$  (time-varying loadings)
- $s_{i,t} = s_{i,t-1} + \exp(h^s_{i,t}/2)\,\eta^s_{i,t}$  (sector trends with SV)
- $h^c_t = h^c_{t-1} + \sigma_{h_c}\eta^{h_c}_t$  (SV of common factor)
- $h^s_{i,t}, h^\varepsilon_{i,t}$  analogous

**Two common-factor variants (estimate both, compare):**

- **Variant A — AR(1) transitory:** $c_t = \rho c_{t-1} + \exp(h^c_t/2)\varepsilon^c_t$,
  $|\rho|<1$. "Common transitory shocks comove across sectors."
- **Variant B — Random-walk common trend (Stock-Watson):**
  $c_t = c_{t-1} + \exp(h^c_t/2)\eta^c_t$. "Sectors share a common drifting trend."

These imply different decompositions of trend vs cycle. Variant A's "trend" is
$\tau_t = \sum_i w_i s_{i,t}$; Variant B's "trend" is
$\tau_t = (\sum_i w_i \lambda_{i,t})\cdot c_t + \sum_i w_i s_{i,t}$.

**Identification:**

- Normalise $\lambda_{1,t} = 1 \forall t$ (Housing as reference sector, largest weight)
- All RWs initialised with weak priors near zero or one as appropriate
- Variant A: stationary initial distribution for $c_1$
- Variant B: $c_1 \sim N(0, 2^2)$ as a diffuse start

## Concrete deliverables

```
mct-aus/
├── CLAUDE.md                    # this file (project brief)
├── README.md                    # how to run, what the outputs mean
├── _targets.R                   # targets pipeline
├── DESCRIPTION                  # treat as a package: deps via renv
├── renv.lock
├── R/
│   ├── data_fetch.R             # readabs + direct downloads
│   ├── data_prep.R              # build (T x N) y matrix and obs_type
│   ├── fit_stan.R               # cmdstanr wrappers, both variants
│   ├── postprocess.R            # extract trend, common_share, sector contribs
│   ├── compare_variants.R       # WAIC / LOO + side-by-side trend chart
│   └── viz.R                    # ggplot helpers for dashboard
├── stan/
│   ├── mct_aus_A.stan           # AR(1) common factor
│   └── mct_aus_B.stan           # RW common trend
├── dashboard/
│   ├── index.qmd                # main dashboard page
│   ├── _quarto.yml
│   └── styles.scss
├── data/
│   ├── raw/                     # cached ABS downloads with vintage timestamp
│   └── processed/
├── outputs/
│   ├── draws/                   # posterior draws (parquet)
│   ├── figures/                 # png exports
│   └── tables/
├── .github/
│   └── workflows/
│       └── update.yml           # monthly re-run on ABS release day
└── tests/
    └── testthat/                # tests on data prep + Stan model on simulated data
```

## Step-by-step plan

Work through these in order. Don't move on until each step works.

### Step 1: project skeleton
- Initialise as an R project with `renv`. Pin: `cmdstanr`, `posterior`,
  `bayesplot`, `targets`, `tarchetypes`, `dplyr`, `tidyr`, `lubridate`,
  `readabs`, `ggplot2`, `arrow`, `quarto`, `loo`, `testthat`.
- Set up `.gitignore` for `data/raw/`, `outputs/draws/`, `.cmdstan/`.
- Write the `README.md` with run instructions.

### Step 2: data ingest
- `R/data_fetch.R`:
  - `fetch_cpi_quarterly()`: pull quarterly CPI by group from ABS catalogue
    6401.0 (table 7 for index numbers by group). Cache by ABS release date.
  - `fetch_cpi_monthly()`: pull the complete Monthly CPI (the new publication
    that started Oct 2025; back-series to April 2024). NOT the old Monthly CPI
    Indicator, which has been discontinued.
  - Both functions return tidy long-format data frames with columns
    `date, group, index, source`.
- `R/data_prep.R`:
  - `build_inflation_series()`: convert index to annualised period-on-period
    log-differenced inflation. Annualisation factor: 4 for quarterly, 12 for
    monthly.
  - `demean_series()`: subtract the 2000–2019 mean per series.
  - `build_stan_data()`: produce the $(T \times N)$ matrices `y` and `obs_type`
    on a monthly latent clock from 1990-01 (or earliest available) to the
    current month. Quarterly observations are placed at the last month of the
    quarter with `obs_type = 2`; monthly observations from April 2024 onward
    have `obs_type = 1`; everything else is 0.
  - Tests: confirm that summing 3 consecutive months in the simulated case
    recovers the corresponding quarterly mean (sanity check on the encoding).

### Step 3: Stan models
- `stan/mct_aus_A.stan` (AR(1) common): start from the spec in this brief.
  Use non-centred parameterisations throughout (random walks as
  `x[t] = x[t-1] + sigma * z[t]` with `z ~ std_normal()`).
- `stan/mct_aus_B.stan`: same as A but with the common factor as a random
  walk (drop `rho_c`, replace `c[t] = rho * c[t-1] + ...` with
  `c[t] = c[t-1] + ...`).
- Both models:
  - Stochastic volatility on common factor, sector trends, and measurement
    noise (log-variance random walks with `mu`, `sigma` hyperparams)
  - Time-varying factor loadings as random walks
  - `lambda[t, 1] = 1` hard-coded for identification (the reference sector)
  - `generated quantities` block: `trend[T]`, `common_share[T]`,
    `sector_trend[T, N]`, `common_transitory[T]`

### Step 4: simulation test
- `tests/testthat/test-stan-simulated.R`:
  - Simulate data from each Stan model with known parameters
  - Fit the model back and check posterior medians fall within 95% intervals
    of the truth for headline `trend` and `common_share` at most $t$
  - Run with small T (e.g. 120 months, 5 sectors) to keep it fast
- This is the single most important step. Don't skip it.

### Step 5: fitting on real data
- `R/fit_stan.R`:
  - `fit_mct(stan_data, variant = c("A", "B"), ...)`:
    `chains = 4`, `parallel_chains = 4`, `iter_warmup = 1500`,
    `iter_sampling = 1500`, `adapt_delta = 0.95`, `max_treedepth = 12`
  - Save the full `CmdStanMCMC` object via `$save_object()` and also export
    a slimmer parquet of the key generated quantities for the dashboard
  - On a development machine, fit 1990–present for ~12 sectors. Budget
    30–90 min per variant per fit; document the runtime in the README.

### Step 6: post-processing
- `R/postprocess.R`:
  - `summarise_path(draws, dates, probs = c(0.16, 0.5, 0.84))`: return tidy
    date/median/lower/upper tibble. Use a wider band (e.g. 0.05/0.95) as
    optional.
  - `sector_contributions(fit, w)`: return $w_i \cdot s_{i,t}$ per date and
    sector for stacked-area charts.
  - `common_share_path(fit)`: extract the headline "how much is common" series.
  - `posterior_revisions(current_fit, prior_fit)`: compare trend estimates
    from this month's fit vs last month's fit at the same dates — feeds the
    "revisions" panel of the dashboard.

### Step 7: variant comparison
- `R/compare_variants.R`:
  - LOO comparison via `loo::loo_compare()` — note these models have many
    states so PSIS-LOO may flag high Pareto-k; report what they say and the
    caveats. Use point-wise log-likelihood from the `generated quantities`
    block (add a `log_lik[T, N]` array to the Stan files).
  - Side-by-side trend chart with credible bands from both variants
  - Side-by-side common-share chart
  - One-paragraph write-up dropped into the dashboard explaining which
    variant fits better and where they disagree

### Step 8: targets pipeline
- `_targets.R`:
  - Targets: `cpi_q_raw`, `cpi_m_raw`, `cpi_q_clean`, `cpi_m_clean`,
    `stan_data`, `fit_A`, `fit_B`, `trend_A`, `trend_B`, `common_share_A`,
    `common_share_B`, `sector_contribs_A`, `sector_contribs_B`,
    `variant_comparison`, `dashboard`
  - Use `tarchetypes::tar_quarto()` for the dashboard target
  - Stan fits cached by hash of input data — they don't re-run if the data
    hasn't changed
- Verify `tar_visnetwork()` shows the right DAG before declaring done.

### Step 9: Quarto dashboard
- `dashboard/index.qmd`:
  - **Top row**: headline trend (median + 68% band) for variant B, with
    current month annotation
  - **Second row**: common-share over time (when does common drive
    inflation, vs sector-specific noise)
  - **Third row**: stacked area of current sector contributions to trend
    inflation (last 24 months)
  - **Fourth row**: revisions panel — how this month's estimate compares
    to last month's estimate at the same dates
  - **Fifth row**: variant comparison (both trends overlaid, with LOO note)
  - **Methodology section** at the bottom: collapsible, with the model
    equations rendered in LaTeX
- Style for readability — modest fonts, clean grid, no chart-junk.
  Inspiration: Liberty Street Economics blog posts.

### Step 10: auto-update
- `.github/workflows/update.yml`:
  - Trigger: scheduled (last Wednesday of each month at 12:30 AEST, ~1 hour
    after the ABS release at 11:30) AND manually dispatchable
  - Steps: checkout, install R + cmdstan, restore renv, `targets::tar_make()`,
    publish the rendered Quarto site to `gh-pages` branch
  - Note Stan compile time + sampling time — likely needs a 2-hour timeout
    and a larger runner. If GitHub-hosted runners are too slow, sketch the
    equivalent setup for a self-hosted runner

## Working style

- Default to small, focused commits. One commit per step above.
- After each Stan model edit, run the simulation test before fitting on real
  data. The state-space is non-trivial; bugs are easier to catch on simulated
  data with known truth.
- Use `posterior::summarise_draws()` to inspect divergences, `rhat`, and
  `ess_bulk` after each fit. Don't move on with `rhat > 1.01` or many
  divergences without understanding why.
- For long Stan compile times, leave the model files alone unless changing
  the spec — recompiling is the slow part of iteration.
- Document choices in `README.md` as you make them. Especially: the
  expenditure weights you settled on, the demeaning window, identification
  choices, and any pre-processing of fuel/food.

## Things to ask me before deciding

- The exact ABS expenditure weights to use (I'll provide the source)
- Whether to drop or downweight any volatile sub-indices (fuel, fruit & veg,
  holiday travel) — the model should handle volatility via SV, but it might
  still be worth a pre-processing step
- Whether to start the data in 1990, or pre-1990 if available
- Any prior choice that ends up driving results — flag it
