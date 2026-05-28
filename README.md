# Australian Multivariate Core Trend Inflation

A Bayesian dynamic factor model for Australian CPI sub-groups, in the
spirit of the NY Fed's [Multivariate Core Trend (MCT)](https://www.newyorkfed.org/research/policy/mct)
inflation model (Stock & Watson 2016). The production model is **Variant B
(random-walk common trend)** fit by an in-house Gibbs sampler on the 33
ABS CPI sub-groups at quarterly frequency.

See [`CLAUDE.md`](CLAUDE.md) for the full project brief.

## Status

Variant B (Stock-Watson RW common trend) fit by Gibbs sampling on the
T=133 quarterly sample (1993-Q1 → 2026-Q1, RBA inflation-targeting era),
N=33 ABS CPI sub-groups. Wall time ~4-5 min on M1.

### Headline finding (2026-05-28)

At the **sub-group quarterly panel**, AU CPI inflation is dominated by a
common drift component. Latest-quarter (March 2026):

- **Trend inflation**: 0.33 (68% CI [-0.51, 1.21]), demeaned vs 1993-2019.
- **Common share of variance**: 81% (68% CI [69%, 89%]) — most of the
  recent inflation move is broad-based across sub-groups, not driven by
  any single category.

The common-share trajectory has declined from ~97% in the mid-1990s to
~81% in the 2020s — sectors have become somewhat more idiosyncratic over
the inflation-targeting era, but the common drift is still the dominant
driver. This is **consistent with the NY Fed finding for US PCE** and
matches the original Stock-Watson story.

#### Previous vintages (archived)

Earlier vintages of this project fit AR(1)-cycle (Variant A), RW common
trend (Variant B) and combined (Variant C) variants on a smaller
**11-group monthly mixed-frequency panel** using Stan. On that coarser
panel Variant A appeared to win LOO (A > B > C) and common-share was
weakly identified at ~14%. That finding was a **granularity artefact**:
group-level aggregation washes out the common signal that sub-group
analysis reveals. The 11-group fits are preserved in
`outputs/draws/legacy_11group/` for reference. The original 1990-onward
fits (before the inflation-targeting-era cutover) are in
`outputs/draws/legacy_1990/`.

## Requirements

- R ≥ 4.3
- [Quarto](https://quarto.org) ≥ 1.4 — needs sudo on macOS; install
  interactively via `brew install --cask quarto` or the .pkg installer
  from quarto.org. The `tar_quarto(dashboard, ...)` target needs this on
  PATH.

R package dependencies are managed with `renv`; see `DESCRIPTION` for the
list. **No CmdStan needed** — the production sampler is the in-house
Gibbs implementation in `R/gibbs/`, which depends only on base R + MASS.
On first checkout:

```r
install.packages("renv")
renv::restore()
```

## Running the pipeline

```r
targets::tar_make()           # full pipeline (data fetch + fit + dashboard)
targets::tar_visnetwork()     # inspect the DAG
targets::tar_load(fit_B)      # pull a specific target into the session
```

The full pipeline runs in ~10 min on a modern laptop (M-series Mac) on a
cold cache — the Gibbs fit itself is ~4-5 min; the rest is data fetch +
dashboard render. The `targets` cache means data-prep changes do not
invalidate the fit; data updates (new ABS release) do.

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

- **CPI panel**: 33 ABS CPI sub-groups (one level below the 11 standard
  groups), e.g. "Rents", "New dwelling purchase by owner-occupiers",
  "Bread and cereal products", etc. Closer in granularity to NY Fed's 17
  PCE sub-indices. See `abs_cpi_subgroups()` in `R/data_fetch.R` for the
  full list. The previous 11-group production model is archived in
  `outputs/draws/legacy_11group/`.
- **Frequency**: quarterly latent clock matching the data. ABS
  sub-group breakdowns are quarterly only (no monthly sub-group
  publication). A monthly latent clock would introduce ~3:1
  states-to-obs underidentification, which is what the failed first
  attempt revealed; quarterly avoids the problem entirely.
- **Reference sub-group for identification**: Rents (w=6.60 in the 2025
  weighting pattern), $\lambda_{\text{Rents},t} = 1 \forall t$. Chosen
  over the larger "New dwelling purchase by owner-occupiers" (w=7.54)
  because the latter only starts in 1998-06, which would break the
  $\lambda_{ref}$ identification on the 1993-onward sample.
- **Sample start**: 1993-Q1 (RBA inflation-targeting era — Fraser's
  March 1993 "2-3% on average over the cycle" speech). Pre-1993 ABS
  quarterly data exists back to 1948 but the 1990-1992 disinflation is
  a one-off regime change that biases the SV processes. The
  `legacy_1990/` folder under `outputs/draws/` retains the prior-vintage
  fits.
- **Demeaning window**: 1993–2019 per series. Full pre-pandemic in-sample
  period. Sensitivity to the exact window endpoints is small (~0.2pp at
  the weighted-aggregate mean).
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
- **Sampler**: in-house Gibbs implementation in `R/gibbs/`, a custom port
  of the NY Fed MCT MATLAB sampler. Sweeps the joint posterior in
  conventional Stock-Watson order: FFBS for the latent (common + sector)
  trends, Kim-Shephard for the stochastic-volatility processes,
  Beta-conjugate for the discrete scale-mixture parameters (when
  outliers are enabled), Gaussian-conjugate for the constant loadings.
  Production config: `n_burn = 1000`, `n_draw = 2000`, single chain
  (~4-5 min wall time at N=33, T=133).
- **Priors**:
  - `sigma_c[1], sigma_s[i,1], sigma_eps[i,1]` initialised at 0.5
    (typical SV scale for annualised AU inflation)
  - `gamma_*` (RW step on log-variance): Inverse-Gamma prior with
    `nu = 60`, `s2 = 1 / (60 * 4)` calibrated for quarterly cadence
    (annual implied volatility-of-log-vol ≈ 1/60)
  - `lambda_i ~ N(1, 0.3)` (sub-group loading)
  - `c[1] ~ N(0, 2^2)` — diffuse Variant B start (Stock-Watson canonical)
  - `s_{i,1} ~ N(0, var_s_init = 0.25)`
- **Initialisation**: deterministic single-chain start (`fit_mct_gibbs()`
  default `init`) — all sector trends and common factor at 0, sigmas at
  0.5, gammas at 0.05, lambda at 1.
- **Sampler dispatch**: Variant A (AR(1) common, with `rho` update) vs
  Variant B (RW common, no `rho` update — production model) selected by
  `config$variant`. See `R/gibbs/gibbs_sweep.R` and
  `R/gibbs/build_ssm.R`.
