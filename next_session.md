# Next session plan — MCT_Australia

Hand-off note for the next working session. Pair with [`CLAUDE.md`](CLAUDE.md)
(project brief), [`README.md`](README.md) (current decisions + headline
findings), and the auto-memory at
`~/.claude/projects/-Users-davidstephan-Documents-MCT-Australia-MCT-Australia/`.

## Session update (2026-05-24)

- Dashboard now renders end-to-end via `Rscript -e 'targets::tar_make()'`.
  Output at `dashboard/_site/index.html`.
- Targets pipeline refactored so the three fits don't duplicate into
  `_targets/objects/`. `fit_{A,B,C}` are `format = "file"` targets returning
  the absolute path to `outputs/draws/fit_X.rds`; downstream targets call
  `readRDS()` inline. End result: `_targets/` stays under 1 MB; the only
  copies of the fits on disk are the canonical ones in `outputs/draws/`.
  To refit a new vintage, see the comment block above the fit_X targets
  in `_targets.R`.
- Deleted as unused: `outputs/draws/fit_{A,B,C}_slim.parquet` (1.5 GB total).
  The slim-parquet exporter (`export_fit_parquet`, `read_slim_parquet`) is
  still present in `R/` but has no live callers — fine to delete next
  cleanup pass.
- `tar_quarto(dashboard, ..., working_directory = "dashboard")` is now
  required because `dashboard/_quarto.yml` makes quarto treat `dashboard/`
  as the project root; without it tar_quarto and quarto disagree on output
  paths. Project-root detection added at the top of `dashboard/index.qmd`
  so both `tar_make` and standalone `quarto render dashboard/index.qmd`
  work from any pwd.
- Disk total dropped from ~14.5 GB → ~6.5 GB after cleanup. Future
  `tar_make` runs won't re-duplicate.
- Remote is already configured (`DavidAStephan/MCT_Australia`); `main`
  pushed. Action #3 in the original handoff is done.

## State at session end (2026-05-23)

- All 10 brief steps landed; Variant C added on top (combined common-trend
  + common-cycle, based on NY Fed MCT structure).
- Three Stan fits cached at `outputs/draws/fit_{A,B,C}.rds` (gitignored).
  C: 108 min on 6 chains × 1500+1500, clean (R̂ ≤ 1.00, ESS_bulk median ~13K,
  2/9000 divergent, ebfmi 0.90–1.04 all chains, watchdog didn't trip).
- 3-way LOO comparison computed: **A ≻ B ≻ C**.
  - A: 0.0 elpd_diff (best)
  - B: -13.0 (SE 4.4) — A > B by ~3σ
  - C: -71.6 (SE 11.2) — A > C by ~6σ
- Item 4 (revisions decomposition by sector) shipped — works on cached fits,
  wired into dashboard. Active once a prior-vintage fit is rotated in.
- Items 2 (MA(3) sector errors) and 3 (outlier scale mixture / Student-t)
  from the NY Fed comparison are **deferred** — pilot tests showed HMC
  cannot mix them. See README and `stan/mct_aus_C.stan` comments for the
  rationale. Path forward: port NY Fed Gibbs sampler to R rather than fight
  HMC further.
- No remote configured yet — `git push -u origin main` is the first step
  once a GitHub repo exists.

## Headline empirical finding (carry into v2 work)

Australian sectoral CPI co-movement is **overwhelmingly transitory**, not a
drifting common trend. This is opposite to NY Fed's finding for US PCE.
Variant C identifies both factors (all 10 free `lambda_eps` loadings have
5-95% intervals excluding zero; cycle log-vol RW step `sigma_hec=0.17` is
~3× the trend's `sigma_hc=0.065`) but over-parameterises and loses LOO
decisively to the simpler AR(1) cycle (A).

**Implication:** don't pursue further combined-factor variations as the
natural next step. The LOO signal is strong.

Recent month medians (March 2026, demeaned annualised %):
- A trend: 0.65 (68% band [-0.08, +1.41])
- B trend: 0.89 (wider, [-0.07, +1.85])
- C trend: 1.02 (widest, [+0.08, +1.96])

Add ~2.5pp back for headline % p.a. (the 2000-2019 mean).

## Immediate actions (priority order)

### 1. Install Quarto (still blocks dashboard rendering)

Quarto CLI still not on PATH. Run in a terminal with sudo:

```sh
brew install --cask quarto    # interactive password prompt
# OR download the .pkg from https://quarto.org/docs/get-started/
quarto --version              # verify ≥ 1.4
```

Then render the dashboard:

```sh
Rscript -e 'targets::tar_make(dashboard)'
```

The dashboard is wired for three variants now: top panels stay on A
(headline), variant-comparison panel overlays A/B/C with the LOO table
and the auto-generated comparison paragraph.

### 2. Reconcile `targets` cache with the manual fit_C

`fit_C` was fitted standalone via `fit_mct(..., variant = "C")`, saved to
`outputs/draws/fit_C.rds` directly. The `targets` cache (under `_targets/`)
doesn't know about it. Running `targets::tar_make()` will see that the
`fit_C` target is uncomputed and try to refit — **a 1.8-hour repeat**.

Options:

- **Best**: import the existing fits into the targets cache via
  `tar_target_raw()` or `tar_load_raw()` so targets sees them as fresh.
  This is a one-time fiddle but avoids the duplicate refit.
- **Pragmatic**: run `targets::tar_make()` once and let it refit
  everything. fit_A is ~93 min, fit_B is ~227 min (or attempt td=14 again,
  abandoned earlier — see memory), fit_C is ~108 min. Total ~7 hours.
- **Cheapest**: skip the `fit_*` targets and just `tar_make()` the
  downstream targets (`trend_*`, `loo_*`, `variant_comparison`, etc.)
  using `tar_target_raw` with manual `command = ...`. Requires editing
  `_targets.R`.

The easiest path may be the pragmatic one if there's no time pressure.

### 3. Push to GitHub when ready

```sh
git remote add origin git@github.com:<owner>/mct-aus.git
git push -u origin main
```

The Actions workflow in `.github/workflows/update.yml` activates on first
push. First cron fires the next last-Wednesday-of-the-month.

### 4. Optional: rotate fit_A_prior.rds for the revisions panel

The revisions and revision-attribution panels in the dashboard go live
once `outputs/draws/fit_A_prior.rds` + `outputs/draws/stan_data_prior.rds`
exist. On the next monthly update:

```r
# After fitting against a new vintage, BEFORE overwriting fit_A.rds:
file.copy("outputs/draws/fit_A.rds",       "outputs/draws/fit_A_prior.rds",
          overwrite = TRUE)
file.copy("outputs/draws/stan_data.rds",   "outputs/draws/stan_data_prior.rds",
          overwrite = TRUE)
# Then refit on new data → outputs/draws/fit_A.rds gets the new vintage.
```

This rotation should be automated in the GH Actions workflow eventually.

## Substantive findings to remember (already in memory + README)

- **3-way LOO ranking**: A ≻ B ≻ C. Australian common factor reads
  transitory, not drift, AND combining trend + cycle doesn't help.
- **Variant C diagnostics on real data**: 6/6 healthy chains, ebfmi
  0.90–1.04, 2/9000 divergent, ESS_bulk median ~13K, R̂ ≤ 1.00. The
  multi-modality from simulation tests on T=120 didn't materialise on
  T=435 — more data, more identification.
- **Latest-month medians** (March 2026, demeaned annualised %, headline=A):
  trend median 0.65, 68% band [-0.08, +1.41]. Common-share (variance ratio)
  median ~7%. Add 2.5pp for raw %.
- **Sampler config that works at T=435**:
  - A: 4 chains × 1500+1500, `sigma_h ~ N(0, 0.1)`, `adapt_delta = 0.97`,
    deterministic init via `default_init()`.
  - B: same as A; warmup spent ~25% of post-warmup transitions at
    treedepth=12. td=14 retry was attempted, ran ~70% before being killed
    for taking too long (~6h projected). Don't retry without a plan.
  - C: 6 chains × 1500+1500, `sigma_h ~ N(0, 0.04)`, deterministic init,
    `adapt_delta = 0.97`. Took 108 min total.
- **Watchdog**: `fit_mct()` now emits a warning if any chain has
  `ebfmi < 0.3` after sampling. Discard those chains before summarising.

## Open decisions deferred for later (v2)

- **NY Fed items 2 & 3 (MA(3) + Student-t)**: deferred. If revisited, port
  the NY Fed Gibbs sampler from MATLAB (`MCT-NYFED/functions/`) to R. The
  Kim-Shephard SV mixture and slice-sampled outlier indicators are
  fundamentally incompatible with HMC at the scale we need. Estimated
  effort: days of R porting work. Run-time would be ~10–25 min single
  chain on T=435 (vs Stan's 100+ min for the simpler models).
- **Vintage-accurate weights**: still using 2025 weights for the whole
  sample. Pre-2017 the ABS used 5-yearly weight links; vintage-appropriate
  weights are a nice-to-have but trend-decomposition story is robust to
  ±2pp shifts.
- **Volatile sub-index treatment**: still no preprocessing. Variant A
  handles the volatile sectors via SV alone. A sensitivity run with
  winsorised fuel / fruit-&-veg / travel is worth doing once for the
  methodology section.
- **Sub-group panel (~30–50 series)**: closer to NY Fed's 17 PCE
  sub-indices in spirit. Better factor identification but heavier compute
  and patchier history. Group-level (N=11) is the v1 baseline.
- **Pareto-k > 0.7 observations**: many state-space LOO comparisons hit
  this (132 / 113 / 76 for A / B / C of 1698). K-fold CV via `loo::kfold()`
  would be more defensible — slow but doesn't rely on the importance-
  sampling assumption that fails for these states.

## Quick orientation for the next session

```sh
# Resume work — branches/commits
git log --oneline | head
git status

# Re-load the cached objects without re-fitting
Rscript -e '
  for (f in list.files("R", "\\.R$", full.names=TRUE)) source(f)
  stan_data <- readRDS("outputs/draws/stan_data.rds")
  fit_A <- readRDS("outputs/draws/fit_A.rds")
  fit_B <- readRDS("outputs/draws/fit_B.rds")
  fit_C <- readRDS("outputs/draws/fit_C.rds")
  # ready to call trend_path(), common_share_path(), compare_*_plot(), etc.
'

# Re-run the 3-way LOO comparison
Rscript -e '
  for (f in list.files("R", "\\.R$", full.names=TRUE)) source(f)
  loo_A <- loo_for_fit(readRDS("outputs/draws/fit_A.rds"))
  loo_B <- loo_for_fit(readRDS("outputs/draws/fit_B.rds"))
  loo_C <- loo_for_fit(readRDS("outputs/draws/fit_C.rds"))
  print(compare_fits_loo(list(A = loo_A, B = loo_B, C = loo_C)))
'

# Run the test suite (fast bits)
Rscript -e 'testthat::test_dir("tests/testthat")'

# Run the simulation test (slow — ~15 min total for all three variants)
RUN_STAN_TESTS=1 Rscript -e 'testthat::test_dir("tests/testthat")'

# Full pipeline once Quarto is installed
Rscript -e 'targets::tar_make()'
```
