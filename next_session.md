# Next session plan — MCT_Australia

Hand-off note for the next working session. Pair with [`CLAUDE.md`](CLAUDE.md)
(project brief), [`README.md`](README.md) (current decisions + headline
findings), and the auto-memory at
`~/.claude/projects/-Users-davidstephan-Documents-MCT-Australia-MCT-Australia/`.

## Session update (2026-05-27 evening) — Steps 12 + 13: MA(q) + outliers ported

Major model-quality additions to the Gibbs port. **120/120 tests pass.**

### Step 12 — MA(q) measurement errors (monthly-only)

State-augmented SSM with (q+1) eps lags per sector. NY Fed convention
(eps state's variance = sigma_eps^2). 3 new files:
[`R/gibbs/update_theta.R`](R/gibbs/update_theta.R),
[`R/gibbs/build_ssm_ma.R`](R/gibbs/build_ssm_ma.R),
[`R/gibbs/gibbs_sweep_ma.R`](R/gibbs/gibbs_sweep_ma.R).

Use via `fit_mct_gibbs(y, q_MA = 3, ...)`. **LIMITATION**:
monthly-only — combined mixed-freq + MA path NOT implemented (state
dim would be ~92 at N=11, q=3). For AU production this means MA
can only fit the post-Apr-2024 monthly slice (T=23). Sim-validated
capability; not yet production-ready for AU.

### Step 13 — Outlier scale mixture (mixed-freq compatible)

Discrete per-obs scales (40-point grid 1→10) inflate the obs noise.
Integrates cleanly with quarterly observations because scales are
per-obs scalars, not state-augmented. 2 new files:
[`R/gibbs/update_scl.R`](R/gibbs/update_scl.R),
[`R/gibbs/update_ps.R`](R/gibbs/update_ps.R).

Use via `fit_mct_gibbs(y, obs_type = ot, config = list(use_outliers = TRUE), ...)`.
Tests verify it flags 3+/5 injected outliers with <15% false-positive
rate. **READY FOR REAL AU DATA** — directly addresses CLAUDE.md's
flagged volatile sub-indices (fuel, fruit & veg, holiday travel).

### Headline result

The Gibbs port is now methodologically AT PARITY with NY Fed's MCT
model on the things that matter for AU (mixed-freq + outliers), with
MA(q) available as a monthly-only capability for future use.

**Next session targets** (in order of value):

1. **Real-data outlier benchmark** — fit `fit_mct_gibbs(use_outliers = TRUE)`
   on the cached `outputs/draws/stan_data.rds` (real AU CPI, T=435,
   N=11). Expected ~25-30 min single chain. Check posterior median
   `s_outlier` for known-volatile sectors (transport, recreation,
   alcohol) at fuel-spike dates (2022-2023) — should flag them.

2. **LOO comparison: baseline vs +outliers** — compute pointwise
   log-likelihood for both fits and compare via loo::loo_compare.
   Tells us empirically whether outliers improve AU model fit.

3. **(Stretch) Extend MA to mixed-freq** — would unlock MA on full
   T=435 AU data instead of just T=23 monthly tail. State dim ~92;
   needs careful implementation but no new conceptual work.

4. **(Stretch) Refit Stan Variant Ac for the apples-to-apples
   comparison** that was deferred at end of Step 11.

## Session update (2026-05-27 morning) — Step 11: Gibbs validated; "rho bias" was a faulty premise

**Correction to previous sessions:** I had claimed Stan's rho on real
ABS data was ~0.7, then spent time diagnosing/fixing a "rho bias" in
the Gibbs sampler. **Stan's actual rho posterior on real data is
0.371.** Gibbs gives 0.360 on the same data. They're statistically
identical.

Long-chain comparison (`scripts/step10b_gibbs_vs_stan_real.R`,
2000 burn + 5000 draws, 22 min):

| Quantity | Stan Variant A | Gibbs | Verdict |
|---|---|---|---|
| **rho median** | **0.371** [0.18, 0.54] | **0.360** [0.17, 0.53] | **identical** ✓ |
| Trend latest-t | 0.647 | 1.059 | differs — **NOT a sampler bug** |
| Trend cor across t | — | 0.678 | shape match |
| Wall time | 93 min | 22 min (MH) / 15 min (cond) | 4-6× faster |

**Trend gap is a model-spec difference**, not a sampler issue:
- Stan Variant A: TIME-VARYING lambda (RW)
- Gibbs: CONSTANT lambda (matches Stan Variant Ac)

With constant λ, less variation absorbed via λ_t * c_t ⇒ more into
sector trends ⇒ trend (sum w_i * s_i) is higher.

**The marginal MH rho update was implemented** anyway
([`R/gibbs/update_rho_mh.R`](R/gibbs/update_rho_mh.R)), default
`use_marginal_mh_rho = TRUE` for safety. Adds ~50% wall-time overhead
(one KF pass per iter). Set FALSE to skip — gives the same posterior.

**Total Gibbs port state: 91 tests passing across `R/gibbs/`.**

**Cleared next-session targets** (in order of value):

1. **Refit Stan Variant Ac to get the apples-to-apples comparison.**
   The Ac fit was deleted in the May speedup-exploration. With Ac
   re-cached, the Stan-vs-Gibbs trend should match exactly (both use
   constant lambda). ~5-10 min on M1 via the cached pipeline.

2. **OR add TVP lambda to the Gibbs sweep** using
   [`R/gibbs/update_tvcoef.R`](R/gibbs/update_tvcoef.R) which we
   already ported but never wired in. Then Stan Variant A vs Gibbs
   should match directly. Effort: ~half day.

3. **Wire Gibbs into `_targets.R`** as an alternative pipeline target
   (e.g. "fit_Ag"). Drop-in via the existing mct_gibbs_fit S3 wrapper.
   Worth doing once #1 or #2 confirms full validation.

4. **Stretch: Kim-Shephard SV update timing.** The Gibbs is already
   18× faster than Stan even without optimisation. Profiling could
   identify whether the FFBS or SV updates dominate, and there's
   room for further speed gains.

## Session update (2026-05-26 late) — Gibbs port done end-to-end on real data; 18× speedup with caveat

**The Gibbs port now works end-to-end on production-scale real ABS data.**
Steps 8, 9a, and 10 of [`gibbs_port_plan.md`](gibbs_port_plan.md) closed
this session. 85+ passing tests across `tests/testthat/gibbs/`.

**Headline numbers (T=435, N=11):**

| | Stan Variant A | Gibbs | Note |
|---|---|---|---|
| Wall time | **93 min** | **5.2 min** | **~18× faster** ✓ |
| rho posterior median | ~0.70 | **0.21** [-0.12, 0.51] | FFBS-noise bias |
| Trend latest-t median | 0.647 | 1.048 | Diverges with rho |
| Trend cor across t | — | 0.685 | Shape similar |
| Trend MAE | — | 0.486 | Level differs |
| Lambda medians | various | 0.91–1.03 | Matches Stan pattern |

**Production recommendation:** Stan stays headline; the monthly cron
and dashboard are **unchanged**. Gibbs is a fast prototyping /
sensitivity tool, not a production estimator — until the rho bias is
fixed via marginal MH (mitigation #1 in `gibbs_port_plan.md`,
estimated 1-2 days).

**Files added this session:**
- [`R/gibbs/mct_gibbs_fit.R`](R/gibbs/mct_gibbs_fit.R) — S3 wrapper mirroring cmdstanr API
- [`R/gibbs/build_ssm_mixed.R`](R/gibbs/build_ssm_mixed.R) — augmented 3(N+1)-dim SSM for mixed-freq
- [`R/gibbs/gibbs_sweep_mixed.R`](R/gibbs/gibbs_sweep_mixed.R) — mixed-freq Gibbs sweep
- [`scripts/step10_gibbs_real_data.R`](scripts/step10_gibbs_real_data.R) — real-data benchmark
- [`scripts/diagnose_rho_bias.R`](scripts/diagnose_rho_bias.R) — bias diagnostic
- [`scripts/diagnose_gibbs_mixing.R`](scripts/diagnose_gibbs_mixing.R) — mixing investigation

**Next session targets** (in order of expected value):

1. **rho bias fix via marginal MH** — propose rho from N(current, prop_sd),
   accept via p(y|rho) integrating out c with the Kalman filter. Should
   recover the true posterior of rho and fix the trend-level gap to Stan.
   Effort: 1-2 days. Path documented in `gibbs_port_plan.md` Step 7
   update block.

2. **Wire the Gibbs path into `_targets.R`** as an alternative variant
   ("Ag" for Gibbs?) so it can be invoked via `tar_make()` and the
   dashboard. Only worth doing AFTER the bias fix.

3. **Validate via the existing LOO comparison helper** in
   `R/compare_variants.R`. Should be a drop-in given the
   mct_gibbs_fit S3 wrapper exposes `$draws(variables = "log_lik")`.

## Session update (2026-05-26 evening) — Gibbs port Step 7 done; end-to-end works

`fit_mct_gibbs()` works end-to-end on simulated MCT Variant A data
(54 passing tests). Three new files:

- [`R/gibbs/build_ssm.R`](R/gibbs/build_ssm.R) — SSM constructor
- [`R/gibbs/gibbs_sweep.R`](R/gibbs/gibbs_sweep.R) — one full Gibbs pass
- [`R/gibbs/fit_mct_gibbs.R`](R/gibbs/fit_mct_gibbs.R) — main fit entry

**Status:** the sampler runs and recovers signal but **mixes slowly**.
With 500 burn + 1000 draws on T=200 N=5 sim, rho posterior median
~0.4-0.5 vs truth 0.7 (the AR(1) common factor competes with RW
sector trends for signal). Tests use loose "is the sampler working"
thresholds; production fits will need longer chains AND/OR a reparam.
See [`gibbs_port_plan.md`](gibbs_port_plan.md) Step 7 update for
mitigations to try in Step 9.

**Next session targets:**
- Step 8: `mct_gibbs_fit` S3 class so it plugs into existing
  downstream code (postprocess, dashboard, etc.) without changes
- Step 9: mixed-frequency observation handling + try center-rotation
  of sector trends to break the c-vs-s additive identification (or
  longer chains, or smarter init) — fix the slow rho mixing
- Step 10: end-to-end validation against the Stan Variant Ac fit on
  real ABS data

## Session update (2026-05-26 morning) — Gibbs port Steps 1-6 done

Steps 1-6 of [`gibbs_port_plan.md`](gibbs_port_plan.md) ported and
tested. R code at [`R/gibbs/`](R/gibbs/), tests at
[`tests/testthat/gibbs/`](tests/testthat/gibbs/) (**42 passing**).

| Step | Port | Status |
|------|------|--------|
| 1 | `simulate_SSM.m` → `sim_ssm.R` | ✅ |
| 2 | `Kalman_filter.m` → `kalman.R` | ✅ |
| 3 | `fast_smoother.m` + `simulation_smoother.m` → `fast_smoother.R` + `sim_smoother.R` | ✅ (MSE recursions skipped — not needed by Gibbs) |
| 4 | `update_vol.m` → `update_vol.R` (Kim-Shephard SV) | ✅ |
| 5 | `update_gam.m` → `update_gam.R` (hyperparams) | ✅ (`update_theta.m` for MA(q) skipped — we don't use MA errors) |
| 6 | `update_tvcoef.m` → `update_tvcoef.R` (TVP loadings) | ✅ |
| 7 | `estimate_MCT.m` → main Gibbs loop | ⏳ pending (BIG — 487 lines) |
| 8 | `mct_gibbs_fit` S3 class for downstream compatibility | ⏳ pending |
| 9 | AR(1) + mixed-frequency adaptation for AU | ⏳ pending |
| 10 | End-to-end validation vs Variant Ac Stan fit | ⏳ pending |

Critical correctness validations in tests:
- KF log-likelihood matches analytic Gaussian density (with + without missing data)
- sim_smoother draws average to fast_smoother point (validates posterior draws)
- Kim-Shephard SV recovers truth path with no bias (|log bias| < 0.20, r > 0.5)
- update_gam concentrates around truth as T grows (T=10K → <2% error)
- update_tvcoef recovers TVP path with r > 0.8

**Investigation also done:** [`investigation_ar1_vs_rw.md`](investigation_ar1_vs_rw.md)
explains why AR(1) beats RW in our model. Likely a regime-change artifact
of the 1990-96 disinflation. Four follow-up experiments proposed.

## Session update (2026-05-25/26) — speedup exploration

Full details in [`plan.md`](plan.md) (Phase 1-3 experiments) and
[`gibbs_port_plan.md`](gibbs_port_plan.md) (Phase 4 / Gibbs port plan).

**TL;DR:** explored 4 speedup variants (constant-λ Ac, threaded Act,
KF-marginalized Akf, sparse-KF Akfs). All math validated; posteriors
match the baseline Variant A within MC noise. **None of them beat A in
production-equivalent settings on this M1 machine.** Best result:
Akfs ≈ 4× faster than dense Akf but still ~8× slower than Ac. The
binding constraint is Stan's autodiff overhead on matrix-shaped state.

**Production unchanged.** Cron still uses Variant A. Local
`tar_make()` still uses Variant A. The four experimental Stan
variants exist but aren't wired into the targets pipeline.

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
