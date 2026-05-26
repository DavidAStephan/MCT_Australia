# Next session plan — MCT_Australia

Hand-off note for the next working session. Pair with
[`CLAUDE.md`](CLAUDE.md) (project brief),
[`README.md`](README.md) (current decisions + headline findings),
[`gibbs_port_plan.md`](gibbs_port_plan.md) (Gibbs port details), and
the auto-memory at
`~/.claude/projects/-Users-davidstephan-Documents-MCT-Australia-MCT-Australia/`.

## TL;DR — what to do next

**Run the outlier-enabled Gibbs on real AU CPI data** and check whether
the model flags known volatile events (fuel spikes 2022-23, COVID
disruptions 2020). Then compare to the baseline (no-outlier) Gibbs via
LOO to settle whether outliers materially improve fit on AU data.

## Concrete next actions (priority order)

### 1. Real-data outlier benchmark (~30 min compute, ~15 min review)

Fit the Gibbs sampler with outliers ON to the cached AU data. Goal:
verify the outlier path works at production scale AND empirically
flags real-world outliers in the AU CPI series.

```r
# In project root, R session:
suppressMessages(library(MASS))
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)

sd_real <- readRDS("outputs/draws/stan_data.rds")
T_ <- sd_real$T; N <- sd_real$N

# Build y/obs_type from the sparse stan_data lists
y        <- matrix(NA_real_, T_, N)
obs_type <- matrix(0L,      T_, N)
for (k in seq_along(sd_real$y_m)) {
  y[sd_real$t_m[k], sd_real$i_m[k]] <- sd_real$y_m[k]
  obs_type[sd_real$t_m[k], sd_real$i_m[k]] <- 1L
}
for (k in seq_along(sd_real$y_q)) {
  y[sd_real$t_q[k], sd_real$i_q[k]] <- sd_real$y_q[k]
  obs_type[sd_real$t_q[k], sd_real$i_q[k]] <- 2L
}

fit_out <- fit_mct_gibbs(
  y, obs_type = obs_type, ref = sd_real$ref,
  n_burn = 1000L, n_draw = 2000L, verbose = TRUE,
  config = list(use_outliers = TRUE,
                use_marginal_mh_rho = FALSE)   # MH not needed, saves ~50% time
)

saveRDS(fit_out, "outputs/draws/fit_gibbs_outliers.rds")
```

Then inspect posterior median `s_outlier`:

```r
s_med <- apply(fit_out$draws$s_outlier, c(1, 2), median)
# For each sector, find the top-5 flagged dates (highest s_med)
for (i in seq_len(N)) {
  top <- order(-s_med[, i])[1:5]
  cat(sprintf("%-30s ", sd_real$groups[i]))
  cat(paste(format(sd_real$dates[top], "%Y-%m"),
            sprintf("(s=%.1f)", s_med[top, i])),
      sep = "  ")
  cat("\n")
}
```

**What to look for:**
- Transport / Recreation / Holiday travel should flag fuel spikes
  around 2022-04 (Russia-Ukraine), 2008-10 (GFC oil), 2022 broadly.
- Food / Alcohol may flag COVID supply-chain dates (2020-04, 2020-Q3).
- A sector with NO flagged dates suggests its SV process is capturing
  everything — that sector doesn't need the outlier model.
- Most monthly obs (24/month × 11 sectors) should have s_med = 1.
  If many obs are flagged, prior `ps_a_prior = 97.5` may be too loose.

Wrap this into `scripts/step14_outlier_benchmark.R` so it's reproducible.

### 2. LOO comparison: baseline Gibbs vs +outliers (~5 min)

If item 1 produces sensible flags, compare LOO between the two Gibbs
fits to quantify the model-fit improvement:

```r
fit_base    <- readRDS("outputs/draws/fit_gibbs_baseline.rds")  # no outliers
fit_outlier <- readRDS("outputs/draws/fit_gibbs_outliers.rds")
# Both expose $draws("log_lik"); plug into loo::loo_compare
ll_base    <- fit_base$draws(variables = "log_lik")
ll_outlier <- fit_outlier$draws(variables = "log_lik")
# ... loo() each, then loo_compare()
```

A baseline Gibbs fit (`use_outliers = FALSE`) isn't cached yet — need
to run that too in item 1 or as a step 0. Document the elpd_diff +
SE in the commit message.

### 3. (Stretch) Extend MA path to mixed-freq

The MA(q) path is currently monthly-only. Combined with mixed-freq it
would unlock MA on the full T=435 AU sample (vs the T=23 monthly tail
today). State dim ~92 at N=11, q=3; needs careful implementation but
no new conceptual work. Defer unless item 2 shows outliers aren't
enough.

### 4. (Stretch) Stan Variant Ac refit for the apples-to-apples bench

Per Step 11 notes, the trend-level gap between Stan Variant A (TVP λ)
and Gibbs (constant λ) is a model-spec difference, not a sampler bug.
The fair Stan comparison is against Variant Ac (constant λ), but the
Ac fit was deleted in the May speedup-exploration. Re-fitting Ac (~5
min on M1) closes that validation question.

## Current project state

**Production (unchanged):**
- Stan Variant A on monthly GitHub Actions cron — last successful run
  2026-05-24, dashboard live at github.io page.
- Local `targets::tar_make()` still uses Variant A.

**Gibbs port** (14 R files in [`R/gibbs/`](R/gibbs/), **120 tests
passing**):

| Capability | Files | Status |
|---|---|---|
| Generic SSM + KF + smoothers | `sim_ssm.R`, `kalman.R`, `fast_smoother.R`, `sim_smoother.R` | ✅ |
| Hyperparam + SV updates (Kim-Shephard) | `update_vol.R`, `update_gam.R`, `update_tvcoef.R` | ✅ |
| MCT-specific spec, monthly-only | `build_ssm.R`, `gibbs_sweep.R`, `fit_mct_gibbs.R` | ✅ |
| Mixed-frequency obs | `build_ssm_mixed.R`, `gibbs_sweep_mixed.R` | ✅ |
| Marginal-MH rho update | `update_rho_mh.R` | ✅ |
| MA(q) measurement errors (monthly-only) | `update_theta.R`, `build_ssm_ma.R`, `gibbs_sweep_ma.R` | ✅ |
| Outlier scale mixture (mixed-freq compat) | `update_scl.R`, `update_ps.R`, sweep updates | ✅ |
| cmdstanr-style S3 wrapper | `mct_gibbs_fit.R` | ✅ |
| Combined MA + mixed-freq | — | ⏳ deferred |

**Headline numbers (real AU data, T=435, N=11):**

| | Stan Variant A | Gibbs (no MA, no outliers) |
|---|---|---|
| Wall time | 93 min | **5.2 min** (4-6× faster) |
| rho posterior | 0.371 [0.18, 0.54] | 0.360 [0.17, 0.53] — **identical** |
| Trend latest-t | 0.647 | 1.059 (model-spec diff: TVP vs constant λ) |

## Session log (recent → older, compressed)

### 2026-05-27 evening — Steps 12 + 13: MA(q) + outliers

Both NY Fed methodology gaps closed. MA(q) needs state augmentation
((q+1) eps lags per sector) and is currently monthly-only. Outliers
needed only per-obs scalars; integrates cleanly with mixed-freq.
14 new tests; total 120 passing. Commits: `612dd14` (MA), `1276281`
(outliers).

### 2026-05-27 morning — Step 11: validated, "rho bias" was faulty premise

Marginal MH for rho implemented. Long-chain comparison vs Stan
showed rho posteriors are statistically identical (Stan 0.371 vs
Gibbs 0.360); the earlier "Stan rho ~0.7" claim was misremembered.
Trend-level gap remains but it's a TVP-λ-vs-constant-λ model-spec
difference, not a sampler bug. Commit: `0013ad8`.

### 2026-05-26 — Gibbs port Steps 1-10 (3 sessions)

Full Gibbs port from NY Fed MATLAB to R/gibbs/. Steps 1-6 (simulator
through TVP loadings, morning), Step 7 (main loop + first
end-to-end, evening), Steps 8-10 (S3 wrapper, mixed-freq, real-data
benchmark, late). 91 tests passing at end. Production untouched.
Per-iter cost on real data: ~0.1 sec (smoke test) → 5.2 min for full
chain.

### 2026-05-25/26 — Stan-side speedup exploration

Built 4 experimental Stan variants (const-λ Ac, threaded Act, KF-marg
Akf, sparse-KF Akfs). All math-validated; **none beat A in production
config**. Best (Akfs) was 4× faster than the dense KF but still 8×
slower than Ac. Binding constraint identified as Stan's autodiff
overhead on matrix-shaped state — motivated the Gibbs port.
Full results in [`plan.md`](plan.md).

### 2026-05-24 — Cron + dashboard

Dashboard renders end-to-end. GitHub Actions monthly cron wired up,
first successful auto-publish that day. `_targets.R` refactored so
fits don't duplicate into `_targets/objects/` (saved ~7 GB). Disk
shrank from ~14.5 GB to ~6.5 GB.

### 2026-05-23 — Initial 10-step model build complete

All 10 brief steps shipped (data, prep, Stan models, sim test,
real-data fits, postprocess, LOO comparison, targets pipeline,
dashboard, GH Actions workflow). Variant C added on top.
3-way LOO: **A ≻ B ≻ C** by 3σ and 6σ respectively.

## Standing v2 backlog (parked items)

- **MA + mixed-freq combined** — see #3 above. State dim ~92 at
  N=11, q=3. Needed before MA can be used on AU production data.
- **TVP λ in Gibbs** — already ported (`R/gibbs/update_tvcoef.R`)
  but not wired into the sweep. Would close the trend-level gap
  to Stan Variant A.
- **Wire Gibbs into `_targets.R`** as alternative fit_Ag target.
  Drop-in via the `mct_gibbs_fit` S3 wrapper.
- **Vintage-accurate weights** (pre-2017 ABS used 5-yearly weight
  links). Robust to ±2pp shifts in the trend story.
- **Sub-group panel** (~30-50 series) closer to NY Fed's 17 PCE
  spirit. Gibbs's speed advantage compounds at higher N.
- **K-fold CV via `loo::kfold()`** to dodge the high Pareto-k
  issue in state-space LOO (132 / 113 / 76 high-k obs for A/B/C).

## Headline empirical finding (kept for context)

Australian sectoral CPI co-movement is **overwhelmingly transitory**,
NOT a drifting common trend. Opposite to NY Fed's finding for US
PCE. Variant A (AR(1) common) decisively beats B (RW common) on LOO
by 3σ.

The agent investigation [`investigation_ar1_vs_rw.md`](investigation_ar1_vs_rw.md)
flagged this is likely a regime-change artifact of the 1990-96
disinflation (a one-off mean-reverting episode that AR(1) captures
naturally). NY Fed prefers RW because their sample is longer + wider
+ has MA(3) + outlier scaffolding that absorbs short-run shocks our
model attributes to the AR(1) cycle.

Latest-month medians (March 2026, demeaned annualised %, headline=A):
trend median 0.65, 68% band [-0.08, +1.41]. Common-share (variance
ratio) median ~7%. Add 2.5pp for raw %.

## Quick orientation commands

```sh
# Resume — branches/commits
git log --oneline | head
git status

# Reload cached objects without re-fitting
Rscript -e '
  for (f in list.files("R", "\\.R$", full.names = TRUE, recursive = TRUE))
    source(f)
  stan_data <- readRDS("outputs/draws/stan_data.rds")
  fit_A     <- readRDS("outputs/draws/fit_A.rds")
  # ready to call trend_path(), common_share_path(), etc.
'

# Run the full Gibbs test suite (~3-5 min)
Rscript -e 'testthat::test_dir("tests/testthat/gibbs")'

# Run the full Stan-side test suite (fast bits only)
Rscript -e 'testthat::test_dir("tests/testthat")'

# Re-run the 3-way LOO comparison on cached Stan fits
Rscript -e '
  for (f in list.files("R", "\\.R$", full.names = TRUE)) source(f)
  loo_A <- loo_for_fit(readRDS("outputs/draws/fit_A.rds"))
  loo_B <- loo_for_fit(readRDS("outputs/draws/fit_B.rds"))
  loo_C <- loo_for_fit(readRDS("outputs/draws/fit_C.rds"))
  print(compare_fits_loo(list(A = loo_A, B = loo_B, C = loo_C)))
'

# Full production pipeline (Stan-side)
Rscript -e 'targets::tar_make()'
```
