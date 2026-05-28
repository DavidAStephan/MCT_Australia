# Next session plan — MCT_Australia

Hand-off note for the next working session. Pair with
[`CLAUDE.md`](CLAUDE.md) (project brief),
[`README.md`](README.md) (current decisions + headline findings),
[`gibbs_port_plan.md`](gibbs_port_plan.md) (Gibbs port details), and
the auto-memory at
`~/.claude/projects/-Users-davidstephan-Documents-MCT-Australia-MCT-Australia/`.

## TL;DR — what to do next

**Investigate the constant-λ trend/common-share gap between Stan Ac
and Gibbs.** Step 14/15 (2026-05-26 afternoon, this session) ran:
(a) baseline Gibbs vs +outliers LOO comparison — baseline wins
**(elpd_diff = −3911, SE 883, 4.4σ; outliers HURT fit on AU data)**;
(b) Stan Variant Ac refit — the previously-suspected "TVP-vs-constant-λ"
explanation for the Gibbs-vs-Stan trend gap is **falsified**. Stan Ac
(constant λ) gives essentially the same trend as Stan A (TVP λ),
both ~0.65 latest-t. Gibbs gives 1.02 latest-t. Common-share is the
shocker: Stan Ac=0.069 vs Gibbs=0.915.

`rho` agrees across both samplers (Stan Ac 0.376 vs Gibbs 0.359 — the
AR(1) machinery is correct). So the model-spec difference is in the
variance decomposition, not the persistence parameter. Suspect (a) a
common_share formula divergence between `mct_gibbs_fit.R`'s
`.compute_derived_gq` and Stan's `generated quantities`, or
(b) a sigma_c / sigma_s prior or scaling discrepancy. Diff the two
formulas first, then compare posterior sigma_c and sigma_s draws.

## Concrete next actions (priority order)

### 0. Investigate trend / common-share gap (Stan Ac vs Gibbs)

The new finding from this session: with the Stan Ac (constant λ) refit
done, the trend-level gap and common-share gap can no longer be blamed
on the TVP-vs-constant λ difference. The rho posteriors match, so it's
not a sampler bug in the AR(1) update. Concrete debugging steps:

```r
# Load both
for (f in list.files("R",       "\\.R$", full.names = TRUE)) source(f)
for (f in list.files("R/gibbs", "\\.R$", full.names = TRUE)) source(f)
sd <- readRDS("outputs/draws/stan_data.rds")
fit_Ac  <- readRDS("outputs/draws/fit_Ac.rds")
fit_gib <- mct_gibbs_fit(readRDS("outputs/draws/fit_gibbs_baseline.rds"),
                         sd)

# 1. Compare posterior sigma_c[T] and sigma_s[T, i]
sc_Ac  <- as.numeric(fit_Ac$draws(paste0("sigma_c[", sd$T, "]"),
                                  format = "matrix"))
sc_gib <- as.numeric(fit_gib$draws(paste0("sigma_c[", sd$T, "]"),
                                   format = "matrix"))
# If sc_Ac / sc_gib differ by ~10x, that's where the common-share gap is.

# 2. Diff common_share formulas:
#    Stan: see stan/mct_aus_Ac.stan generated_quantities block
#    Gibbs: R/gibbs/mct_gibbs_fit.R, .compute_derived_gq()
```

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

| | Stan Variant A (TVP λ) | Stan Variant Ac (const λ) | Gibbs (const λ, no outliers) | Gibbs +outliers |
|---|---|---|---|---|
| Wall time | 93 min | 103 min | **5.4 min** | 5.3 min |
| rho posterior median | 0.371 [0.18, 0.54] | 0.376 [0.18, 0.55] | 0.359 [0.17, 0.53] | (similar) |
| Trend latest-t median | 0.647 | 0.652 | 1.023 | — |
| Common-share latest-t | (similar) | 0.069 | 0.915 | — |
| LOO vs Gibbs baseline | — | — | (best of the two below) | elpd −3911, SE 883 |

The trend/common-share gap between Stan and Gibbs is NOT a TVP-vs-const
artifact (Stan A and Ac agree). See action #0 for the investigation
plan. rho matches across all three so the AR(1) sampler is correct.

## Session log (recent → older, compressed)

### 2026-05-26 afternoon — Steps 14 + 15: outlier benchmark + Ac refit

**Step 14: real-data outlier benchmark.**
Ran baseline Gibbs (no outliers, 5.4 min) and +outliers Gibbs (5.3 min)
on cached AU data. Initial outlier run had a bug — 99.9% of obs flagged
because the mixed-freq path's update_ps was counting missing/quarterly
cells as "outliers", collapsing ps from prior 0.975 to ~0.19. **Fixed
in `R/gibbs/gibbs_sweep_mixed.R`**: v1 outlier mixture now applies to
monthly obs only (quarterly + missing cells get s_outlier=1; ps update
excludes them via NA mask). Re-run: 6/1698 obs flagged (0.4%), ps stays
near 0.97 across sectors. Flags are sensible: Health 2025-01/2026-01,
Transport 2026-03, Education 2025-02/2026-02 — all monthly tail spikes.

Also extended `R/gibbs/mct_gibbs_fit.R` `.compute_derived_gq()` to
compute log_lik for quarterly observations (avg-of-3-months mean +
1/9-of-sum-of-3 variance, matching `build_mct_ssm_mixed`'s observation
row for `ot==2`). Without this, LOO would only score the 253 monthly
obs; with this, all 1698 obs contribute. With-outlier path also
multiplies sigma_eps by s_outlier in the log_lik computation.

**LOO comparison (Step 14b):** baseline preferred by **elpd_diff = −3911
(SE 883, |Δ|/SE = 4.43)**. Outliers HURT fit on AU data — too few real
outliers to justify the extra parameters. p_loo nearly doubled (3203
→ 7772). Pareto-k diagnostics flag ~30-50% of obs as bad in both fits
(state-space LOO is inherently high-k; treat the comparison as
indicative). Combined-MA+mixed-freq stretch goal is now deferred per
the plan's own gating logic ("Defer unless item 2 shows outliers
aren't enough").

**Step 15: Stan Variant Ac refit (~103 min, 9/6000 divergences — fine).**
The Ac fit was deleted in the May speedup exploration; regenerated as
`outputs/draws/fit_Ac.rds` for the apples-to-apples Gibbs comparison.
Result is genuinely surprising:

|                       | Stan Variant A (TVP λ) | Stan Variant Ac (const λ) | Gibbs (const λ)  |
|---                    |---                     |---                        |---               |
| rho median            | 0.371 [0.18, 0.54]     | **0.376 [0.18, 0.55]**    | 0.359 [0.17, 0.53] |
| latest-t trend median | 0.647                  | **0.652**                 | 1.023            |
| latest-t common-share | (similar)              | **0.069**                 | 0.915            |

Stan A and Stan Ac give essentially the same trend (~0.65); Gibbs
gives 1.02. So the trend-level gap is **not** the TVP-vs-constant-λ
model spec — it's something else. The common-share gap (0.07 vs 0.92)
is the bigger smoking gun. rho matches across all three, so the AR(1)
machinery is correct. Suspect a `common_share` formula divergence or a
sigma_c / sigma_s prior/scaling difference. See action #0 above.

Files written:
`outputs/draws/fit_gibbs_baseline.rds`,
`outputs/draws/fit_gibbs_outliers.rds`,
`outputs/draws/fit_Ac.rds`,
`outputs/draws/loo_baseline_vs_outliers.rds`,
`outputs/tables/outlier_top5_per_sector.csv`.
Scripts: `scripts/step14_outlier_benchmark.R`,
`scripts/step14b_loo_baseline_vs_outliers.R`,
`scripts/step15_refit_Ac.R`, `scripts/step15b_compare_Ac_vs_gibbs.R`.

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
