# Why does Variant A (AR(1)) beat Variant B (RW) on Australian CPI?

LOO: `A = 0.0`, `B = −13.0` (SE 4.4) — a ~3σ rejection of B. NY Fed's analogous
US PCE model uses RW. This note tries to attribute the divergence and proposes
specific follow-ups.

## 1. Theoretical contrast

- **RW common factor (B / NY Fed):** `c_t = c_{t−1} + ση_t`. Shocks are
  permanent. The cross-sector co-movement contributes to *trend* inflation.
  In our Variant B the headline trend formula is
  `τ_t = (Σ w_i λ_{i,t})·c_t + Σ w_i s_{i,t}` — `c_t` enters the trend.
- **AR(1) common factor (A):** `c_t = ρ c_{t−1} + ση_t`, `|ρ|<1`. Shocks die
  out at rate `1−ρ`. Co-movement is a *cycle*; trend is recovered purely from
  sector-specific RW drifts: `τ_t = Σ w_i s_{i,t}`.

These are observationally distinct only at low frequencies. At T=435 with SV
on `s_{i,t}` already absorbing slow drift sector-by-sector, the RW common
factor's permanent component is largely redundant — sector RWs can already
produce a common low-frequency direction whenever the loadings happen to
align. The AR(1) cycle then has unique work to do: it must capture
short-to-medium horizon cross-sector co-movement (energy/tradables shocks,
GST in 2000, COVID in 2020-22) that would otherwise hit each sector's idio
noise. With `sigma_lambda ~ N(0, 0.01)` the loadings are nearly time-constant,
so a single shared cycle is genuinely informative about contemporaneous
co-movement.

## 2. Why NY Fed picks RW for US PCE

Three structural differences favour RW for them:

1. **Sample length**: ~750 monthly obs from 1960 spans the Great Inflation,
   Volcker disinflation, Great Moderation, GFC, and post-COVID. The
   identifying variance in *long-run* common drift is huge. AU's 435 obs
   from 1990 has very little long-run drift left to identify once the
   disinflation completes by 1996.
2. **Sectoral granularity**: 17 PCE sub-indices vs our 11 ABS groups. More
   sectors give the RW common factor more cross-sectional traction (per-t
   signal-to-noise scales as √N).
3. **Noise scaffolding**: NY Fed has MA(3) errors and outlier scale mixture.
   Those absorb short-run idiosyncratic shocks that, in our model, hit
   measurement-error SV. Without that scaffolding the AR(1) cycle in A
   "picks up" what NY Fed's MA(3)+outliers would absorb, giving A an
   unfair advantage on LOO over B.

## 3. The 1996 regime change matters a lot

RBA adopted 2-3% inflation targeting in 1993; disinflation completed ~1996.
Pre-1996, Australian inflation drifted from ~8% → ~2%; post-1996 it has been
mean-reverting around target. **Our 1990-onward sample contains exactly one
big mean-reverting episode and then 30 years of stationarity.** Demeaning on
2000-2019 leaves the 1990-1995 segment with a large positive level deviation
that decays toward zero — *the textbook AR(1) signature*. An AR(1) cycle
with `ρ` around 0.85-0.95 fits this beautifully; a RW common trend has to
explain it as a single huge permanent shock that, awkwardly, then reverts.

So the AR(1) win is almost certainly **partly artifact** of fitting one big
1990-96 disinflation episode with a stationary cycle. That doesn't mean A is
wrong — it means we are not learning something deep about
US-vs-AU sectoral physics; we are learning that AU's sample is dominated by a
single regime change that looks like mean reversion.

## 4. Are sector RW trends already doing the work?

In B, `s_{i,t}` has its own SV-driven RW. The likelihood can put the 1990-96
disinflation entirely into the sector trends — leaving `c_t` with little
work and a large penalty for its extra parameters. The 13-elpd LOO gap could
plausibly be the marginal-likelihood penalty for an under-utilised RW common
factor, NOT evidence that the data wants AR(1) per se. This is consistent
with Variant C losing by 71 elpd to A: adding the RW factor on top of the
AR(1) cycle doesn't help because the sector RWs already eat its variance.

## 5. Other artifacts to rule out

- **Mixed-frequency design**: quarterly obs (1990-2024) have averaged
  measurement-noise variance `(Σ exp(h_e))/9`. Pre-2024 the effective
  observation noise is √3 smaller than monthly, which makes the model
  *over-confident* on the pre-2024 fit and may amplify the LOO weight of
  pre-2000 transitions — exactly the regime-change segment.
- **Demean window (2000-2019)**: this puts the 1990-96 disinflation
  ENTIRELY on the positive-residual side. A symmetric demean using
  1990-present would attenuate the AR(1) cycle's level signal.
- **`sigma_lambda ~ N(0, 0.01)`**: very tight. RW loadings can barely move,
  so the AR(1) cycle's interpretation is forced to be the same shock loading
  across the full sample, which suits a mean-reverting cycle better than a
  drifting trend (since a drifting trend with fixed loadings would imply a
  proportional drift across sectors that AU data does not show).

## Recommended follow-ups (ranked)

1. **Re-fit A and B on post-1996 only (T≈360).** If the elpd_diff collapses
   toward zero or flips, the AR(1) win is a regime-change artifact, not a
   structural property of AU sectoral co-movement. Cheap (~1.2× the current
   single-fit cost, since post-1996 only drops ~17% of obs).
2. **A with a deterministic break dummy at 1996Q1** (or a level-shift on
   `s_init`). Compare elpd to vanilla A. If the break absorbs most of the
   AR(1)'s explanatory load, `ρ` posterior should fall sharply and A vs B
   should narrow. This is the cleanest test of "AR(1) = regime transition
   in disguise".
3. **Re-fit on the post-2000 demean window only (T=315).** Eliminates the
   pre-demean-window segment entirely; if A still beats B by >2σ here, the
   AR(1) story has real legs even within the inflation-targeting regime.
4. **Tighten `rho` prior to favour permanence** (e.g. `ρ ~ N(0.98, 0.02)`
   truncated to (0,1)). If A's elpd is insensitive, the data identifies `ρ`
   away from 1 on its own; if it drops sharply, the RW limit is poorly
   identified rather than rejected — which would weaken the "AU is
   transitory" claim.

A clean experiment-1+2 combination would let us write: "Variant A beats B
mainly because the 1990-96 disinflation is mean-reverting; conditional on a
1996 break dummy and on post-1996 data, the AR(1) vs RW comparison is
[X]". That is a defensible and substantive finding either way.
