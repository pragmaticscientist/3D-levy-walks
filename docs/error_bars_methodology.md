# Uncertainty Bands in the Published Figures

This note documents how the uncertainty bands in the figures under
`plots/pdf_figures/` are computed. All of the relevant code lives in
`data-analysis.ipynb`; this document describes what that code does and records
the known limitations of the current choice of estimator.

## 1. Summary

Every band is a **±1 standard error of the mean (SEM)** interval, rendered as a
shaded region (`fill_between`, `alpha=0.2`, no edge line) rather than as capped
error bars. No confidence level is claimed and no distributional assumption is
stated in the figures themselves.

The SEM is obtained directly from the sample of independent simulation runs in
each plotted bin:

```
sem = std / sqrt(count)
```

where `std` is the pandas default (Bessel-corrected, ddof=1) sample standard
deviation of `detection_time` over the runs in the bin, and `count` is the
number of such runs.

## 2. Generation procedure

Each plotting function is invoked twice: once with `error_bands=False` and once
with `error_bands=True`. The two calls write to distinct files, the second
carrying an `_errorbands` suffix — for example

```
plots/pdf_figures/detection_time_vs_surface.pdf
plots/pdf_figures/detection_time_vs_surface_errorbands.pdf
```

Consequently, **whether a given figure carries uncertainty bands is determined
by which of the two files is used downstream**, not by a property of the
underlying analysis. Both files are always produced.

The common pipeline is:

1. Filter the raw per-run records to the fixed parameters of the panel
   (`n_walkers`, `n_targets`, `n_volume`, `probability`, `TargetShape`, …).
2. Group by the plotted independent variable and any distinguishing series
   keys.
3. Aggregate `detection_time` to a central statistic together with `std` and
   `count` in a single `.agg()` call.
4. Form `sem = std / sqrt(count)`.
5. Shade `centre - sem` to `centre + sem`.

## 3. Direct plots

For figures whose ordinate is a detection time rather than a ratio, the central
statistic is the **mean** and the band is `mean ± sem`.

| Figure family | Function | Grouping |
| --- | --- | --- |
| `detection_time_vs_surface` | `plot_detection_vs_surface` | `surface` |
| `detection_time_vs_mu_ball` | `plot_detection_vs_mu_for_all_surface` | `mu`, `surface`, `TargetShape`, `n_walkers` (+ fixed keys) |
| `detection_time_vs_delta_by_mu`, `…_small` | `plot_detection_time_vs_delta_by_mu_facets` | `mu`, `delta`, `D` |

## 4. Ratio plots

For figures whose ordinate is a ratio of two detection times, the two SEMs are
propagated through the quotient by first-order (delta-method) error
propagation, i.e. relative errors added in quadrature:

```python
rel_err = np.sqrt((sem_num / mean_num)**2 + (sem_den / mean_den)**2)
err = ratio * rel_err
```

This treats numerator and denominator as independent, which holds for the
Ball/Line comparison (separate simulation runs) but **not** for the ratio to the
μ=2 baseline in `plot_detection_ratio_vs_mu_for_all_surface`, where the same
μ=2 estimate is reused as the denominator across the whole curve. The bands on
such a curve are therefore correlated, and the μ=2 point is fixed at exactly 1
by construction.

| Figure family | Function | Ratio |
| --- | --- | --- |
| `detection_time_ratio_vs_mu_{ball,disk,line}` | `plot_detection_ratio_vs_mu_for_all_surface` | t(X^μ) / t(X^2), per shape and surface |
| `detection_time_ratio_ball_line`, `…_surface` | `plot_detection_time_ratio` | t_Ball / t_Line at matched geometry |

## 5. Known limitations

These are recorded so that the bands are not over-interpreted.

**5.1 Mismatched estimator in `plot_detection_time_ratio` — resolved.**
That function previously aggregated the central statistic as a **median** while
still building the band from `std / sqrt(count)`, the standard error of the
*mean*. The central statistic is now the mean, exposed as a `statistic`
parameter defaulting to `"mean"`, so the band matches the plotted quantity.
Passing `statistic="median"` restores the former behaviour, with the caveat
above.

The change is cosmetic for the central curve: with ~10,000 runs per bin the
ratio of means and the ratio of medians agree to within about 2% at μ = 1, 2
and 3, which is the scale of the run-to-run scatter. No conclusion drawn from
`detection_time_ratio_ball_line*.pdf` changes.

**5.2 Heavy tails.** Detection times under a Lévy walk with μ ≤ 2 can have
infinite or very slowly converging variance. Where that is the case the sample
`std` is not a stable quantity: it grows with `count` instead of settling, so
`std / sqrt(count)` neither converges nor bounds the error of the plotted
statistic in the usual way. This affects the low-μ end of every ratio-vs-μ
figure. The mean itself is subject to the same concern.

**5.3 Double-counted input in the fixed-surface figure — resolved.** The cell
producing `detection_time_ratio_ball_line_surface*.pdf` globbed
`results/detection_time_fixed_surface/*` and concatenated every file matched.
That directory holds both the aggregate `detection_time_fixed_surface.csv` and
the 50 per-thread shards it was built from; the shards were verified to
reproduce the aggregate exactly (1,540,000 rows, identical line sets). Every
sample was therefore loaded twice. Means and medians are unaffected by exact
duplication, but `count` doubled, so the published bands were narrower than
warranted by a factor of sqrt(2). The cell now reads the aggregate CSV alone.

Any other analysis cell that globs a results directory is exposed to the same
hazard wherever aggregate and shard files are stored side by side.

**5.4 No confidence level.** A ±1 SEM band corresponds to roughly 68% coverage
only under approximate normality of the sampling distribution of the estimator.
That approximation is weakest exactly where the tails are heaviest.

## 6. Suggested replacement

For any figure intended for publication, a nonparametric bootstrap over the
per-run samples — resampling within each bin and taking percentile intervals of
the *plotted statistic* (mean, median, or ratio as appropriate) — would remove
both the estimator mismatch of §5.1 and the reliance on a finite variance in
§5.2, at the cost of one resampling pass per bin. This has not been implemented.
