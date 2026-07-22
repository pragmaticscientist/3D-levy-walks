# Empirical trajectory test of a finite Lévy-walk surrogate in three immune-cell conditions

## Abstract

We tested whether a finite flat-to-power-law relocation model can reproduce
search statistics derived from measured three-dimensional immune-cell
trajectories. Ground-truth tracks were obtained from the Leukocyte Tracking
Database (LTDB), segmented into directed runs at large changes in direction,
and fitted separately by biological and sampling condition. The main-text
analysis retained a prespecified natural-killer (NK) cell condition and the two
remaining eligible conditions with the largest numbers of complete directed
runs: influenza-associated NK cells in popliteal lymph node (27 tracks, 517
runs), influenza-associated neutrophils in popliteal lymph node (224 tracks,
5,829 runs), and HIV-model T cells in popliteal lymph node (163 tracks, 3,899
runs). Every complete run received equal likelihood weight.

The fitted effective exponents were \(2.744\), \(2.655\), and \(2.560\),
respectively. Formal goodness of fit was acceptable only for the NK condition;
the model is therefore interpreted as an effective surrogate rather than
evidence for an exact Lévy law. We replayed each measured track against Ball,
Disk, and Line targets in a periodic three-dimensional domain and compared five
finite trajectory constructions. A single random rotation of an intact track
did not systematically change detection. In contrast, independently rotating
every directed-run vector while retaining the exact track's ordered lengths and
travel budget increased detection in all 27 condition–target settings and
significantly increased it in 21. Subsequent replacement by runs sampled from
the condition pool produced a smaller change. The fitted model closely
reproduced this independently oriented empirical process in 26 of 27 finite
settings and 22 of 27 unbounded first-passage settings. Thus the marginal
run-length model transfers well to randomized empirical searches in these
examples, while directional structure in intact tracks remains consequential.

## 1. Study design

The complete analysis was run for all 12 LTDB conditions containing at least
100 complete directed runs. To obtain figures that fit in the main paper
without selecting conditions according to their outcomes, we used the
following nested rule:

1. retain the NK condition chosen before the cross-condition comparison; and
2. add the two remaining eligible conditions with the largest complete-run
   samples.

This selected NK cells at 30 s, neutrophils at 15 s, and HIV-model T cells at
15 s. The complete 12-condition results, including all exclusions and
sensitivity analyses, remain in the
[all-condition report](equal_run_all_conditions_analysis.md) and should be
provided as Supplementary Information. The selected conditions are
illustrative empirical tests, not a controlled comparison of cell types:
sampling interval and experimental context still differ.

## 2. Dataset

### 2.1 Source and format

The trajectories come from the manually curated consensus ground truth in the
Leukocyte Tracking Database described by
[Pizzagalli *et al.* (2018)](https://doi.org/10.1038/sdata.2018.129); the
corresponding data collection is archived on
[Figshare](https://doi.org/10.6084/m9.figshare.c.3827890). LTDB combines
intravital two-photon microscopy experiments spanning several leukocyte types,
organs, stimuli, and acquisition intervals.

Each semicolon-delimited ground-truth file contains a header

```text
VideoID ; dx ; dy ; dz ; dt
```

followed by a channel-visibility row and records of the form

```text
TrackID ; x ; y ; z ; frame
```

The \(x,y,z\) coordinates are already in micrometres. Physical time is the
frame index multiplied by the file-specific interval \(dt\).

Across the full corpus, 44 biological files contained 44,722 annotations from
728 raw tracks. After quality control, 44,657 annotations from 710 tracks
remained, yielding 43,947 frame-to-frame displacements and 17,981 candidate
directed runs. The three selected conditions comprised 28,240 raw annotations
from 437 tracks. After preprocessing, 28,218 annotations from 431 tracks
remained; 414 tracks supplied at least one complete run and entered fitting and
replay.

| Condition | Contributing videos | Raw tracks | Replay tracks | Complete runs used for fitting | Exact-track run vectors |
|---|---:|---:|---:|---:|---:|
| NK · PLN · influenza · 30 s | 3 | 29 | 27 | 517 | 571 |
| Neutrophils · PLN · influenza · 15 s | 10 | 234 | 224 | 5,829 | 6,277 |
| T · PLN · HIV model · 15 s | 6 | 174 | 163 | 3,899 | 4,225 |
| **Total** | **19** | **437** | **414** | **10,245** | **11,073** |

“PLN” denotes popliteal lymph node. Exact-track vectors include the recorded
first and last runs; complete fitting runs exclude them as described below.

### 2.2 Preprocessing

Observations were sorted within tracks. Tracks with fewer than six
observations were removed, duplicate observations were prohibited, missing
frames split a track into separate fragments, and no interpolation or
trajectory smoothing was applied. For consecutive positions, the displacement
vector was

\[
\mathbf u_i=\mathbf x_{i+1}-\mathbf x_i ,
\]

and the three-dimensional turning angle was calculated as

\[
\theta_i=
\operatorname{atan2}
\!\left(
\lVert\mathbf u_{i-1}\times\mathbf u_i\rVert,\,
\mathbf u_{i-1}\cdot\mathbf u_i
\right).
\]

The `atan2` construction is numerically stable near both \(0\) and
\(\pi\). Zero-length displacements terminated a run and were not assigned a
direction. There was no positive speed cutoff in the primary analysis.

### 2.3 Directed-run definition and the \(90^\circ\) threshold

Consecutive non-zero displacements were concatenated until the angle between
successive vectors was strictly greater than \(90^\circ\). The displacement
at which the threshold was crossed started the next run. Each relocation
length was the end-to-end norm of its directed-run vector, rather than the
sum of the underlying frame displacements.

The \(90^\circ\) cutoff was a prespecified, literature-informed operational
definition. In an in-vivo neutrophil study,
[Georgantzoglou *et al.* (2022)](https://doi.org/10.1083/jcb.202103207)
classified changes below \(45^\circ\) as persistent, changes from
\(45^\circ\) to \(90^\circ\) as small turns, and changes from \(90^\circ\) to
\(180^\circ\) as large U-turns. Our rule therefore ends a directed run when
the change enters their U-turn category. Their study did not propose our exact
three-dimensional segmentation algorithm, so \(90^\circ\) should not be
described as a universal biological boundary. Sampling interval can itself
alter inferred speed–turning relationships
([Ganusov *et al.*, 2023](https://doi.org/10.1088/1478-3975/acb18c)).
Accordingly, conditions with different \(dt\) were not pooled, and
\(45^\circ\) and \(60^\circ\) rules were prespecified sensitivity settings.

The first and last run of each track were excluded from fitting because the
recording window can censor their full length. They were restored for finite
replay because they form part of the observed path and travel budget.

## 3. Run-length model and inference

For each condition, directed-run length \(\ell\) was fitted to

\[
f(\ell\mid\mu,c,L_{\max})
=
a(\mu,c,L_{\max})
\min\!\left[
1,\left(\frac{\ell}{c}\right)^{-\mu}
\right],
\qquad
0\leq\ell\leq L_{\max},
\]

where \(a\) normalizes the density, \(c\) is the crossover from the flat core
to the power-law tail, and \(L_{\max}\) is a finite upper cutoff. We set
\(L_{\max}\) to the condition-specific empirical maximum and estimated
\(\mu\) and \(c\) jointly by maximum likelihood.

Every complete run had equal likelihood weight,

\[
w_i=\frac{1}{N_k}
\]

within condition \(k\). Thus tracks contributing more complete runs also
contributed more likelihood. Biological clustering was retained in uncertainty
calculations: 95% intervals used 1,000 whole-track bootstrap resamples.
Goodness of fit used 500 parametric-bootstrap Kolmogorov–Smirnov replicates,
refitting \(\mu\) and \(c\) each time. Predictive comparisons against
truncated lognormal, exponential, Weibull, and cutoff-power-law alternatives
held out whole tracks in up to ten folds.

| Condition | Tracks / runs | \(\mu\) [95% CI] | \(c\), \(\mu\)m [95% CI] | \(L_{\max}\), \(\mu\)m | KS \(p\) | Best held-out model |
|---|---:|---:|---:|---:|---:|---|
| NK · PLN · influenza · 30 s | 27 / 517 | 2.744 [2.175, 3.645] | 2.332 [1.653, 2.941] | 40.043 | 0.0719 | truncated lognormal |
| Neutrophils · PLN · influenza · 15 s | 224 / 5,829 | 2.655 [2.533, 2.767] | 3.138 [3.070, 3.239] | 133.446 | 0.0020 | cutoff power law |
| T · PLN · HIV model · 15 s | 163 / 3,899 | 2.560 [2.378, 2.732] | 4.187 [4.134, 4.261] | 73.797 | 0.0020 | fitted mixture |

Only the NK fit was not rejected at the 5% level, and its held-out winner was
still the truncated lognormal. The fitted mixture ranked first for T cells,
but its held-out advantage over the cutoff power law was negligible. The three
values of \(\mu\) are consequently **effective exponents**, not evidence that
the underlying cells execute an exact Lévy walk.

Equal-run weighting also changed the estimand relative to the earlier
cell-balanced analysis. For NK, neutrophils, and T cells, respectively,
\(\mu\) changed from \(2.325\) to \(2.744\), \(2.475\) to \(2.655\), and
\(1.928\) to \(2.560\). Equal-run weighting answers what is observed after
selecting a complete run uniformly; cell-balanced weighting first selects a
track uniformly.

### 3.1 Segmentation sensitivity

Changing the angle threshold altered both the number of complete runs and the
effective exponent:

| Condition | \(45^\circ\): runs / \(\mu\) | \(60^\circ\): runs / \(\mu\) | \(90^\circ\): runs / \(\mu\) |
|---|---:|---:|---:|
| NK | 862 / 2.523 | 757 / 2.515 | 517 / 2.744 |
| Neutrophils 15 s | 9,462 / 3.247 | 8,339 / 3.037 | 5,829 / 2.655 |
| T HIV | 6,397 / 2.957 | 5,577 / 2.774 | 3,899 / 2.560 |

This sensitivity supports treating \(\mu\) as an operational summary. The
\(90^\circ\) rule was chosen from the prespecified biological interpretation,
not because it optimized fit. Search simulations were performed only for the
primary \(90^\circ\) segmentation.

## 4. Search simulations

### 4.1 Dimensionless geometry

Lengths in each condition were divided by its fitted crossover:

\[
\ell^*=\frac{\ell}{c},\qquad c^*=1.
\]

The detection radius was \(d=1\), and the cubic torus side was

\[
W=\frac{2L_{\max}}{c}.
\]

Targets were centred Ball, Disk, and Line objects. Their dimensions were
chosen so that the maximum face-on projected area of the corresponding
detection neighbourhood was

\[
A/c^2\in\{4,6,8\}.
\]

All 27 selected condition–shape–area targets fit without periodic
self-overlap. A common normalized area is not a common physical area:
\(A_{\mu\mathrm m^2}=Ac^2\). Likewise, both \(W\) and the distribution of
recorded travel budgets differ among conditions, so absolute detection
probabilities should not be compared causally across cell types.

| Condition | \(c\), \(\mu\)m | \(W=2L_{\max}/c\) | Median budget [IQR], \(c\) units | Physical target-area range, \(\mu\mathrm m^2\) |
|---|---:|---:|---:|---:|
| NK | 2.332 | 34.339 | 23.63 [12.84, 35.39] | 21.8–43.5 |
| Neutrophils 15 s | 3.138 | 85.061 | 25.58 [12.13, 44.79] | 39.4–78.8 |
| T HIV | 4.187 | 35.249 | 28.33 [17.83, 47.99] | 70.1–140.3 |

Detection was evaluated at relocation endpoints. The uniformly sampled initial
position was not itself counted as a detection. A relocation was accepted only
if its complete length fit within the remaining budget; the first overshooting
relocation was discarded rather than shortened.

### 4.2 Finite recorded-budget replay

Each track retained the sum of its recorded directed-run lengths as its travel
budget. Five constructions were compared:

| Class | Construction | Structure retained |
|---|---|---|
| Exact orientation | Original ordered directed-run vectors | Lengths, order, directions, turns, and laboratory orientation |
| Global rotation | One independent Haar-uniform 3-D rotation of the complete exact track per replay | All internal geometry, but not laboratory orientation |
| Per-vector rotation | Independently Haar-rotate every exact directed-run vector per replay | Exact ordered lengths, boundary runs, endpoint count, and budget; no angular correlations |
| Condition pool | Sample complete runs uniformly with replacement and orient every draw isotropically | Empirical condition-level length distribution; not focal-track lengths or endpoint count |
| Fitted model | Sample the fitted finite run-length model and an isotropic direction | Fitted marginal length law only |

The condition pool includes complete runs from the focal track and other tracks
in the same condition; it is not leave-one-track-out. There were 10,000
replays per track, class, target shape, and area. The selected finite analysis
therefore contains 186.30 million trajectory–target trials. Condition-level
probabilities average tracks equally, and paired 95% intervals were obtained
from 1,000 whole-track bootstrap resamples.

### 4.3 Unbounded first passage

An exact recorded path cannot be continued indefinitely without imposing an
arbitrary recycling rule. The unbounded experiment therefore compared only
the renewable condition-pool and fitted-model processes. Each continued until
first endpoint detection, with 2,000 independent trials per condition, class,
shape, and area (108,000 selected-condition trials). The outcome was cumulative
path length at detection in crossover units. Ratio intervals used the
independent Monte Carlo standard errors of the two processes.

## 5. Results

### 5.1 Fitted run-length distributions

The empirical distributions span broad but finite ranges. The fitted model
captures their central scale and tail trend, but formal adequacy does not hold
uniformly: only the NK fit was not rejected by the parametric-bootstrap test.
Search transfer was therefore assessed against an empirical resampling control
in addition to the distributional diagnostic. The selected-condition
complementary cumulative distributions, fit intervals, and segmentation
schematic are retained as a Supplementary diagnostic rather than occupying a
main-text figure.

### 5.2 Finite recorded-budget search

For a compact main-text summary, we collapsed the prespecified target grid
using an explicitly defined estimand. For condition \(k\), trajectory
construction \(m\), track \(j\), and target setting \(g\), let
\(\widehat p_{kjmg}\) be the detection fraction among the 10,000 replays. We
first average equally over the nine target shape–area settings within each
track and then average tracks equally:

\[
\overline p_{km}
=
\frac{1}{J_k}
\sum_{j=1}^{J_k}
\left(
\frac{1}{9}\sum_{g=1}^{9}\widehat p_{kjmg}
\right).
\]

The compact figure reports

\[
R_{km}=\frac{\overline p_{km}}{\overline p_{k,\mathrm{global}}},
\]

so each condition's globally rotated intact-track control equals one. Ratios
were recomputed in 10,000 paired whole-track bootstrap resamples. This is a
ratio of two grid-averaged probabilities, not an average of nine
target-specific ratios. The nine target settings define a fixed evaluation
grid and are not treated as biological replicates.

[![Compact finite-search comparison](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure1_finite_search_compact.png)](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure1_finite_search_compact.pdf)

**Figure 1 | Finite search under five trajectory constructions.**
Bars give the equal-track, target-grid-averaged detection probability relative
to the globally rotated intact-track control within each condition. Thus one
bar represents each trajectory experiment: exact orientation, global
rotation, independent rotation of every recorded run vector, condition-pool
resampling with isotropic directions, and the fitted finite run-length model
with isotropic directions. Error bars are 95% whole-track bootstrap intervals
for the normalized estimand; the dashed line marks 100% (ratio \(1\)) for the
global control. Target-specific probabilities and paired contrasts are
provided in the Supplementary diagnostic figure and source tables.

| Condition | Per-vector / global | Pool / per-vector | Fitted / pool |
|---|---:|---:|---:|
| NK | 1.435 (9 higher) | 1.049 (0) | 1.057 (0) |
| Neutrophils 15 s | 1.160 (3 higher) | 1.056 (0) | 0.914 (1 lower) |
| T HIV | 1.338 (9 higher) | 1.067 (7 higher) | 0.993 (0) |

Entries in this target-resolved companion table are descriptive median ratios
over the nine shape–area settings. Parentheses give the number of settings
whose unadjusted paired 95% interval excluded the null in the stated
direction. They are not replicate counts and are not the aggregation used for
the bar heights in Figure 1.

Global rotation minus exact orientation was compatible with zero in all 27
selected target settings. Laboratory-frame orientation therefore does not
explain the main difference between intact and randomized paths.

Independently rotating every vector of the same exact track raised the point
estimate in all 27 settings. Twenty-one increases were significant: all nine
for NK, three for neutrophils, and all nine for T cells. The median
per-vector/global probability ratios were \(1.435\), \(1.160\), and \(1.338\),
respectively. Because this intervention preserves exact ordered lengths,
boundary runs, endpoint count, and budget, it localizes the systematic effect
to directional and turning structure within the recorded paths.

Replacing the focal track's lengths with condition-pool draws produced a
smaller additional increase. The median pool/per-vector ratios were \(1.049\)
for NK, \(1.056\) for neutrophils, and \(1.067\) for T cells; only the seven
T-cell settings were significant. This contrast is composite because pool
sampling also changes the length multiset, removes boundary runs, and can
change the number of accepted endpoints.

The fitted model was close to the condition-pool control in the finite
experiment. Twenty-six of 27 paired intervals included zero. Median
fitted/pool probability ratios were \(1.057\) for NK, \(0.914\) for
neutrophils, and \(0.993\) for T cells; only neutrophil Line at \(A/c^2=8\)
showed a significant difference.

### 5.3 Unbounded first passage

The unbounded experiment contains only the two renewable processes. To produce
the corresponding compact summary, mean first-passage distance was averaged
equally over the same nine prespecified target settings:

\[
\overline D_{km}
=
\frac{1}{9}\sum_{g=1}^{9}\widehat D_{kmg},
\qquad
Q_{km}
=
\frac{\overline D_{km}}{\overline D_{k,\mathrm{pool}}}.
\]

The empirical condition-pool process is therefore one within each condition,
and values above one mean that the fitted process travels farther before
detection. Uncertainty is propagated from the independent Monte Carlo trial
estimates; the nine target settings again define the evaluation grid rather
than an inferential sample.

[![Compact unbounded-search comparison](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure2_unbounded_search_compact.png)](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure2_unbounded_search_compact.pdf)

**Figure 2 | Unbounded first-passage distance for renewable searches.**
Bars give target-grid-averaged mean detection distance relative to the
empirical condition-pool process within each condition. The empirical bar is
therefore 100% (ratio \(1\)), and the fitted-model bar tests transfer of the
fitted marginal run-length law. Error bars are 95% Monte Carlo ratio intervals;
the dashed line marks equality. Lower distance denotes more efficient
detection. The target-resolved shape–area results are retained as a
Supplementary diagnostic.

| Condition | Median fitted / empirical distance | Target settings differing from one |
|---|---:|---:|
| NK | 1.004 | 0 of 9 |
| Neutrophils 15 s | 1.039 | 4 of 9 longer |
| T HIV | 1.006 | 1 of 9 shorter |

The table reports descriptive medians of the nine target-specific ratios, not
the ratio of grid-averaged distances plotted in Figure 2. Target settings are
not biological replicates.

In unbounded first passage, 22 of 27 ratio intervals included one. The
condition-median fitted/empirical distance ratios were \(1.004\) for NK,
\(1.039\) for neutrophils, and \(1.006\) for T cells. All nine NK intervals
included one. Four neutrophil settings required significantly longer fitted
paths, while one T-cell setting required a significantly shorter fitted path.
Thus transfer was excellent for NK and T cells and modestly less accurate in
the largest neutrophil dataset.

Across the selected conditions, finite detection increased with target area.
At matched projected area, Disk targets were detected approximately 9% more
often than Ball targets and Line targets approximately 13% less often under
the randomized empirical control. These target-resolved patterns are shown in
the Supplementary diagnostic rather than the compact main-text summaries.
Equal projected area does not imply equal three-dimensional detection volume,
so these are geometric effects rather than evidence for shape-specific cell
behavior.

## 6. Interpretation

Three conclusions are supported.

First, the empirical search test is not merely a comparison of fitted
exponents. It uses each condition's measured length distribution, exact track
budgets, and exact higher-order trajectory geometry in a common
dimensionless search experiment.

Second, the new per-vector control changes the interpretation of the earlier
randomization result. Global rotation leaves search performance essentially
unchanged, while independently rotating each exact directed-run vector
increases detection without changing the track's lengths or endpoint count.
The dominant finite-search difference is therefore associated with
within-track angular structure—turning, persistence, and higher-order
directional correlations—rather than absolute orientation. Because all
angles are randomized together, the experiment does not identify one
particular turning statistic as causal.

Third, the fitted marginal run-length law generally reproduces the search
performance of independently oriented empirical runs, including in the two
conditions where formal goodness of fit rejects the distribution. Functional
agreement in a specified search task can therefore coexist with detectable
distributional misspecification. The model is useful as an effective surrogate
for the randomized process, but it does not reproduce intact biological paths
and does not establish that immune cells implement a Lévy mechanism.

The NK condition provides the cleanest positive example: its distributional
fit is not rejected, its directional-structure effect is strong, and fitted
and empirical randomized searches agree in both finite and unbounded
experiments. T cells provide a large cross-condition replication with similarly
good search transfer despite a rejected distributional fit. The neutrophil
condition is the largest-sample stress test and exposes a small but coherent
loss of unbounded efficiency under the fitted model. This combination avoids
presenting only best-case conditions.

## 7. Limitations

1. LTDB contains trajectories but no matched target coordinates or observed
   detection events. All encounter outcomes are counterfactual simulations.
2. The three conditions differ in cell identity, experimental context, and
   sampling interval. Their absolute probabilities and fitted parameters
   cannot be interpreted as causal cell-type differences.
3. The \(90^\circ\) threshold is literature informed but operational.
   Effective exponents vary under \(45^\circ\) and \(60^\circ\) segmentations,
   and the replay was not repeated for those thresholds.
4. The first and last runs are censored for fitting. Restoring their recorded
   segments in finite replay preserves observed budgets but does not recover
   the missing parts outside the imaging window.
5. Detection occurs only at relocation endpoints. Continuous swept-path
   detection could alter the effect of long runs.
6. The per-vector control removes all angular correlations simultaneously; it
   is not a single-turning-angle intervention.
7. Condition-pool sampling is not leave-one-track-out and changes more than
   length order alone.
8. Geometry is matched in normalized projected area, not common physical area
   or common three-dimensional capture volume.
9. Finite intervals condition on fitted \(c\), \(L_{\max}\), target geometry,
   and empirical pools. They do not propagate every source of experimental or
   between-video uncertainty. No multiplicity correction was applied to the
   27 descriptive target-cell contrasts.

## 8. Data, figures, and reproducibility

The selected-condition, publication-sized artifacts are stored under
[`publication_selected_conditions`](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions).
They are derived without altering the signed all-condition replay or the
additive per-vector-rotation stage. The two main-text figures are exported as
vector PDFs on an 8.7-cm one-column canvas, with PNG previews:

- [finite recorded-budget comparison](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure1_finite_search_compact.pdf);
- [unbounded first-passage comparison](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure2_unbounded_search_compact.pdf).

Useful source tables are:

- the [selected fit summary](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_fit_summary.csv);
- the [compact finite bar estimates](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_compact_finite_ratios.csv);
- the [compact unbounded bar estimates](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_compact_unbounded_ratios.csv);
- the [five-class finite summary](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_finite_summary.csv);
- the [paired finite contrasts](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_finite_contrasts.csv);
- the [unbounded summary](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_unbounded_summary.csv);
- the [unbounded contrasts](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_unbounded_contrasts.csv); and
- the [selected track inputs](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/selected_track_inputs.csv).

The larger diagnostics are deliberately not main-text figures. The
[directed-run and fitted-distribution panel](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure1_data_and_fit.pdf)
and the
[target-resolved finite and unbounded grid](../outputs/equal_run_weight_run_eligible_conditions/publication_selected_conditions/figures/figure2_search_performance.pdf)
are retained for Supplementary Information and auditability. The latter shows
both the absolute empirical and fitted unbounded detection distances and their
fitted-to-empirical ratios.

After the full replay and additive per-vector control have been generated, the
publication artifacts can be reproduced with:

```bash
PYTHONPATH=. .venv/bin/python -m ltdb_levy.cli \
  selected-publication-figures \
  --config configs/equal_run_weight_run_eligible_conditions.yaml
```

The complete simulations and all-condition figures are documented in the
[supplementary all-condition analysis](equal_run_all_conditions_analysis.md).

## References

1. Pizzagalli DU, Farsakoglu Y, Palomino-Segura M, *et al.* Leukocyte Tracking
   Database, a collection of immune cell tracks from intravital 2-photon
   microscopy videos. *Scientific Data*. 2018;5:180129.
   [doi:10.1038/sdata.2018.129](https://doi.org/10.1038/sdata.2018.129).
2. Georgantzoglou A, Poplimont H, Walker HA, Lämmermann T, Sarris M. A
   two-step search and run response to gradients shapes leukocyte navigation
   in vivo. *Journal of Cell Biology*. 2022;221(8):e202103207.
   [doi:10.1083/jcb.202103207](https://doi.org/10.1083/jcb.202103207).
3. Ganusov VV, Zenkov VS, Majumder B. Correlation between speed and turning
   naturally arises for sparsely sampled cell movements. *Physical Biology*.
   2023;20:016002.
   [doi:10.1088/1478-3975/acb18c](https://doi.org/10.1088/1478-3975/acb18c).
