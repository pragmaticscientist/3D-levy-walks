# NK-cell directed-run fit and search experiment

> **Weighting update.** This document records the earlier, NK-focused
> cell-balanced analysis. The requested equal-weight-per-run analysis across
> all run-eligible conditions is reported in
> [equal_run_all_conditions_analysis.md](equal_run_all_conditions_analysis.md).
> The older results below are retained as a weighting and eligibility
> sensitivity analysis, not as the new primary result.

## Scope and main conclusion

This document summarizes the analysis used to compare the theoretical
flat-to-power-law movement model with real Natural Killer (NK) cell
trajectories. It covers:

1. the LTDB data format and the derived directed-run data;
2. the fitted movement model;
3. why some runs and tracks were excluded, and which tracks were retained;
4. the statistical and simulation tests;
5. the results, plots, interpretation, and limitations.

The focused condition was:

> Natural Killer Cells, popliteal lymph node, Influenza Vaccine,
> sampling interval \(dt=30\) s.

The principal result is that the fitted model reproduces the search performance
of an independently resampled empirical run control reasonably well. However,
the original ordered NK trajectories detect targets less often than the
independent-step controls. Thus, the fitted run-length distribution captures
much of the one-step behavior, while the original tracks contain additional
sequence or correlation structure not represented by the independent model.

The fitted model is a useful surrogate for comparison with the theory, but the
fit should be described as an **effective truncated power-law model**, not as
proof that the cells execute an exact Lévy walk.

## 1. Data format

### 1.1 Raw LTDB files

The input files are the LTDB ground-truth trajectory files in
[`datasets/GT_TRACKS`](../datasets/GT_TRACKS). They are semicolon-delimited and
have a positional format:

```text
line 1: VideoID ; dx ; dy ; dz ; dt
line 2: ch0 ; ch1 ; ch2 ; ch3 ; ch4
line 3+: TrackID ; x ; y ; z ; t
```

- Line 1 contains the video identifier, voxel calibration, and frame interval.
- Line 2 contains channel-visibility flags and is not an observation.
- Lines 3 onward contain a track identifier, three-dimensional position, and
  frame index.
- The \(x,y,z\) coordinates are already in micrometres. They were not
  multiplied by \(dx,dy,dz\), which would scale them twice.
- The variable \(t\) is a frame index. Physical time is \(t\,dt\).

The parser assigns globally unambiguous identifiers such as
`CS018_a#1`, where the part before `#` identifies the video and the part after
it is the file-local track number.

Biological metadata were joined by video identifier. Conditions were kept
separate by

\[
(\text{cell type},\ \text{organ},\ \text{stimulus},\ dt).
\]

Including \(dt\) is important because changing the sampling interval changes
the observed displacement lengths and turning angles.

### 1.2 Corpus-level audit

The complete audited corpus contained:

| Stage | Count |
|---|---:|
| Candidate files | 45 |
| Selected biological files | 44 |
| Parse failures | 0 |
| Raw observations | 44,722 |
| Raw tracks | 728 |
| Retained observations | 44,657 |
| Retained tracks | 710 |
| Frame-to-frame displacements | 43,947 |
| Primary directed runs | 17,981 |
| Runs included in condition pools | 16,578 |
| Conditions | 16 |
| Conditions meeting the fitting threshold | 8 |

The synthetic file `SQUARE_GT..csv` was excluded by configuration. Eighteen
tracks with fewer than six observations were excluded, accounting for 65
observations. There were no duplicate rows, duplicate frames, or missing-frame
fragments in the selected corpus.

All exclusions are recorded rather than silently deleted. The relevant audit
files are:

- [file selection audit](../outputs/primary/inspect/file_selection_audit.csv);
- [preprocessing exclusions](../outputs/primary/preprocess/preprocess_exclusions.csv);
- [segmentation audit](../outputs/primary/preprocess/segment_audit.csv);
- [condition-pool audit](../outputs/primary/preprocess/pool_audit.csv).

### 1.3 Focused NK condition

The 30-s NK condition contains three videos:

| Video | Observations | Raw tracks | Fit-contributing tracks | Complete fit runs |
|---|---:|---:|---:|---:|
| `CS018_a` | 633 | 15 | 15 | 306 |
| `LTDB019_a` | 509 | 12 | 11 | 209 |
| `LTDB020_a` | 35 | 2 | 1 | 2 |
| **Total** | **1,177** | **29** | **27** | **517** |

The NK video `LTDB018_a` was not pooled with these data because it was sampled
at 42 s rather than 30 s. Its 31 complete runs from four contributing tracks
also failed the prespecified minimum of 15 tracks and 100 runs.

### 1.4 From raw positions to directed-run vectors

The model was **not fitted to raw positions or directly to individual frame
steps**. It was fitted to directed-run end-to-end lengths.

For consecutive positions, the pipeline first calculated:

- the displacement vector
  \(\mathbf u_i=(\Delta x_i,\Delta y_i,\Delta z_i)\);
- its Euclidean length;
- elapsed time and speed;
- the turning angle

\[
\theta_i=\operatorname{atan2}
\left(
\lVert\mathbf u_{i-1}\times\mathbf u_i\rVert,
\mathbf u_{i-1}\cdot\mathbf u_i
\right).
\]

The primary segmentation merged consecutive moving displacements until:

- the turning angle was strictly greater than \(90^\circ\);
- a zero-length displacement was encountered; or
- the track fragment ended.

No positive speed cutoff was imposed in the primary analysis. The displacement
that exceeded \(90^\circ\) started the new run. For every directed run we saved:

- its net three-dimensional vector;
- its end-to-end length, used for fitting;
- its cumulative within-run path length;
- duration, mean speed, straightness, and boundary-censoring flags.

The complete derived run table is [runs.csv](../outputs/primary/preprocess/runs.csv).

## 2. Model fitting

### 2.1 Fitted distribution

For directed-run length \(\ell\), we fitted

\[
f(\ell\mid\mu,c,L_{\max})
=
a(\mu,c,L_{\max})
\min\left[
1,\left(\frac{\ell}{c}\right)^{-\mu}
\right],
\qquad 0\leq\ell\leq L_{\max},
\]

where:

- \(c\) is the physical crossover between the flat core and power-law tail;
- \(\mu\) is the tail exponent;
- \(L_{\max}\) is the finite upper cutoff;
- \(a\) is the normalization constant.

\(L_{\max}\) was set to the empirical maximum within the condition. The
crossover and exponent were estimated jointly by maximum likelihood. The
configured value \(3.75\,\mu\text{m}\), half a reference \(7.5\,\mu\text{m}\)
cell diameter, was only the optimizer starting point.

The fit used cell-balanced weights. If track \(j\) supplied \(n_j\) complete
runs and there were \(J\) tracks, each run from that track received weight

\[
w_{ij}=\frac{1}{J n_j}.
\]

Every track therefore contributed equal total weight, preventing long tracks
from dominating the fit.

### 2.2 Primary fitted values

| Quantity | Estimate | Track-bootstrap 95% interval |
|---|---:|---:|
| Number of tracks | 27 | — |
| Number of complete runs | 517 | — |
| Effective sample size | 227.6 | — |
| \(\mu\) | 2.3253 | 1.8199–3.1505 |
| \(c\) | \(2.2637\,\mu\text{m}\) | 1.5264–3.2346 \(\mu\text{m}\) |
| \(L_{\max}\) | \(40.0429\,\mu\text{m}\) | fixed at empirical maximum |
| KS statistic | 0.04478 | — |
| Parametric-bootstrap KS \(p\) | 0.3593 | 500 bootstrap replicates |

The primary goodness-of-fit test did not reject the fitted mixture. However,
the truncated lognormal ranked first in track-held-out predictive comparison,
with the mixture second:

| Model | Mean held-out log density | Rank |
|---|---:|---:|
| Truncated lognormal | -2.0528 | 1 |
| Fitted C-mixture | -2.1144 | 2 |

The mixture was therefore not uniquely supported by the data. We report
\(\mu\) as an effective exponent suitable for a theoretical comparison.

### 2.3 Why the crossover was estimated

Fixing \(c\) at different arbitrary values changed the fitted exponent
substantially:

| Fixed \(c\) (\(\mu\)m) | Fitted \(\mu\) |
|---:|---:|
| 0.9375 | 1.509 |
| 1.875 | 2.081 |
| 3.750 | 3.343 |
| 7.500 | 6.316 |
| 15.000 | 10.000, at the upper fit bound |

This strong coupling is why \(c\) was estimated from the data rather than
chosen solely from a reference cell size.

### 2.4 Normalization used in the simulations

After fitting, all lengths were divided by the fitted crossover:

\[
\ell^*=\frac{\ell}{c}.
\]

Consequently:

- one model unit is \(2.26373585\,\mu\text{m}\);
- \(c^*=1\);
- the detection radius is \(d=1\), or \(2.26373585\,\mu\text{m}\);
- \(L_{\max}^*=17.6889\);
- the torus width is
  \(W=2L_{\max}^*=35.3777\), or \(80.0859\,\mu\text{m}\).

The effective projected capture areas \(4,8,16,32\) model units\(^2\)
correspond to \(20.498,40.996,81.992,163.984\,\mu\text{m}^2\).

## 3. Why runs and tracks were excluded

### 3.1 General track-quality exclusions

Tracks with fewer than six observations were excluded because they cannot
support a stable displacement and turning-angle sequence. This removed 18 of
728 tracks across the full corpus. None of the 29 tracks in the focused 30-s
NK condition failed this rule.

Tracks would have been split at missing frames to avoid creating artificial
long jumps, and duplicate track-frame rows would have triggered an error.
Neither issue occurred in the selected data.

### 3.2 Boundary-censored runs

The first and last directed runs of every recorded track were excluded from the
distribution fit. A cell normally enters the field of view after its first
true run has already started and leaves before its last true run has ended.
Their observed lengths are therefore right- or left-censored and tend to
underestimate the complete biological run length.

In the focused NK condition:

| Stage | Runs | Tracks represented |
|---|---:|---:|
| All segmented runs | 574 | 29 |
| Boundary runs excluded from fitting | 57 | — |
| Complete runs retained for fitting | 517 | 27 |

Two tracks had no complete interior run:

| Excluded track | Observations | Segmented runs | Reason |
|---|---:|---:|---|
| `LTDB019_a#12` | 30 | 2 | one first and one last run; both censored |
| `LTDB020_a#2` | 19 | 1 | the only run was both first and last |

These tracks were not corrupt or too short at the raw-data level. They were
excluded because they supplied no complete directed run from which to estimate
the run-length distribution. They were also omitted from the focused replay so
that every replayed real trajectory belonged to the same 27-track population
that contributed to the fit.

### 3.3 Tracks retained

The 27 retained track identifiers were:

```text
CS018_a#1     CS018_a#2     CS018_a#3     CS018_a#4
CS018_a#5     CS018_a#6     CS018_a#7     CS018_a#8
CS018_a#9     CS018_a#10    CS018_a#11    CS018_a#12
CS018_a#13    CS018_a#14    CS018_a#15

LTDB019_a#1   LTDB019_a#2   LTDB019_a#3   LTDB019_a#4
LTDB019_a#5   LTDB019_a#6   LTDB019_a#7   LTDB019_a#8
LTDB019_a#9   LTDB019_a#10  LTDB019_a#11

LTDB020_a#1
```

For fitting and empirical resampling, only their 517 complete runs were used.
For the exact finite-trajectory replay, their observed first and last vectors
were restored because the goal there was to reproduce the complete recorded
path. This produced 571 ordered, positive-length directed-run vectors across
the 27 tracks. Thus:

- boundary runs were excluded when estimating the distribution because their
  intended lengths are censored;
- the same observed vectors were retained when replaying the finite recorded
  path because they are part of the observed trajectory and its distance
  budget.

No raw data were deleted. The inclusion flag and exclusion reason for every run
remain available in the derived tables.

## 4. Tests that were performed

### 4.1 Fit validation

1. **Whole-track clustered bootstrap.** Whole tracks were resampled 1,000 times
   to obtain uncertainty intervals for both \(\mu\) and \(c\). This preserves
   dependence among runs from the same cell.
2. **Parametric-bootstrap KS test.** Five hundred samples were generated from
   the fitted model. Both \(\mu\) and \(c\) were refitted in every replicate
   before recomputing the KS statistic.
3. **Track-held-out model comparison.** Ten-fold cross-validation compared the
   fitted mixture with truncated lognormal, exponential, Weibull, and
   power-law-with-cutoff alternatives. Entire tracks were held out.
4. **Fixed-crossover sensitivity.** The model was refitted at five prescribed
   crossover scales to measure the dependence of \(\mu\) on \(c\).
5. **Segmentation sensitivity.** The fit was repeated for turning thresholds
   \(45^\circ,60^\circ,90^\circ\) and speed cutoffs of none, 1, 2, and
   \(4\,\mu\text{m}\,\text{min}^{-1}\).

The segmentation sensitivity used 100 KS bootstrap replicates per setting, so
its \(p\)-values have coarser resolution than the 500-replicate primary fit:

| Angle | Speed cutoff (\(\mu\)m/min) | Tracks | Runs | \(c\) (\(\mu\)m) | \(\mu\) | \(L_{\max}\) (\(\mu\)m) | KS \(p\) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 45° | none | 29 | 862 | 1.487 | 2.072 | 43.727 | 0.0495 |
| 45° | 1 | 29 | 716 | 1.755 | 2.173 | 43.727 | 0.0099 |
| 45° | 2 | 28 | 497 | 2.788 | 2.694 | 43.727 | 0.0099 |
| 45° | 4 | 24 | 241 | 3.377 | 2.649 | 43.727 | 0.0099 |
| 60° | none | 28 | 757 | 1.526 | 2.058 | 43.727 | 0.0495 |
| 60° | 1 | 28 | 634 | 1.683 | 2.084 | 43.727 | 0.0099 |
| 60° | 2 | 28 | 447 | 2.532 | 2.280 | 43.727 | 0.0099 |
| 60° | 4 | 24 | 217 | 3.283 | 2.337 | 43.727 | 0.0099 |
| 90° | none | 27 | 517 | 2.264 | 2.325 | 40.043 | 0.3663 |
| 90° | 1 | 27 | 466 | 2.394 | 2.421 | 40.043 | 0.0099 |
| 90° | 2 | 27 | 369 | 2.846 | 2.482 | 26.866 | 0.0099 |
| 90° | 4 | 23 | 193 | 3.377 | 2.392 | 26.460 | 0.0099 |

The primary no-speed, \(90^\circ\) setting was the only clearly non-rejected
setting. This demonstrates that model adequacy depends on how directed runs and
pauses are defined; it must be reported as a sensitivity, not hidden.

### 4.2 Finite-trajectory search experiment

Every retained track supplied its observed cumulative directed-run distance as
the path budget. The mean budget was 25.633 model units, or
\(58.025\,\mu\text{m}\). A uniformly random starting position in the torus was
used, and an overshooting final synthetic step was discarded.

Four trajectory classes were compared:

| Class | Construction | What it tests |
|---|---|---|
| Exact orientation | Original ordered directed-run vectors | Search behavior of the recorded path |
| Exact + global rotation | One random 3-D rotation of the complete path | Effect of absolute orientation while preserving all run lengths and turning angles |
| Uniform empirical vector + rotation | Choose one of the 517 vectors uniformly with replacement and independently rotate every draw | Pooled empirical run lengths with order and directional correlations removed |
| Fitted \(\mu\) | Draw a fitted run length and an isotropic direction independently | The theoretical fitted model |

There were 10,000 replays per track, trajectory class, target shape, and target
area. This gave 270,000 trials per plotted cell and 12.96 million finite
trajectory/target evaluations in total. Detection probabilities were averaged
equally over tracks, and uncertainty intervals used 1,000 whole-track bootstrap
replicates.

The uniform-vector experiment gives every one of the 517 runs equal sampling
probability, as requested. This is a step-weighted empirical control, whereas
the parametric fit used cell-balanced weights. That difference is retained
deliberately but should be remembered when comparing the orange and green
curves.

### 4.3 Unbounded first-passage experiment

Exact observed paths terminate, so a complete first-passage time cannot be
measured from them without inventing a recycling rule. The unbounded experiment
therefore compared only the two renewable classes:

- uniformly sampled empirical run vector plus independent rotation;
- fitted \(\mu\) length plus isotropic direction.

Each path continued until its first target detection. We ran 2,000 independent
trials for each trajectory class, target shape, and target area, for 48,000
trials in total. The recorded outcomes were:

- cumulative path length at detection;
- number of directed-run steps to detection.

### 4.4 Geometry and implementation checks

- Target shapes were Ball, Disk, and Line.
- Shapes were matched by maximum face-on projected area including \(d=1\).
- Detection was checked at directed-run endpoints, not continuously along the
  swept segment.
- The initial position was not counted as a detection.
- All 12 shape/area detection neighbourhoods fit inside the torus without
  periodic self-overlap.
- Even the longest line target retained 9.474 model units of clearance to the
  torus half-width.
- Random rotations were Haar-uniform in three dimensions and preserved vector
  lengths.
- The full software test suite passed: **59 tests passed**.

The focused experiment used the Python process-pool implementation on 32 cores.
Its target geometry matches the existing C implementation and the normalization
matches the C convention \(c=1,d=1\). The existing C sampler was not used for
this replay because it truncates \(L_{\max}\) and the domain side to integers
and cannot consume empirical vector pools or finite recorded paths.

## 5. Results and interpretation

### 5.1 Finite detection probabilities

![Finite-budget detection probability](../outputs/primary/focused_replay/figures/finite_detection_probability.png)

**Figure 1.** Probability of detecting the target before exhausting the
observed track-specific distance budget. Points show the track-averaged
probability; bars are whole-track bootstrap intervals.

For the largest projected area, the estimated detection probabilities were:

| Shape | Exact orientation | Global rotation | Empirical independent | Fitted \(\mu\) |
|---|---:|---:|---:|---:|
| Ball | 0.940% | 0.914% | 1.209% | 1.271% |
| Disk | 0.577% | 0.570% | 0.802% | 0.790% |
| Line | 0.599% | 0.572% | 0.813% | 0.772% |

The finite experiment supports three conclusions.

#### Absolute orientation contributes little

The exact and globally rotated paths were very similar. Most of the 12
whole-track bootstrap intervals for their difference included zero. Therefore,
the absolute orientation of the recorded trajectory relative to the laboratory
coordinate system is not the principal driver of detection.

#### Independent empirical steps detect more often than ordered tracks

The independently sampled empirical control had a higher detection probability
than the globally rotated exact path in all 12 shape/area cells. Every
whole-track bootstrap interval for this contrast excluded zero.

This shows that the pooled one-step distribution alone does not reproduce the
search behavior of the complete recorded tracks. The original paths contain
structure that reduces detection over the same distance budget.

This contrast does **not** isolate turning angles perfectly. It also replaces
each track's run-length sequence by draws from the pooled step-weighted
distribution and can change the number of endpoints generated within a fixed
distance budget. A cleaner turning-only control would preserve the exact
sequence of observed run lengths while independently randomizing each run
direction.

#### The fitted model is close to the independent empirical control

For the fitted-model minus empirical-control contrast, 7 of 12 bootstrap
intervals included zero. The five intervals excluding zero had relative
differences of approximately 5–13%, and the fitted model was not systematically
above or below the empirical control.

Thus, the fitted model is a reasonable approximation to the independent
empirical-run experiment, but it is not an exact replacement.

### 5.2 Unbounded detection distance

![Unbounded detection path length](../outputs/primary/focused_replay/figures/unbounded_detection_distance.png)

**Figure 2.** Mean cumulative path length required for detection when the
renewable empirical and fitted trajectories are allowed to continue until
detection. The vertical axis is logarithmic.

At projected area 32:

| Shape | Empirical mean distance | Fitted mean distance | Fitted/empirical ratio (95% MC interval) |
|---|---:|---:|---:|
| Ball | 3424 model units | 2781 model units | 0.812 (0.761–0.866) |
| Disk | 4727 model units | 4104 model units | 0.868 (0.816–0.924) |
| Line | 4183 model units | 4036 model units | 0.965 (0.906–1.027) |

Across the complete grid, 10 of 12 distance-ratio intervals included one. The
exceptions were the largest Ball and Disk targets, for which the fitted model
detected the target in approximately 19% and 13% less cumulative distance,
respectively. The largest Line result remained compatible with equal mean
distance.

The fitted model has a larger mean run length:

| Run generator | Mean length (model units) | Mean length (\(\mu\)m) |
|---|---:|---:|
| Uniform empirical vector | 1.1048 | 2.5010 |
| Fitted model | 1.3619 | 3.0830 |

It consequently uses fewer run endpoints before detection. The recorded
first-hit step is a count of directed runs, not physical time in seconds.
Converting it to time would require retaining or modelling run durations.

### 5.3 Effect of target area and shape

Detection probability increased, and unbounded detection distance decreased,
as projected capture area increased. This is the expected direction for a
larger target.

Shapes with the same projected area do not have the same three-dimensional
capture volume. For example, at area 32 the Ball capture neighbourhood has a
larger volume than the Disk or Line, explaining its higher detection
probability. Shape differences in the plots therefore reflect geometry as well
as movement.

### 5.4 Overall interpretation

The experiments separate three levels of description:

1. **Absolute orientation:** contributes little.
2. **Ordered track structure:** matters; exact tracks detect less often than
   independent empirical steps.
3. **Marginal run-length model:** the fitted distribution reproduces the
   independent empirical control in most, but not all, target cells.

A concise paper-level interpretation is:

> Using directed runs extracted from 27 NK-cell tracks, we fitted a truncated
> flat-to-power-law distribution with effective exponent
> \(\mu=2.33\), crossover \(c=2.26\,\mu\text{m}\), and cutoff
> \(L_{\max}=40.04\,\mu\text{m}\). After normalizing by \(c\), the fitted model
> reproduced the detection performance of independently resampled empirical
> runs across most target geometries and sizes. Random global rotation of the
> complete observed trajectories had little effect, whereas replacing the
> ordered tracks by independent empirical runs increased detection
> consistently. The fitted run-length law therefore captures much of the
> one-step search behavior, while higher-order trajectory structure outside
> the independent model reduces the efficiency of the observed paths.

## 6. Limitations

1. LTDB contains cell trajectories but no biological targets or observed
   detection events. All target-search outcomes are counterfactual simulations.
2. Directed-run construction depends on turning-angle, speed, and sampling
   choices. The sensitivity grid shows that model adequacy is not invariant to
   these choices.
3. The goodness-of-fit result does not establish an exact Lévy law. A truncated
   lognormal had better held-out prediction for this condition.
4. The finite empirical-independent control does not isolate turning angles
   alone because it also replaces run lengths and changes endpoint counts.
5. Uniform empirical-vector sampling is step-weighted, while the primary fit is
   cell-balanced.
6. Detection is endpoint-only. A continuous swept-path encounter rule could
   produce different results, especially for long steps.
7. Equal projected area does not imply equal three-dimensional target volume.
8. Unbounded Monte Carlo intervals quantify simulation error conditional on
   the fitted parameters and empirical pool. They do not propagate uncertainty
   in \(\mu,c,L_{\max}\) or biological sampling.
9. A directed-run step is not a fixed physical time interval. First-hit step
   counts should not be labelled as seconds without a duration model.
10. Cell, video, and tissue heterogeneity may not be fully represented by 27
    fit-contributing tracks from three videos.

## 7. Reproducibility and artifacts

The analysis configuration is [configs/primary.yaml](../configs/primary.yaml).
The focused replay was run with:

```bash
.venv/bin/python -m ltdb_levy focused-replay \
  --config configs/primary.yaml \
  --projected-areas 4 8 16 32 \
  --finite-replays 10000 \
  --unbounded-trials 2000 \
  --workers 32 \
  --finite-batch-size 1000 \
  --unbounded-chunk-trials 25 \
  --unbounded-step-block 256
```

Key machine-readable outputs are:

- [primary fit](../outputs/primary/fit/mixture_fits.csv);
- [segmentation fit sensitivity](../outputs/primary/fit/segmentation_fit_sensitivity.csv);
- [model comparison](../outputs/primary/fit/model_comparison.csv);
- [focused replay settings](../outputs/primary/focused_replay/settings.json);
- [retained replay tracks](../outputs/primary/focused_replay/track_inputs.csv);
- [target geometry checks](../outputs/primary/focused_replay/targets.csv);
- [finite summary](../outputs/primary/focused_replay/finite_summary.csv);
- [finite paired contrasts](../outputs/primary/focused_replay/finite_contrasts.csv);
- [unbounded summary](../outputs/primary/focused_replay/unbounded_summary.csv);
- [unbounded contrasts](../outputs/primary/focused_replay/unbounded_contrasts.csv).
