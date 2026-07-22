# Finite endpoint detection by empirical and fitted immune-cell trajectories

## Summary

We compared a finite run-length search model with three-dimensional immune-cell
trajectories from the Leukocyte Tracking Database (LTDB). The analysis used
natural killer (NK) cells and neutrophils imaged in popliteal lymph nodes after
influenza vaccination, and T cells imaged in popliteal lymph nodes in an HIV
humanized-mouse model. Recorded tracks were segmented into directed runs at
three-dimensional turns greater than \(90^\circ\). Complete interior runs were
then fitted, with equal weight per run, to a bounded flat-core/power-law-tail
density. The crossover was not inferred from the data in this sensitivity
analysis: it was fixed to the displacement corresponding to the conventional
\(2~\mu\mathrm m\,\mathrm{min}^{-1}\) leukocyte arrest threshold over one
sampling interval.

The fitted model was tested as a search surrogate by placing Ball, Disk, and
Line targets in a periodic cube and imposing the recorded path-length budget of
each source track. The principal comparison is between (i) the exact ordered
run lengths of a track with every run vector independently rotated in three
dimensions and (ii) independent isotropic relocations drawn from the fitted
run-length model. Target dimensions were selected so that all three shapes had
the same normalized face-on projected area, denoted by \(\Delta\), of their full
detection neighborhood. Detection was intermittent: only relocation endpoints
were tested.

Across the 72 displayed condition-by-shape-by-\(\Delta\) settings, the fitted
process had a lower point estimate of detection probability than the
independently rotated empirical process. This is a useful transfer test, but it
is not evidence that the measured run lengths follow the proposed distribution:
the conditional IID-run parametric-bootstrap KS check rejected the
fixed-crossover model in all three conditions. The fitted exponents are
therefore reported as effective fixed-scale sensitivity parameters.

## 1. Dataset and biological conditions

### 1.1 Source and file format

We used the manually curated consensus ground truth in the LTDB, a collection
of leukocyte-centroid tracks obtained by intravital two-photon microscopy in
living mice [1]. The corresponding data collection is deposited on Figshare
[2]. LTDB reports \(x,y,z\) coordinates directly in micrometers; the integer
frame index was converted to time using the file-specific acquisition interval
\(dt\). Each semicolon-delimited ground-truth file begins with

```text
VideoID ; dx ; dy ; dz ; dt
```

followed by channel-visibility flags and records of the form

```text
TrackID ; x ; y ; z ; frame
```

The complete analyzed LTDB corpus contained 44 channel-specific trajectory
files from 38 acquisition videos, 44,722 raw annotations, and 728 raw tracks.
After quality control it contained 44,657 annotations from 710 tracks. The
three conditions used here comprised 20 trajectory files from 15 acquisitions,
28,240 raw annotations, and 437 tracks. After quality control, 28,218
annotations from 431 tracks remained; 414 tracks supplied at least one complete
interior run and entered the fit and replay.

We used the following reproducible subset rule for the compact display: retain
the NK example used in the preceding analysis and, among the remaining
conditions with at least 100 complete runs, retain the two conditions with the
largest fitting samples.

| Condition | Raw files / acquisitions / tracks | Post-QC tracks | Replay files / acquisitions / tracks | Complete fitting runs | Replay run vectors |
|---|---:|---:|---:|---:|---:|
| NK, popliteal LN, influenza vaccine, 30 s | 3 / 3 / 29 | 29 | 3 / 3 / 27 | 517 | 571 |
| Neutrophils, popliteal LN, influenza vaccine, 15 s | 10 / 6 / 234 | 229 | 10 / 6 / 224 | 5,829 | 6,277 |
| T cells, popliteal LN, HIV model, 15 s | 7 / 6 / 174 | 173 | 6 / 5 / 163 | 3,899 | 4,225 |
| **Total** | **20 / 15 / 437** | **431** | **19 / 14 / 414** | **10,245** | **11,073** |

The raw trajectory-file IDs were `CS018_a`, `LTDB019_a`, and `LTDB020_a` for NK
cells; `CS014_a`, `LTDB004_a/b`, `LTDB005_a/b`, `LTDB006_a`, `LTDB007_a/b`,
and `LTDB017_a/b` for neutrophils; and `CS001_a`, `CS002_a`, `CS003_a`,
`CS006_a`, `LTDB012_a/b`, and `LTDB013_a` for T cells. `CS001_a` supplied no
complete run and therefore did not enter fitting or replay. LTDB distinguishes
standard `LTDB` videos, in which all expected visible cells were tracked, from
`CS` challenge videos containing selected or difficult cells [1]. The combined
sample should consequently not be interpreted as a random census of cells in
the imaged tissue.

A later LTDB methods paper reports that `LTDB018` and `LTDB020` were excluded
from its validation because errors were detected in the supplied ground truth
[11]. The current NK result includes one replay track from `LTDB020_a`,
contributing two of 517 complete fitting runs and four of 571 replay vectors.
As a diagnostic deletion, removing those two fitting runs changes
\(\hat\mu\) from 1.673053 to 1.675458 while leaving \(L_{\max}\) unchanged.
Re-averaging the already simulated outcomes over the other 26 tracks changes
any displayed NK probability by at most \(2.57\times10^{-4}\) in absolute
terms (3.85% relative), and independent rotations retain the larger point
estimate in all 18 NK target settings. This is not a complete exclusion rerun,
because the fitted-model paths were not regenerated at the updated exponent.
A submission analysis should either rerun after excluding `LTDB020_a` or
retain it with this limitation stated explicitly.

The influenza-vaccine context of NK-cell and neutrophil motility in the
popliteal lymph node is described in Refs. 3 and 4. The T-cell recordings derive
from the humanized-mouse HIV setting described by Murooka *et al.* [5]. The
authors reported that infected lymph-node cells had an average two-dimensional
velocity of approximately \(7~\mu\mathrm m\,\mathrm{min}^{-1}\) [5]. That biological motility
rate is distinct from the \(2~\mu\mathrm m\,\mathrm{min}^{-1}\) *arrest
threshold* used below.

The pooled T-cell condition requires one qualification. Of the 163 analyzed
tracks, 159 were HIV-GFP tracks and four were CMTMR-labeled comparator tracks
from `LTDB012_a`. These comparator tracks contributed 121 of the 3,899 complete
fitting runs (3.1%) and 129 of the 4,225 replay vectors. They were pooled by the
condition key used in this analysis—cell type, organ, stimulus, and sampling
interval—so
this panel should be described as the “T-cell HIV-model condition,” not as a
sample consisting exclusively of infected T cells. One of the seven raw T-cell
videos supplied no complete run, explaining the six replay videos in the table.

The selected conditions are illustrative empirical examples rather than a
controlled three-cell-type comparison. Their sampling intervals, experimental
contexts, track-length distributions, and numbers of recorded cells differ.
Moreover, LTDB contains trajectories but no paired target locations or observed
encounters. All detections described below are therefore counterfactual replay
outcomes.

## 2. Preprocessing and directed-run construction

Observations were stably sorted by track and frame. Duplicate track-frame
records were treated as errors; a gap of more than one frame split a track into
separate fragments. Fragments containing fewer than six positions were
discarded. We applied no interpolation, spatial smoothing, or positive speed
filter.

For successive positions \(\mathbf x_i\), the frame displacement was

\[
\mathbf u_i=\mathbf x_{i+1}-\mathbf x_i .
\]

For consecutive nonzero displacements, the three-dimensional turning angle was
computed as

\[
\theta_i=
\operatorname{atan2}\!\left(
\lVert\mathbf u_{i-1}\times\mathbf u_i\rVert,
\mathbf u_{i-1}\cdot\mathbf u_i
\right).
\]

The `atan2` form remains stable near both \(0\) and \(\pi\). Consecutive
displacements were merged while \(\theta_i\leq90^\circ\); a displacement for
which \(\theta_i>90^\circ\) began the next directed run. An angle exactly equal
to \(90^\circ\) therefore did not split a run. A zero-length displacement
terminated the current run and supplied no direction.

Georgantzoglou *et al.* classified \(90^\circ\)–\(180^\circ\) reorientations
of neutrophils in vivo as large U-turns, providing a biological anchor for the
\(90^\circ\) boundary [6]. That paper did not propose the present run-merging
algorithm and did not validate the threshold for NK or T cells. We therefore
treat \(90^\circ\) as an operational, literature-informed definition, not a
universal biological constant. This distinction matters because sampling
frequency can itself alter measured speed–turning relationships [7]. Conditions
with different \(dt\) were consequently not pooled.

The vector assigned to a run was the vector sum of its frame displacements,

\[
\mathbf r_k=\sum_{i\in k}\mathbf u_i,
\qquad
\ell_k=\lVert\mathbf r_k\rVert .
\]

Thus \(\ell_k\) is the end-to-end length of the directed-run vector, not the
sum of framewise path lengths within that run. The first and last runs of each
continuous track fragment were excluded from fitting because the observation
window or a frame gap can censor their full duration. They were restored for
replay because they are part of the recorded trajectory and its finite travel
budget. Every selected replay track contained one fragment, so in this dataset
the distinction between track and fragment does not change the counts. This
accounts for the larger number of replay vectors than fitting runs in Table 1.

## 3. Fixed-arrest-scale fitting

### 3.1 Why \(c_{\mathrm{arrest}}\) was fixed

The sensitivity analysis fixed the model crossover to the distance traveled
at the conventional arrest-coefficient cutoff of
\(2~\mu\mathrm m\,\mathrm{min}^{-1}\) during one frame [4,8,9]. Writing
\(dt_{\min}=dt/(60~\mathrm{s\,min}^{-1})\) for the sampling interval expressed
in minutes,

\[
c_{\mathrm{arrest}}
=
\left(2~\mu\mathrm m\,\mathrm{min}^{-1}\right)dt_{\min}.
\]

This gives \(c_{\mathrm{arrest}}=1.0~\mu\mathrm m\) at \(dt=30\) s for NK
cells and \(c_{\mathrm{arrest}}=0.5~\mu\mathrm m\) at \(dt=15\) s for
neutrophils and T cells. The cited literature uses
\(2~\mu\mathrm m\,\mathrm{min}^{-1}\) to classify arrest; it does **not**
identify this one-frame displacement as the crossover of our length density.
Accordingly, \(c_{\mathrm{arrest}}\) is a fixed operational sensitivity scale. It
is not a measured mean cell speed, cell diameter, localization-noise estimate,
or threshold used to remove slow displacements during preprocessing.

### 3.2 Probability model and likelihood

For complete run lengths \(\ell\), we fitted

\[
f(\ell\mid\mu,c,L_{\max})
=
\frac{1}{Z(\mu,c,L_{\max})}
\min\!\left\{1,
\left(\frac{\ell}{c}\right)^{-\mu}\right\},
\qquad 0\leq\ell\leq L_{\max},
\]

where \(Z\) is the exact normalizing integral. The density has a uniform core
on \([0,c]\), a power-law-shaped tail on \([c,L_{\max}]\), and finite support.
It is therefore not a pure, unbounded power law.

Within each condition, \(c\) was fixed to \(c_{\mathrm{arrest}}\), and
\(L_{\max}\) was set to the empirical maximum complete-run length. We estimated
only \(\mu\in[1.001,10]\) by maximizing the exactly normalized likelihood.
Every complete run received equal likelihood weight. Tracks with more complete
runs consequently contributed more information to the point estimate, whereas
tracks—not individual runs—were the resampling units for uncertainty.

This use of maximum likelihood and a simulation-calibrated KS diagnostic
follows the general fitting principles advocated by Clauset, Shalizi, and
Newman [10], but the bounded flat-core model and its likelihood are specific to
the present theory. A 95% interval for \(\mu\) was obtained from 1,000
whole-track percentile-bootstrap samples. In each track-bootstrap replicate,
all runs again had equal weight and \(L_{\max}\) was recomputed as the empirical
maximum. Model adequacy was checked with 500 conditional IID-run
parametric-bootstrap Kolmogorov–Smirnov replicates. Each synthetic replicate
had the observed run count; \(\mu\) was refitted with the observed \(c\) and
data-derived \(L_{\max}\) held fixed. This calibration does not preserve
within-track dependence and is therefore a conditional run-level diagnostic,
not an unconditional clustered goodness-of-fit test.

| Condition | \(c_{\mathrm{arrest}}\), \(\mu\mathrm m\) | \(\hat\mu\) [track-bootstrap 95% CI] | \(L_{\max}\), \(\mu\mathrm m\) | \(L'_{\max}=L_{\max}/c\) | KS statistic | Bootstrap \(p\) |
|---|---:|---:|---:|---:|---:|---:|
| NK | 1.0 | 1.673 [1.456, 1.806] | 40.043 | 40.043 | 0.138 | 0.001996 |
| Neutrophils | 0.5 | 1.274 [1.134, 1.300] | 133.446 | 266.893 | 0.240 | 0.001996 |
| T cells | 0.5 | 1.072 [1.001, 1.097] | 73.797 | 147.594 | 0.214 | 0.001996 |

All three \(p\)-values equal \(1/(500+1)\), the smallest attainable value with
500 bootstrap replicates. The fixed-crossover model is therefore rejected for
all three conditions. The T-cell bootstrap interval also reaches the imposed
lower bound of 1.001. We consequently call \(\hat\mu\) an *effective
fixed-arrest-scale exponent* and use the fitted process as a deliberately
low-dimensional search surrogate, not as evidence that the cells execute an
exact Lévy walk.

## 4. Finite-search simulations

### 4.1 Normalization and periodic domain

All vectors and lengths in a condition were divided by
\(c_{\mathrm{arrest}}\):

\[
\mathbf r'=\frac{\mathbf r}{c_{\mathrm{arrest}}},
\qquad
\ell'=\frac{\ell}{c_{\mathrm{arrest}}}.
\]

The normalized detection radius was \(d'=1\), and the cubic torus side was

\[
W'=2L'_{\max}=\frac{2L_{\max}}{c_{\mathrm{arrest}}}.
\]

Thus the physical torus side was \(W=2L_{\max}\). This is a constructed model
domain, not an estimate of the microscope field of view, lymph-node size, or
in vivo territory available to a cell.

| Condition | \(d'\) | Physical \(d\), \(\mu\mathrm m\) | \(W'\) | Physical \(W\), \(\mu\mathrm m\) | Median normalized track budget [IQR] |
|---|---:|---:|---:|---:|---:|
| NK | 1 | 1.0 | 80.086 | 80.086 | 55.105 [29.940, 82.537] |
| Neutrophils | 1 | 0.5 | 533.786 | 266.893 | 160.541 [76.145, 281.071] |
| T cells | 1 | 0.5 | 295.189 | 147.594 | 237.233 [149.353, 401.867] |

For source track \(j\), the budget was the sum of the norms of its ordered
replay vectors, including the first and last runs:

\[
B'_j=\sum_k\lVert\mathbf r'_{jk}\rVert.
\]

### 4.2 Equal-\(\Delta\) targets

One target was centered in the periodic cube. The C-compatible target-membership
rules define a Ball detection region as a sphere of radius \(D'/2+d'\), a Disk
detection region as a right circular cylinder of radius \(D'/2+d'\) and
half-thickness \(d'\) in the \(z\) direction, and a Line detection region as a
radius-\(d'\) capsule around a segment oriented along \(y\). We set \(d'=1\).
The Disk rule is the cylinder implemented by the simulator, rather than the
Euclidean tubular neighborhood of an ideal planar disk with a rounded rim.
We use \(\Delta\) for the maximum face-on projected area of this **full
detection region**. It is not the intrinsic surface area of the nominal target.

For nominal target dimension \(D'\)—diameter for the Ball and Disk, centre-line
length for the Line—the matching equations at \(d'=1\) are

\[
\Delta_{\mathrm{Ball}}
=\Delta_{\mathrm{Disk}}
=\pi\left(\frac{D'}{2}+1\right)^2,
\qquad
\Delta_{\mathrm{Line}}=2D'+\pi.
\]

The dimensions used in simulation were therefore

\[
D'_{\mathrm{Ball}}=D'_{\mathrm{Disk}}
=2\left(\sqrt{\frac{\Delta}{\pi}}-1\right),
\qquad
D'_{\mathrm{Line}}=\frac{\Delta-\pi}{2}.
\]

The physical projected area is condition dependent,

\[
\Delta_{\mathrm{physical}}
=\Delta\,c_{\mathrm{arrest}}^2.
\]

Consequently, equal horizontal coordinates in the three panels are equal in
model units but not always equal in square micrometers. Each target had to
satisfy \(D'+2d'\leq W'\), preventing overlap with its periodic images. The
Line is the limiting shape, so the candidate grid was truncated separately for
each condition.

| Condition | Retained normalized \(\Delta\) | Physical area range, \(\mu\mathrm m^2\) | Maximum admissible line \(\Delta\) |
|---|---|---:|---:|
| NK | 4, 8, 16, 32, 64, 125 | 4–125 | 159.313 |
| Neutrophils | 4, 8, 16, 32, 64, 125, 250, 375, 500 | 1–125 | 1,066.713 |
| T cells | 4, 8, 16, 32, 64, 125, 250, 375, 500 | 1–125 | 589.519 |

Matching \(\Delta\) does not match detection-neighborhood volume. In
particular, Ball, Disk, and Line curves at a common \(\Delta\) compare shapes at
equal projected detection area, not at equal capture volume.

### 4.3 Trajectory experiments

Four finite trajectory constructions were simulated and saved. The paper-facing
figure displays the last two, which directly compare the randomized empirical
track with the fitted model.

| Construction | Simulation rule | Structure retained | Displayed? |
|---|---|---|---:|
| Exact orientation | Replay the recorded ordered directed-run vectors | Lengths, order, directions, turns, endpoint count, and laboratory orientation | No |
| Whole-track rotation | Apply one independent Haar-uniform 3-D rotation to the complete track in each replay | All internal track geometry, but not laboratory orientation | No |
| Independent rotations | Apply an independent Haar-uniform 3-D rotation to every recorded run vector in each replay | Exact ordered lengths, boundary runs, endpoint count, and focal-track budget; no directional correlations | Yes, solid |
| Fitted model | Draw run lengths IID from the condition-specific fitted density and directions IID uniformly on the sphere | The fixed-\(c\) marginal length model and the focal-track budget | Yes, dashed |

The solid-versus-dashed comparison does not isolate \(\mu\) alone. The fitted
construction also replaces track-specific length heterogeneity, length order,
and the recorded endpoint count by IID condition-level draws. The two
constructions additionally differ at the terminal budget boundary. The
independently rotated path uses all recorded vectors and reaches its budget.
For the fitted process, synthetic runs were generated until their cumulative
length covered the budget, but a final relocation that exceeded the remaining
budget was discarded rather than shortened, redrawn, or tested. The synthetic
path can therefore stop short of \(B'_j\), potentially by a substantial
amount for a short budget.

The replay used the vectorized Python implementation with target-membership
rules matched to the existing C geometry; the legacy C executable was not used
for this empirical-path analysis.

### 4.4 Detection rule, Monte Carlo sampling, and uncertainty

For every source track, trajectory construction, target shape, and value of
\(\Delta\), we performed 10,000 replays. Initial positions were drawn uniformly
from the torus, and the same generated path was tested against all relevant
targets. Positions were wrapped periodically, and target distance used the
minimum-image convention.

Detection was endpoint-only. The initial position was not tested, and a path
crossing a target between endpoints without ending within its detection
neighborhood was not detected. A replay was a success if at least one
budget-admissible endpoint lay in the target neighborhood.

Let \(\widehat p_j\) be the detected fraction among the 10,000 replays of source
track \(j\). The plotted condition estimate is

\[
\widehat p
=\frac{1}{J}\sum_{j=1}^{J}\widehat p_j.
\]

Thus every source track has equal weight in the detection probability even
though every run has equal weight in the fit. Shaded bands are pointwise 95%
percentile intervals from 1,000 whole-track bootstrap resamples. They are not
simultaneous confidence bands and are not adjusted for the multiple target
grid points. The bootstrap treats each per-track Monte Carlo fraction as fixed
and does not propagate its finite-replay binomial uncertainty; a separate
Monte Carlo standard-error column is retained in the summary table but is not
shown in the figure. A plotted point represents 270,000 trials for NK cells,
2,240,000 for neutrophils, or 1,630,000 for T cells. The two displayed
trajectory classes comprise 218.7 million path–target evaluations in total;
paths were reused across target shapes and sizes.

## 5. Figure and caption

[![Finite detection probability with three compact panels in one row](../outputs/arrest_scale_selected_area_volume/figures/finite_detection_probability_equal_area_compact_row.png)](../outputs/arrest_scale_selected_area_volume/figures/finite_detection_probability_equal_area_compact_row.pdf)

**Figure 1 | Finite endpoint-detection probability for targets matched by
projected detection-neighborhood area.** Colors and markers identify Ball
(blue circles), Disk (orange squares), and Line (green triangles) targets.
Solid curves show the empirical construction in which every recorded
directed-run vector was independently Haar-rotated, preserving the focal
track's ordered run lengths, number of endpoints, and finite travel budget.
Dashed curves show IID isotropic relocations sampled from the corresponding
fixed-\(c_{\mathrm{arrest}}\) fitted model under that track's budget. Points are
equal-track means from 10,000 uniformly initialized replays per track and
target; shading gives pointwise 95% intervals from 1,000 whole-track bootstrap
resamples. Detection was evaluated only at relocation endpoints in a periodic
cube with \(d'=1\) and \(W'=2L'_{\max}\); the initial position was not tested.
The horizontal variable \(\Delta\) is the normalized face-on projected area of
the full detection neighborhood, with
\(\Delta_{\mathrm{physical}}=\Delta c_{\mathrm{arrest}}^2\). Axes are linear,
and each panel has its own horizontal and vertical range to make the
within-condition comparison legible.

The vector PDF is available by clicking the figure. A common-y-axis version,
useful when the absolute scale across panels is important, is shown below.

[![Finite detection probability with a common y axis in a one-column layout](../outputs/arrest_scale_selected_area_volume/figures/finite_detection_probability_equal_area_one_column_shared_y.png)](../outputs/arrest_scale_selected_area_volume/figures/finite_detection_probability_equal_area_one_column_shared_y.pdf)

**Figure S1 | The data in Figure 1 shown with one common vertical scale.** All
symbols, trajectory constructions, target definitions, simulations, and
intervals are identical to Figure 1. The common y-axis range shows that the
absolute simulated probabilities are largest for the NK setup. This visual
difference should not be interpreted as a causal cell-type effect because the
conditions have different torus sizes, sampling intervals, physical target
scales, run-length distributions, and recorded-budget distributions.

**Layout note.** Figure 1 is a compact \(1\times3\) composition rendered
directly at the PNAS single-column width of 8.7 cm (3.425 inches), rather than
a reduced full-width figure. Every panel retains its own linear horizontal and
vertical scales; a shared vertical title, three ticks per axis, inset
scientific multipliers, and two legend rows keep the approximately 2.2-cm-wide
panels legible. All final-size text is at least 6 pt. A
[stacked panel-specific version](../outputs/arrest_scale_selected_area_volume/figures/finite_detection_probability_equal_area_one_column.pdf)
is available if a larger plotting area is preferred. The original
[panel-specific](../outputs/arrest_scale_selected_area_volume/figures/finite_detection_probability_equal_area.pdf)
and [common-axis](../outputs/arrest_scale_selected_area_volume/figures/finite_detection_probability_equal_area_shared_y.pdf)
full-width PDFs are also retained.

## 6. Results and interpretation

Detection probability increased monotonically with \(\Delta\), as expected when
the target neighborhood becomes larger. At a common projected area, Ball
targets generally had larger probabilities than Disk or Line targets. This is
not a pure shape-at-fixed-volume result: equal \(\Delta\) gives the Ball a
different three-dimensional detection-neighborhood volume from the
lower-dimensional Disk and Line.

The independently rotated empirical construction produced a larger point
estimate than the fitted model in all 72 displayed
condition-by-shape-by-\(\Delta\) cells. In paired whole-track bootstrap contrasts
of fitted minus independently rotated probability, the 95% interval lay
entirely below zero in 10 of 18 NK cells, 22 of 27 neutrophil cells, and all 27
T-cell cells (59 of 72 overall). These counts summarize an unadjusted target
grid; the grid cells are neither independent biological replicates nor a
familywise-error-controlled hypothesis test.

At the largest displayed \(\Delta\), the estimated simulated probabilities and their
ratios were:

| Condition and \(\Delta\) | Shape | Independent rotations | Fitted model | Independent / fitted |
|---|---|---:|---:|---:|
| NK, \(\Delta=125\) | Ball | \(9.1074\times10^{-3}\) | \(8.5704\times10^{-3}\) | 1.063 |
|  | Disk | \(3.8519\times10^{-3}\) | \(3.3556\times10^{-3}\) | 1.148 |
|  | Line | \(4.4815\times10^{-3}\) | \(3.4778\times10^{-3}\) | 1.289 |
| Neutrophils, \(\Delta=500\) | Ball | \(4.1027\times10^{-4}\) | \(2.0759\times10^{-4}\) | 1.976 |
|  | Disk | \(1.1741\times10^{-4}\) | \(4.5536\times10^{-5}\) | 2.578 |
|  | Line | \(1.2054\times10^{-4}\) | \(4.5536\times10^{-5}\) | 2.647 |
| T cells, \(\Delta=500\) | Ball | \(2.9337\times10^{-3}\) | \(1.9067\times10^{-3}\) | 1.539 |
|  | Disk | \(7.1104\times10^{-4}\) | \(3.2883\times10^{-4}\) | 2.162 |
|  | Line | \(7.4785\times10^{-4}\) | \(3.3067\times10^{-4}\) | 2.262 |

The direction of the discrepancy is therefore consistent: replacing the exact
ordered empirical lengths and endpoint count by the condition-level fitted
process reduces endpoint encounter probability under the chosen finite-budget
protocol. The discrepancy is modest for NK and more pronounced for the
neutrophil and T-cell panels, especially for Disk and Line targets. This is
consistent with, but not uniquely attributable to, the rejected marginal
run-length fit. Track-to-track heterogeneity, length ordering, the altered
number of endpoints, and terminal-budget undershoot can all contribute.

The strongest defensible conclusion is consequently a *model-checking*
conclusion: the theoretical fitted process makes a qualitative prediction that
detection rises with projected target area and usually remains on the same
order of magnitude in these finite replays, but it has lower point estimates
than the randomized-empirical construction in every displayed cell. The experiment
grounds the theory in measured displacement scales and budgets while also
identifying where the one-parameter marginal surrogate loses empirical search
performance.

## 7. Limitations

1. **Counterfactual targets.** LTDB contains cell trajectories but no paired
   target coordinates or observed detection events. The analysis tests the
   model on empirical movement paths, not on biological encounter outcomes.
2. **Rejected run-length model.** The fixed-\(c_{\mathrm{arrest}}\) density
   failed all three conditional IID-run KS checks. The checks do not reproduce
   within-track dependence, and the reported exponents are sensitivity
   summaries rather than validated Lévy exponents.
3. **Operational scales.** Neither \(c_{\mathrm{arrest}}\) nor \(W=2L_{\max}\)
   is a measured tissue scale. The former derives from an arrest cutoff; the
   latter is a constructed finite domain.
4. **Bundled trajectory contrast.** Solid versus dashed curves differ in more
   than \(\mu\): they also differ in track-specific length heterogeneity, length
   order, endpoint count, and terminal use of the finite budget.
5. **Intermittent detection.** Continuous segment intersections are not tested.
   Conclusions apply to endpoint-only detection and need not transfer to
   continuous sensing.
6. **Cross-condition comparability.** Absolute probabilities should not be
   interpreted as intrinsic NK-versus-neutrophil-versus-T-cell performance.
   The experiments differ in cadence, stimulus, torus side, physical target
   area, and budget distribution.
7. **Uncertainty hierarchy.** Whole-track bootstrapping accounts for repeated
   runs from a cell, but not possible dependence among cells recorded in the
   same video or animal or finite Monte Carlo error within each track estimate.
   Very small probabilities, particularly at small \(\Delta\) for
   neutrophils, are also supported by few detected replays.
8. **T-cell pooling.** Four CMTMR comparator tracks were pooled with 159 HIV-GFP
   tracks. Their contribution is small but should be disclosed; an
   infected-only rerun is an appropriate sensitivity analysis if this panel is
   used for a biological claim specific to infected cells.
9. **Corrected LTDB ground truth.** `LTDB020_a` contributes one NK track even
   though later LTDB work reports errors in the `LTDB020` ground truth [11].
   The diagnostic deletion above is small, but the publication analysis should
   exclude this file in a formal rerun or justify its retention.

## 8. Reproducibility record

The figure is generated from the saved, versioned artifacts below:

- [Fixed-crossover fit table](../outputs/arrest_scale_selected_area_volume/fit/fixed_arrest_fits.csv)
- [Target geometry manifest](../outputs/arrest_scale_selected_area_volume/targets/target_manifest.csv)
- [Per-track inputs and budgets](../outputs/arrest_scale_selected_area_volume/inputs/track_inputs.csv)
- [Per-track replay summaries](../outputs/arrest_scale_selected_area_volume/finite_replay/finite_track_summaries.csv)
- [Condition-level finite summary](../outputs/arrest_scale_selected_area_volume/finite_replay/finite_summary.csv)
- [Paired finite contrasts](../outputs/arrest_scale_selected_area_volume/finite_replay/finite_contrasts.csv)
- [Exact simulation settings](../outputs/arrest_scale_selected_area_volume/00_plan/settings.json)
- [Analysis and plotting implementation](../ltdb_levy/arrest_scale_experiment.py)

The random seed was 731,905. Each output table has a companion metadata JSON,
and the completion record stores the configuration signature. The full replay
used 64 worker processes. The plotting function reads the saved finite summary
and fit table, so cosmetic figure changes do not require rerunning the Monte
Carlo trajectories. Both PDF and PNG variants can be regenerated with

~~~bash
PYTHONPATH=. .venv/bin/python -c \
  "from ltdb_levy.arrest_scale_experiment import make_arrest_scale_figures; make_arrest_scale_figures()"
~~~

## References

1. Pizzagalli DU, *et al.* Leukocyte Tracking Database, a collection of immune
   cell tracks from intravital 2-photon microscopy videos. *Scientific Data*.
   2018;5:180129. [doi:10.1038/sdata.2018.129](https://doi.org/10.1038/sdata.2018.129).
2. Pizzagalli DU, *et al.* Leukocyte Tracking Database data collection.
   Figshare. [doi:10.6084/m9.figshare.c.3827890](https://doi.org/10.6084/m9.figshare.c.3827890).
3. Farsakoglu Y, *et al.* Influenza vaccination induces NK-cell-mediated type-II
   IFN response that regulates humoral immunity in an IL-6-dependent manner.
   *Cell Reports*. 2019;26:2307–2315.e5.
   [doi:10.1016/j.celrep.2019.01.104](https://doi.org/10.1016/j.celrep.2019.01.104).
4. Pizzagalli DU, *et al.* Characterization of the dynamic behavior of
   neutrophils following influenza vaccination. *Frontiers in Immunology*.
   2019;10:2621.
   [doi:10.3389/fimmu.2019.02621](https://doi.org/10.3389/fimmu.2019.02621).
5. Murooka TT, *et al.* HIV-infected T cells are migratory vehicles for viral
   dissemination. *Nature*. 2012;490:283–287.
   [doi:10.1038/nature11398](https://doi.org/10.1038/nature11398).
6. Georgantzoglou A, *et al.* A two-step search and run response to gradients
   shapes leukocyte navigation in vivo. *Journal of Cell Biology*.
   2022;221:e202103207.
   [doi:10.1083/jcb.202103207](https://doi.org/10.1083/jcb.202103207).
7. Ganusov VV, Zenkov VS, Majumder B. Correlation between speed and turning
   naturally arises for sparsely sampled cell movements. *Physical Biology*.
   2023;20:016002.
   [doi:10.1088/1478-3975/acb18c](https://doi.org/10.1088/1478-3975/acb18c).
8. Lindquist RL, *et al.* Visualizing dendritic cell networks in vivo.
   *Nature Immunology*. 2004;5:1243–1250.
   [doi:10.1038/ni1139](https://doi.org/10.1038/ni1139).
9. Lee BJ, Hegewisch Solloa E, Shannon MJ, Mace EM. Generation of cell-derived
    matrices that support human NK cell migration and differentiation.
    *Journal of Leukocyte Biology*. 2020;108:1369–1378.
    [doi:10.1002/JLB.1MA0420-635R](https://doi.org/10.1002/JLB.1MA0420-635R).
10. Clauset A, Shalizi CR, Newman MEJ. Power-law distributions in empirical
    data. *SIAM Review*. 2009;51:661–703.
    [doi:10.1137/070710111](https://doi.org/10.1137/070710111).
11. Pizzagalli DU, *et al.* CANCOL, a computer-assisted annotation tool to
    facilitate colocalization and tracking of immune cells in intravital
    microscopy. *Journal of Immunology*. 2022;208:1493–1499.
    [doi:10.4049/jimmunol.2100811](https://doi.org/10.4049/jimmunol.2100811).
