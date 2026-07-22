# Empirical Assessment of Lévy-Walk Surrogates in LTDB Leukocyte Trajectories

> **Superseded fit snapshot.** The numerical fit tables in this document were
> generated with the former 1 µm crossover. The primary configuration now
> estimates one crossover per condition, using 3.75 µm (half of a 7.5 µm
> reference cell diameter) only as the optimizer start. Current fit results are
> in `outputs/primary/fit/mixture_fits.csv`. This report will be regenerated
> after the corresponding replay-geometry scale is fixed.

## Abstract

This report describes an empirical analysis of leukocyte trajectories from the
Leukocyte Tracking Database (LTDB). The primary objective was to determine
whether directed relocations extracted from the trajectories are adequately
described by the flat-plus-power-law step distribution used by the accompanying
C simulator, and whether fitted exponents could subsequently be evaluated in a
shape-controlled search task.

The analysed corpus comprised 44 ground-truth trajectory files, 44,722
space-time observations, and 728 raw tracks. After prespecified quality control,
44,657 observations and 710 tracks were retained. Under the primary
90-degree turning rule, these produced 17,981 candidate directed runs, of
which 16,578 were complete after boundary-run exclusion. Sixteen biological
and sampling conditions were identified; eight met the minimum requirement of
15 tracks and 100 directed runs and together contributed 15,773 runs.

The fitted mixture distribution was rejected by each separate, unadjusted
condition-level parametric-bootstrap Kolmogorov–Smirnov test
(\(p=0.001996\) in all eight conditions). A truncated lognormal distribution
achieved the best track-held-out predictive score in all eight conditions.
Consequently, the estimated power-law parameters are reported only as
**effective exponents**; they do not constitute evidence for an exact Lévy
walk. This conclusion was preserved across the viable
segmentation-sensitivity analyses.

The primary fit imposed no positive speed cutoff. Sensitivity analyses treated
steps below 1, 2, or 4 \(\mu\)m min\(^{-1}\) as pauses. The conventional
2-\(\mu\)m min\(^{-1}\) arrest threshold changed the mean effective exponent
from 1.349 to 1.272 across the eight eligible conditions, but it did not change
the model-adequacy conclusion: every fit was rejected. The simulation and
finite-budget search stages are deliberately deferred until the empirical fit
has been reviewed.

## 1. Content of the data

### 1.1 Source and scientific content

The data originate from the Leukocyte Tracking Database described by
Pizzagalli *et al.* [1]. LTDB contains manually curated leukocyte trajectories
obtained from intravital two-photon microscopy. The present analysis used the
consensus ground-truth tracking files in `datasets/GT_TRACKS/`.

The selected corpus contains:

- 44 biological ground-truth files;
- 44,722 instantaneous position annotations;
- 728 raw cell tracks;
- four tracked cell classes: B cells, T cells, neutrophils, and natural killer
  cells;
- two imaging sites: popliteal lymph node and spleen;
- five experimental states or stimuli: influenza vaccine, ovalbumin, vaccinia
  virus, steady state, and HIV-infected humanized T-cell experiments; and
- file-specific sampling intervals ranging from 12 to 61 seconds.

The `CS` files are selected case-study regions, whereas the `LTDB` files are the
standard collection. Filename suffixes `_a` and `_b` denote distinct tracked
cell populations and are not channel identifiers.

Biological metadata were obtained from Tables 3 and 4 of the LTDB data
descriptor and cross-checked against the deposited SQL database [2]. All 44
metadata rows used here have explicit provenance and are marked as verified.
The distribution of files across the principal biological categories was:

| Cell type | Imaging site | Stimulus/state | Files |
|---|---|---|---:|
| B cells | Popliteal lymph node | Ovalbumin | 2 |
| B cells | Spleen | Ovalbumin | 1 |
| Natural killer cells | Popliteal lymph node | Influenza vaccine | 4 |
| Neutrophils | Popliteal lymph node | Influenza vaccine | 23 |
| Neutrophils | Spleen | Ovalbumin | 3 |
| Neutrophils | Spleen | Vaccinia virus | 2 |
| T cells | Popliteal lymph node | HIV-infected humanized T cell | 7 |
| T cells | Popliteal lymph node | Steady state | 2 |

### 1.2 File structure and units

Each semicolon-delimited ground-truth file was interpreted positionally:

1. line 1 contains `VideoID;dx;dy;dz;dt`;
2. line 2 contains channel-visibility flags and is not an observation row; and
3. line 3 onward contains `TrackID;x;y;z;t`.

The spatial coordinates \(x\), \(y\), and \(z\) are already expressed in
micrometres and were therefore not rescaled by the voxel dimensions in the
header. The variable \(t\) is a frame index. Physical time was calculated as
\(t \times dt\), where \(dt\) is the file-specific frame interval in seconds.

Raw track duration was heterogeneous: tracks contained between 2 and 241
observations, with a median of 45; 66 raw tracks contained fewer than 10
observations. This heterogeneity motivated both a minimum track-length rule and
track-balanced statistical weighting.

### 1.3 Audited source anomalies

Source anomalies were surfaced rather than silently corrected:

- `SQUARE_GT..csv` is a synthetic fixture and was excluded from the biological
  corpus by configuration;
- the internal identifier in `LTDB004_b_GT.csv` is `LTDB005_b`, whereas its
  filename identifies it as `LTDB004_b`; the filename was used and the mismatch
  was recorded; and
- the channel row in `CS003_a` is entirely zero, but this does not alter the
  positional interpretation of its trajectory rows.

No parse failures occurred among the 44 selected biological files.

## 2. Preprocessing and construction of directed relocations

### 2.1 Track-level quality control

Observations were sorted by track and frame using a stable ordering. Duplicate
track-frame records would have triggered the configured error policy; none were
present. Tracks would also have been split at every frame gap greater than one
to avoid converting missing observations into artificial long relocations; no
such gaps occurred in the selected files. No interpolation was performed.

Tracks with fewer than six observations were excluded. This removed 18 tracks
and 65 observations, leaving:

- 44,657 observations;
- 710 retained tracks/fragments; and
- 43,947 frame-to-frame displacement vectors.

Among the retained tracks, the number of observations ranged from 6 to 241,
with a median of 47; 48 retained tracks contained fewer than 10 observations.

For each frame displacement \(\mathbf{u}_i\), the pipeline calculated its
three-dimensional vector, Euclidean length, elapsed time, and speed. Turning
angles were evaluated as

\[
\theta_i =
\operatorname{atan2}
\left(
\lVert \mathbf{u}_{i-1}\times\mathbf{u}_i\rVert,
\mathbf{u}_{i-1}\cdot\mathbf{u}_i
\right),
\]

which is numerically stable near both zero and \(\pi\). There were 5,776
zero-length displacements; these were retained in the audit and treated as
nonmoving boundaries during run segmentation.

### 2.2 Directed-run segmentation

The primary relocation unit was a directed run rather than an individual frame
displacement. The primary segmentation used a turning-angle threshold of
90 degrees and no positive speed cutoff. This choice was anchored to an in vivo
leukocyte study by Georgantzoglou *et al.* [3], which classified successive
neutrophil movements below 45 degrees as persistent, movements between 45 and
90 degrees as small turns, and movements between 90 and 180 degrees as large
turns or U-turns. Thus, the primary rule ends a directed run only when the
measured change in direction enters that study's “large-turn” category. The
published angle construction is related to, but not identical with, the
consecutive-vector angle used here. Accordingly, 90 degrees is treated as a
literature-informed operational definition rather than a universal biological
constant.

This qualification is important because turning-angle distributions depend on
the sampling interval [4,5], and another cell-migration segmentation study that
used a 45-degree threshold explicitly described it as heuristic and potentially
cell-type specific [6]. Conditions with different frame intervals were
therefore kept separate, and 45- and 60-degree thresholds were retained as
prespecified sensitivity analyses. No smoothing was applied before angle
calculation.

A displacement whose turning angle strictly exceeded the threshold started the
next run. This convention assigns the turning displacement deterministically
to the new run. Nonmoving steps, fragment boundaries, and—if present—frame gaps
also terminated a run.

For each run, two distances were retained:

- **run length:** the Euclidean norm of the net displacement from the first to
  the last endpoint; and
- **path length:** the left-to-right cumulative sum of its constituent
  displacement lengths.

Straightness was stored as run length divided by path length. The primary
segmentation generated 17,981 runs. The first and last run of each track
fragment were treated as potentially field-of-view censored and excluded from
the primary fit. This excluded 1,403 boundary runs and left 16,578 included
runs. No runs were excluded for failing the configured minimum run length or
frame count. Because one-step runs were permitted, 9,620 included runs
(58.0%) consisted of a single frame-to-frame displacement; the possible
contribution of frame-scale positional noise must therefore be considered when
interpreting the fitted distributions.

### 2.3 Conditions, pools, and statistical weighting

A condition was defined by the tuple

\[
(\text{cell type},\ \text{imaging site},\ \text{stimulus},\ dt).
\]

Including \(dt\) prevents trajectories sampled at materially different time
intervals from being pooled as though they had the same observational
resolution. Sixteen conditions were present. A condition was eligible for
distribution fitting only if it contributed at least 15 tracks and at least
100 included runs; eight conditions met both requirements and contributed
15,773 runs. The neutrophil/spleen/vaccinia-virus condition at 32 seconds,
which was eligible under the more aggressive 45-degree splitting rule, was
ineligible at 90 degrees because only 13 tracks contributed complete directed
runs, below the prespecified minimum of 15.

Within each eligible condition, the primary empirical relocation pool used
cell-balanced weights. If track \(j\) contributed \(n_j\) runs and the condition
contained \(J\) tracks, each run from that track received weight

\[
w_{ij} = \frac{1}{J n_j}.
\]

Thus, every track contributed equal total probability mass irrespective of its
duration. This avoids allowing a small number of long tracks to dominate the
estimated relocation distribution.

### 2.4 Fitted step model

The fitted parametric model reproduced the step distribution used by the
existing C simulator:

\[
f(\ell\mid\mu,L,c)
=
a(\mu,L,c)\,
\min\left[1,\left(\frac{\ell}{c}\right)^{-\mu}\right],
\qquad 0\leq\ell\leq L.
\]

The crossover \(c\) was fixed at 1 \(\mu\)m and was not estimated. The upper
support \(L\) was the empirical maximum within each eligible condition.
Only the exponent \(\mu\) was fitted. Because the crossover is imposed rather
than measured, the analysis also evaluated the sensitivity of the fitted
exponent to alternative crossover values.

## 3. Tests and validation procedures

### 3.1 Parametric goodness-of-fit test

Model adequacy was assessed independently in every eligible condition with a
parametric-bootstrap Kolmogorov–Smirnov test. The observed KS statistic compared
the weighted empirical distribution with the fitted C-mixture distribution.
For each of 500 bootstrap replicates, a new sample was generated from the
fitted model, \(\mu\) was refitted, and the KS statistic was recalculated. The
reported \(p\)-value therefore accounts for parameter estimation. With 500
replicates, the smallest attainable corrected \(p\)-value was
\(1/(500+1)=0.001996\).

A KS power analysis against the best-fitting alternative family was also
generated at multiple sample sizes. Its purpose was to prevent a failure to
reject from being interpreted as support for the mixture when statistical power
was low. In the present analysis, however, the mixture was rejected rather than
merely not rejected.

### 3.2 Uncertainty in the fitted exponent

Uncertainty intervals for \(\mu\) were obtained from 1,000 clustered bootstrap
replicates. Whole tracks—not individual runs—were resampled, and all runs from a
selected track were retained together. This preserves within-track dependence
and treats the cell trajectory as the independent biological unit.

### 3.3 Predictive model comparison

The C-mixture was compared with four alternative distributions on the same
support:

- truncated exponential;
- truncated lognormal;
- truncated Weibull; and
- power law with an exponential cutoff.

Comparison used 10-fold grouped cross-validation, with entire tracks assigned
to held-out folds. Mean held-out log density was the scoring rule. No run from a
held-out track was allowed to appear in its training data.

### 3.4 Segmentation sensitivity

Segmentation choices were specified before examining the fitted exponents. The
analysis repeated segmentation and fitting across turning-angle thresholds of
45, 60, and 90 degrees and speed cutoffs of none, 1, 2, and 4
\(\mu\)m min\(^{-1}\). In the software configuration, where speed is expressed
in \(\mu\)m s\(^{-1}\), these values are respectively `null`, 0.0167, 0.0333,
and 0.0667. Eligibility thresholds were reapplied independently in each grid
cell, and every eligible fit received its own refitted parametric-bootstrap KS
test.

The central 2-\(\mu\)m min\(^{-1}\) value is widely used in intravital
immune-cell studies to define an “arrest coefficient”: the fraction of
instantaneous observations below that speed [7]. That convention does not, by
itself, establish that every slow step is a boundary between biological runs.
The speed-aware analyses here are therefore explicitly operational
sensitivities. A step is retained as moving when its measured speed is greater
than or equal to the selected cutoff. A slower step is omitted from the run,
terminates any run before it, and prevents connection to a run after it. The
1- and 4-\(\mu\)m min\(^{-1}\) values bracket the conventional cutoff by a
factor of two. The primary analysis retains every non-zero step and uses no
positive speed threshold.

This scale is compatible with the observed measurements. Of 43,947
frame-to-frame displacements, 5,776 were exactly zero. Among the 38,171
non-zero displacements, the median speed was 0.1003 \(\mu\)m s\(^{-1}\)
(6.02 \(\mu\)m min\(^{-1}\)); 5.1%, 14.4%, and 33.5% were below 1, 2,
and 4 \(\mu\)m min\(^{-1}\), respectively. By contrast, the earlier candidate
settings of 0.5 and 1.0 \(\mu\)m s\(^{-1}\)—30 and 60 \(\mu\)m
min\(^{-1}\)—would have labelled 98.0% and 99.9% of non-zero displacements as
slow. They were therefore removed before this fit rather than presented as
biologically meaningful sensitivities.

With no speed cutoff, 9, 9, and 8 conditions were eligible at 45, 60, and
90 degrees, respectively. Across the 1-, 2-, and 4-\(\mu\)m min\(^{-1}\)
analyses, the corresponding eligible-condition counts were 9/9/8 at
45 degrees, 8/8/8 at 60 degrees, and 8/8/8 at 90 degrees.

### 3.5 Replay-engine validation

The finite-budget replay engine was validated independently of the biological
analysis.

1. **Analytic one-step probability.** For a single endpoint generated from a
   uniform starting position on the torus, translation invariance gives the
   exact probability
   \[
   P(\text{detection}) =
   \frac{\text{volume of the target neighbourhood}}
        {\text{domain volume}}.
   \]
   Monte Carlo estimates for ball, disk, and line targets were required to fall
   within four binomial standard errors of this analytic value.
2. **Geometry and boundary tests.** Minimum-image wrapping was tested at exact
   side-length multiples, half-box boundaries, and small negative coordinates.
   Ball, disk, and line membership included adversarial boundary points.
3. **Endpoint and budget semantics.** Tests verified endpoint-only detection,
   omission of the initial position, left-to-right cumulative distance, and the
   distinction between discarding and truncating an overshooting step.
4. **Control construction.** Tests verified empirical length preservation,
   block preservation, per-vector random rotation, isotropic direction
   sampling, pool-support containment, and the intended cell-balanced and
   step-weighted probabilities.
5. **MLE recovery.** Synthetic samples of at least 20,000 draws were generated
   at known values of \(\mu\), including values near the analytic
   \(\mu=1\) limit, and the fitted exponent was required to recover the
   generating value within its expected numerical/statistical tolerance.

The complete automated test suite contains 48 tests; all 48 passed. Python
bytecode compilation and whitespace checks on the files changed for this
analysis also completed successfully.

### 3.6 Deferred simulation and search validation

The implemented replay compares six trajectory classes:

| Class | Relocation lengths | Directions/order retained |
|---|---|---|
| Empirical | Observed | Fully observed trajectory |
| Block resampled | Observed blocks | Short-range order within blocks |
| Ordered-length rotated | Observed and ordered | Length order only; directions independently rotated |
| Pool isotropic | Empirical pool | Neither direction nor temporal order |
| Fitted \(\hat{\mu}\) | C-mixture | Isotropic and unordered |
| \(\mu=2\) | C-mixture | Isotropic and unordered |

Every source track will supply its own finite path-length budget. Detection
will be checked only at relocation endpoints, and an overshooting step will be
discarded. Common random numbers will be shared across trajectory classes
within a replicate to improve the precision of paired contrasts.

No replay or C-simulator result is included in the present fit-focused
analysis. Earlier exploratory replay artifacts were generated under a
different primary segmentation and are intentionally not interpreted here.
The target grid must be reviewed and the simulations rerun using the accepted
90-degree fit before any search-performance or Python–C claim is made.

## 4. Outcomes and interpretation

### 4.1 Distributional fits

All eight eligible conditions fitted numerically. Results are shown below.
“PLN” denotes popliteal lymph node.

| Cell type; site; stimulus; \(dt\) | Tracks | Runs | \(\hat{\mu}\) | Clustered 95% CI | KS \(p\) | Best held-out model |
|---|---:|---:|---:|---:|---:|---|
| B cells; PLN; ovalbumin; 15 s | 64 | 2,167 | 1.281 | 1.009–1.353 | 0.0020 | Truncated lognormal |
| B cells; spleen; ovalbumin; 30 s | 26 | 1,226 | 1.678 | 1.411–1.817 | 0.0020 | Truncated lognormal |
| Natural killer cells; PLN; influenza vaccine; 30 s | 27 | 517 | 1.548 | 1.290–1.718 | 0.0020 | Truncated lognormal |
| Neutrophils; PLN; influenza vaccine; 12 s | 46 | 1,367 | 1.391 | 1.162–1.463 | 0.0020 | Truncated lognormal |
| Neutrophils; PLN; influenza vaccine; 15 s | 224 | 5,829 | 1.487 | 1.347–1.515 | 0.0020 | Truncated lognormal |
| Neutrophils; PLN; influenza vaccine; 60 s | 44 | 569 | 1.122 | 1.001–1.184 | 0.0020 | Truncated lognormal |
| T cells; PLN; HIV-infected humanized T cell; 15 s | 163 | 3,899 | 1.138 | 1.035–1.181 | 0.0020 | Truncated lognormal |
| T cells; PLN; steady state; 20 s | 22 | 199 | 1.142 | 1.001–1.203 | 0.0020 | Truncated lognormal |

The point estimates ranged from 1.122 to 1.678. This numerical range must not
be interpreted in isolation: the fitted distribution was rejected in every
condition. The estimated \(\mu\) values therefore summarize the best
one-parameter approximation *within an inadequate model*. They are not
estimates of a validated biological Lévy exponent.

The predictive comparison reinforced the goodness-of-fit result. The truncated
lognormal distribution ranked first in every condition. Relative to the best
alternative, the C-mixture's held-out mean log-density difference was negative
in all eight conditions, ranging from \(-0.336\) to \(-0.123\) log-density
units per held-out run. Thus, there was
no condition in which the mixture provided superior held-out predictive
performance.

The power diagnostic did not indicate that rejection depended on using the
full data volume. For every condition, all 200 simulated datasets drawn from
the fitted truncated lognormal alternative were rejected by the mixture-model
KS test at each examined sample size (25%, 50%, and 100% of the observed run
count). Thus, the estimated rejection probability was 1.00 even for the
smallest examined sample of 49 runs. This is a finite Monte Carlo estimate, not
a mathematical statement that power is exactly one.

### 4.2 Interpretation with respect to Lévy-walk behaviour

The results do not support the claim that the observed directed-run lengths are
generated by the exact flat-plus-power-law mixture implemented in the C
simulator. In particular:

- all goodness-of-fit tests rejected that distribution;
- an alternative truncated lognormal model predicted held-out tracks better in
  every eligible condition;
- the effective exponent varied materially across cell types, conditions, and
  sampling intervals; and
- several confidence intervals approached the lower fitting bound, indicating
  that the imposed parametric form was being stretched toward its boundary.

None of the eight clustered 95% intervals included \(\mu=2\). Even that
observation is secondary to model adequacy: a confidence interval for
\(\mu\) has a literal generative interpretation only if the assumed mixture is
an adequate description of the data. More generally, a fitted exponent near 2
could not establish Lévy optimality if the underlying step distribution were
misspecified.

The segmentation analysis altered the numerical estimates but not this central
conclusion. Without a speed cutoff, effective exponents ranged from
1.036 to 1.905 at 45 degrees, 1.006 to 1.882 at 60 degrees, and
1.122 to 1.678 at 90 degrees. All 100 eligible fits across the complete
angle-by-speed grid fitted numerically, and every one rejected the mixture
(\(p \leq 0.00599\)).

At the primary 90-degree angle, the speed results were:

| Minimum moving speed | Included runs, all conditions | Eligible conditions | Runs in eligible conditions | Mean \(\hat{\mu}\) | Range of \(\hat{\mu}\) | Fits at lower bound |
|---|---:|---:|---:|---:|---:|---:|
| None (primary) | 16,578 | 8 | 15,773 | 1.349 | 1.122–1.678 | 0 |
| 1 \(\mu\)m min\(^{-1}\) | 16,430 | 8 | 15,627 | 1.335 | 1.083–1.647 | 0 |
| 2 \(\mu\)m min\(^{-1}\) | 15,697 | 8 | 15,054 | 1.272 | 1.028–1.478 | 0 |
| 4 \(\mu\)m min\(^{-1}\) | 13,747 | 8 | 13,348 | 1.156 | 1.001–1.429 | 3 |

The mean in this table is an unweighted descriptive average of the eight
condition-specific estimates; it is not a pooled biological estimate. The
1-\(\mu\)m min\(^{-1}\) rule had little effect. The conventional
2-\(\mu\)m min\(^{-1}\) rule lowered the descriptive mean by 0.077, with
larger changes in the B-cell/spleen and natural-killer-cell conditions. The
4-\(\mu\)m min\(^{-1}\) rule had a qualitatively stronger effect and pushed
three estimates to the permitted lower bound of 1.001. This pattern means that
aggressive removal of slow displacements changes the numerical effective
exponent and can make the one-parameter fit unstable. It does not rescue the
mixture: the goodness-of-fit test rejected it in every speed analysis.

### 4.3 Dependence on the imposed crossover

The fitted model assumes that the density is flat below a crossover length
\(c\), and follows a power law above \(c\). The primary value \(c=1\)
\(\mu\)m came from the C simulator; it was not estimated from LTDB. Refitting
the same 90-degree runs at alternative crossover values produced:

| Crossover \(c\) (\(\mu\)m) | Mean \(\hat{\mu}\) | Range of \(\hat{\mu}\) | Fits at lower bound |
|---:|---:|---:|---:|
| 0.25 | 1.019 | 1.001–1.104 | 6 |
| 0.50 | 1.108 | 1.001–1.247 | 3 |
| 1.00 (primary) | 1.349 | 1.122–1.678 | 0 |
| 2.00 | 1.877 | 1.400–3.094 | 0 |
| 4.00 | 3.038 | 1.869–6.618 | 0 |

These are again unweighted descriptive averages across conditions. The strong,
monotonic dependence is expected because changing \(c\) changes which observed
lengths are assigned to the flat and power-law portions of the model. It shows
that \(\hat{\mu}\) is not a standalone property of the tracks: it is conditional
on the simulator's imposed crossover as well as on the segmentation rule. This
is an additional reason to call it an effective exponent.

### 4.4 Simulation status

The simulation stage has not been rerun under the 90-degree primary
segmentation. It is therefore outside the evidential scope of this version of
the report. No search-performance contrast, target-shape conclusion, fitted-
versus-\(\mu=2\) comparison, or Python–C agreement claim should be inferred
from the fit results. Those questions require a separately documented
simulation run based on the accepted fit configuration.

### 4.5 Overall conclusion

The completed distributional analysis provides evidence **against** the exact
C-mixture as a generative model for LTDB directed-run lengths. The estimated
parameters are useful as effective summaries and as inputs for controlled
simulation, but they should not be described as validated Lévy exponents.
Observed deviations could arise from persistence, pauses, correlations,
heterogeneous tissue structure, measurement resolution, field-of-view
censoring, or a different relocation-length family; the present analysis does
not identify a unique causal mechanism.

The search-performance question remains open. Resolving it requires a
scientifically justified domain and target grid that yields adequate event
counts without periodic self-overlap, followed by a new pilot, a frozen
confirmatory replay, and clustered uncertainty estimation. C experiments at
the accepted fitted parameters must also be executed before the independent
Python–C first-passage comparison can be reported.

### 4.6 Scope and limitations

The following limitations govern interpretation:

- LTDB contains trajectories but no biological targets or detection events;
  every search outcome is counterfactual.
- Sampling interval affects which displacements and turns can be observed.
- The 90-degree primary threshold is grounded in a related leukocyte analysis
  of large turns but is not a universal or independently validated biological
  boundary; the published angle definition is not identical to the one used
  here.
- No trajectory smoothing was applied, and 58.0% of included runs contained
  only one frame-to-frame displacement, making the segmentation potentially
  sensitive to frame-scale positional noise.
- The 2-\(\mu\)m min\(^{-1}\) value is established in the cited literature as
  an arrest-coefficient cutoff, not as a run-segmentation law; its use as a
  hard run boundary is only a sensitivity analysis.
- The 1 \(\mu\)m flat-to-power-law crossover is imposed rather than estimated.
- First and last directed runs are subject to field-of-view censoring.
- Pooling can conceal between-track and between-video heterogeneity, although
  conditions were stratified and tracks were equally weighted.
- Empirical pools cannot generate relocations outside their observed support.
- The confined \(z\)-range makes empirical directions non-isotropic by
  construction.
- Face-on projected-area matching is not orientation-fair: a randomly oriented
  disk has a smaller mean projection than its face-on area, so a future shape
  contrast would partly reflect this convention.
- The C engine estimates unbounded first passage, whereas the biological replay
  uses finite track-specific budgets. These are different estimands and should
  be compared only with that distinction stated explicitly.

## Reproducibility

The current fit used
[`configs/primary.yaml`](../configs/primary.yaml). Its source and configuration
hashes are recorded separately for preprocessing and fitting in
[`outputs/primary/preprocess/run_context.json`](../outputs/primary/preprocess/run_context.json)
and
[`outputs/primary/fit/run_context.json`](../outputs/primary/fit/run_context.json).
Detailed machine-readable fit results are available in:

- [`outputs/primary/preprocess/`](../outputs/primary/preprocess/);
- [`outputs/primary/fit/mixture_fits.csv`](../outputs/primary/fit/mixture_fits.csv);
- [`outputs/primary/fit/model_comparison.csv`](../outputs/primary/fit/model_comparison.csv);
  and
- [`outputs/primary/fit/segmentation_fit_sensitivity.csv`](../outputs/primary/fit/segmentation_fit_sensitivity.csv).

Replay, C-comparison, and generated-report directories may contain exploratory
artifacts from the earlier segmentation configuration. They are not evidence
for, and should not be combined with, the 90-degree results reported here.

## References

1. Pizzagalli DU, Farsakoglu Y, Palomino-Segura M, *et al.* Leukocyte Tracking
   Database, a collection of immune cell tracks from intravital 2-photon
   microscopy videos. *Scientific Data*. 2018;5:180129.
   [doi:10.1038/sdata.2018.129](https://doi.org/10.1038/sdata.2018.129).
2. Pizzagalli DU, *et al.* SQL Dump of the Database. Figshare; 2018.
   [doi:10.6084/m9.figshare.5946472.v1](https://doi.org/10.6084/m9.figshare.5946472.v1).
3. Georgantzoglou A, Poplimont H, Walker HA, Lämmermann T, Sarris M. A
   two-step search and run response to gradients shapes leukocyte navigation
   in vivo. *Journal of Cell Biology*. 2022;221(8):e202103207.
   [doi:10.1083/jcb.202103207](https://doi.org/10.1083/jcb.202103207).
4. Beltman JB, Marée AFM, de Boer RJ. Analysing immune cell migration.
   *Nature Reviews Immunology*. 2009;9:789–798.
   [doi:10.1038/nri2638](https://doi.org/10.1038/nri2638).
5. Ganusov VV, *et al.* Correlation between speed and turning naturally arises
   for sparsely sampled cell movements. *Physical Biology*. 2023;20:016002.
   [doi:10.1088/1478-3975/acb18c](https://doi.org/10.1088/1478-3975/acb18c).
6. Potdar AA, Lu J, Jeon J, Weaver AM, Cummings PT. Bimodal analysis of
   mammary epithelial cell migration in two dimensions. *Annals of Biomedical
   Engineering*. 2009;37:230–245.
   [doi:10.1007/s10439-008-9592-y](https://doi.org/10.1007/s10439-008-9592-y).
7. Patel D, Lin R, Majumder B, Ganusov VV. Brain-localized CD4 and CD8 T cells
   perform correlated random walks and not Lévy walks. *F1000Research*.
   2023;12:87.
   [doi:10.12688/f1000research.129923.2](https://doi.org/10.12688/f1000research.129923.2).
