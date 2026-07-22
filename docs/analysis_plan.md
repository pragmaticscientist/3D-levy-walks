# LTDB Lévy-walk surrogate analysis pipeline (`ltdb_levy`)

## Context

`fast_simulations` holds a C Monte-Carlo simulator for Lévy-walk search on a 3-D torus, driven by `.conf` files and shell scripts. It answers a *theoretical* question: how detection time scales with the Lévy exponent μ across target shapes.

We now want an *empirical* question answered: do real leukocyte trajectories (LTDB, 45 ground-truth files in `datasets/GT_TRACKS/`) behave like Lévy walks, and does a fitted exponent — or μ=2 specifically — reproduce their search performance on shape-controlled targets? That needs machinery the repo lacks: LTDB parsing, directed-run segmentation, exponent fitting with goodness-of-fit and model comparison, nonparametric relocation-length pools, a ladder of surrogate/control trajectory classes, and a finite-budget censored replay.

The power-law hypothesis is to be **tested, not assumed**, and nothing is tuned toward μ≈2.

**All new code is Python.** The only C involvement is a run of the *existing, unmodified* simulator at the fitted μ̂, which you will run. No new C is written and `func.c` / `simulation.c` are not edited.

---

## Decisions taken

| Topic | Decision |
|---|---|
| Language | **All new code in Python.** C used only for your own run at μ̂ with the existing binary. |
| Relocation definition | **Directed runs (Analysis B) primary**; frame displacements secondary. |
| Fitted model | **Your C mixture**, not a pure power law. **μ and one crossover per condition are fitted**; `lmax` is the empirical condition maximum. |
| Detection radius | **d = 1**, fixed. |
| Disk membership | **C's cylinder test**, not the exact rounded neighbourhood. |
| Line orientation | Along **y**, as in C. |
| Area matching | **Face-on maximum projection** — C's existing convention. |
| Direction sampling | Isotropic (`cos φ ~ U[−1,1]`) — **you already fixed this in `func.c`**. |
| Deterministic seeding | Not a priority; one master seed, kept simple. |
| Extra control | **Per-vector random rotation** of empirical displacement vectors. |
| Coordinate units | x, y, z **already in µm**; time = `t·dt` seconds. |
| Conditions | `(cell_type, organ, stimulus, dt)` from the paper's tables; fallback per `(video, population)`. |
| z-slab | Implement 3-D as specified; quantify slab anisotropy in diagnostics + limitations. |

---

## The model being fitted

Your C sampler ([func.c:119-131](func.c#L119-L131)) draws a Uniform[0,1) step with probability `a`, else a power-law step on `[1, lmax]`. Working out the implied density, the two branches join exactly:

> **f(ℓ) = a(μ, lmax) · min(1, ℓ^−μ)**, for ℓ ∈ [0, lmax]

flat below 1, power law above, continuous at the kink, where `a` is precisely the existing `get_normalization_constant` ([func.c:45-51](func.c#L45-L51)) — because `1/a = 1 + ∫₁^lmax ℓ^−μ dℓ`. Verified independently: substituting `u→1` into the tail inverse-CDF returns exactly `lmax`.

For a fixed crossover this gives the conditional log likelihood:

> **log L(μ) = n·log a(μ, lmax) − μ · Σ_{i : ℓᵢ ≥ 1} log ℓᵢ**

with the `μ=1` branch (`1/a = 1 + ln lmax`) handled analytically and tested at `μ = 1 ± 1e-9`.

The primary fit now maximizes this likelihood jointly over \((\mu,c)\), with
one crossover per eligible biological condition. Physical lengths remain
stored in µm; **3.75 µm**, half of a 7.5 µm reference leukocyte diameter, is
only the optimizer start and an external biological reference. `lmax` remains
the empirical condition maximum. The pipeline also reports the original fixed-
crossover sensitivity curve so the strong \(\hat\mu\)-versus-\(c\) dependence
stays visible. The `lmin` grid-selection machinery from the original spec is
dropped: this model has no separate lower cutoff.

Still performed: parametric-bootstrap KS goodness-of-fit (refitting μ each replicate), clustered bootstrap CIs, and model comparison against truncated exponential / lognormal / Weibull / power-law-with-cutoff over the same support via **track-held-out grouped cross-validation**. Where the mixture is not clearly supported, report an **"effective exponent"**.

---

## Data format — resolved from the source paper

Pizzagalli et al., *Scientific Data* 5:180129 (2018), Table 8:

- **Line 1** — `VideoID;dx;dy;dz;dt`; `dx,dy` voxel size (µm), `dz` z-spacing (µm), `dt` frame interval (**seconds**).
- **Line 2** — `ch0..ch4` channel-visibility flags. **Skip** (identifiable only positionally; the value varies by file).
- **Line 3+** — `TrackID;x;y;z;t`, positional, no column names.
- **x, y, z are already in micrometres** — "position with respect to the top-up-left most corner of the z-stack". **Do not scale by dx/dy/dz.** `t` is a frame index; time = `t·dt`.
- **`_a`/`_b` are distinct tracked cell populations**, not channels.

Measured corpus: **44 files (excluding the `SQUARE` fixture), 44,722 rows, 728 tracks**; track length min 2, median 45, max 241; 66 tracks under 10 points.

**Anomalies — audited, never silently fixed:** `LTDB004_b_GT.csv` declares `LTDB005_b` internally → prefer the filename, record `id_mismatch`; `SQUARE_GT..csv` excluded via config (and reused as a parser fixture); `CS003_a` has an all-zero channel row.

**Biological metadata** comes from the paper: Table 3 (`VideoID, CH0..CH3`) → cell type per population via `a:`/`b:` prefixes and stain notation; Table 4 (`VideoID, Site of Imaging, Immune stimulus, …`) → organ and stimulus. Cell types and sites are cross-checked against the deposited SQL database. Joined on VideoID + population, with `dt` from each file's own header. Written to `configs/ltdb_metadata.csv` with `provenance` and `verified` columns; **unverified rows are not pooled** without an explicit override, and **nothing is ever fabricated**.

---

## Conventions inherited from C

Reproduced exactly so your μ̂ run stays comparable. Recorded in `outputs/report/conventions.md`.

- Step law: the mixture above — **not** a pure truncated power law.
- Directions: isotropic, per your fix at [func.c:299](func.c#L299) and [func.c:215](func.c#L215).
- Detection radius **d = 1**; ball `‖q‖ ≤ D/2 + 1`.
- Disk: **cylinder** `ρ ≤ D/2 + 1` **and** `|Δz| ≤ 1` ([func.c:379-389](func.c#L379-L389)).
- Line: capsule along **y**, length `D`, radius 1 ([func.c:327-342](func.c#L327-L342)).
- Torus minimum-image ([func.c:53-73](func.c#L53-L73)); **endpoint-only detection**; target at box centre; single walker from a uniform random start (`nest` with `n_walkers=1` is identical to "uniform in torus").
- Area matching: face-on maximum projection, C's existing formulas ([func.c:21-43](func.c#L21-L43)).

**Two documented caveats** (recorded, not silently corrected):

1. **Face-on matching under random rotation.** A ball projects πR² from every direction; a randomly-oriented flat disk projects on average *half* its face-on area. Matching face-on therefore leaves the disk effectively smaller on average, so the shape-sensitivity score G partly reflects the matching convention. Stated plainly in `limitations.md`.
2. **No finite budget in C.** The C engine runs to detection; the primary analysis here is finite-budget and censored. These are different estimands (see Validation).

---

## Architecture

```
src/ltdb_levy/
  config.py  errors.py  cli.py            # YAML -> frozen dataclasses; argparse CLI
  dataio/   ltdb_reader.py schemas.py artifacts.py
  preprocess.py                           # sort/dedup/gap-split/displacements/turning angles
  segment.py                              # directed-run segmentation (primary analysis)
  pools.py                                # relocation-length pools + audit
  dist/     mixture.py families.py mle.py gof.py compare.py
  bootstrap.py  directions.py  classes.py  targets.py
  replay/   engine.py api.py              # vectorized reference engine (budgeted + unbounded)
  outcomes.py  stats.py  figures.py  report.py
configs/  primary.yaml  ltdb_metadata.csv
tests/
```

**Environment constraints (measured):** Python **3.9.2** → `from __future__ import annotations` in every module, no PEP-604 unions at runtime. **typer is not installed → use `argparse`** (zero new dependencies). numpy 1.26.4, scipy, pandas, matplotlib, pyyaml 5.3.1 (`safe_load` only) all present.

Cross-cutting: type annotations, docstrings, no global state, RNGs passed explicitly, structured logging, and **every exclusion counted into an audit table**.

### Replay engine — the key structural point

Because detection is endpoint-only, stateless, and the budget is a pure prefix-sum condition, a whole trial is branch-free:

```
vectors → lengths → cumsum → mask = cumsum ≤ B
        → endpoints = start + cumsum(rotated vectors), wrapped
        → hit = contains(endpoints[mask]) → first = argmax(hit)
```

No `while` loop anywhere. This is what makes pure Python fast enough to be the production engine — estimated **~60 s single-core** for the full grid (6 classes × 3 shapes × 5 sizes × 2 budgets × 10⁴ trials ≈ 1.8×10⁶ trials). The C engine's `while(1)` ([func.c:257](func.c#L257)) exists only for multi-walker asynchronous first passage, irrelevant at `n_walkers=1`.

Use `np.cumsum` (left-to-right), never `np.sum` (pairwise), for path lengths. Random rotations via **Shoemake's quaternion method** (three uniforms, reproducible).

**Common random numbers:** reuse the same initial positions and rotations across trajectory classes at a given trial index. Class contrasts become paired and their variance drops substantially — free, and this study's entire output is between-class differences.

---

## Work plan

**Phase 1 — Data spine.** `config.py`, `dataio/`, `preprocess.py`, artifacts + manifest; `inspect` and `preprocess` commands. Sort by time, resolve duplicates, **split fragments at any frame gap > 1** (never bridge a gap into one long relocation), compute displacements, speeds, and turning angles via `atan2(‖u×v‖, u·v)` — never `arccos`, which loses precision exactly at the thresholds that matter. Deliverable: tidy tracks table, per-file inspection report, anomalies surfaced. Build synthetic fixtures here.

**Phase 2 — Segmentation, pools, directions.** Directed runs end on: frame gap, fragment end, speed below threshold, or turning angle above threshold. **Deterministic turning-point convention:** the displacement whose turning angle exceeds the threshold *starts* the new run; ties break to the lowest index; documented and tested. `run_length = ‖p_b − p_a‖` (end-to-end) is the relocation length, with path length stored alongside so straightness is available. First and last run per fragment excluded from the primary fit (field-of-view censoring), included as sensitivity.

The primary turning-angle threshold is **90 degrees**, anchored to
Georgantzoglou *et al.* (2022, doi:10.1083/jcb.202103207), who classified
90--180-degree changes in leukocyte movement as large U-turns. Their
45-degree boundary between persistent steps and small turns, and the original
60-degree intermediate rule, remain prespecified sensitivities. This is a
literature anchor, not a claim that 90 degrees is a universal biological
constant: their angle construction is related but not identical to ours, and
immune-cell turning angles depend on the sampling interval (Beltman *et al.*,
2009, doi:10.1038/nri2638; Ganusov *et al.*, 2023,
doi:10.1088/1478-3975/acb18c). Conditions with different frame intervals
therefore remain separate.

The primary fit has **no positive speed cutoff**. Speed-aware sensitivities
treat steps below 1, 2, or 4 µm/min as pauses that terminate a run. The central
2-µm/min threshold is widely used to calculate the in-vivo lymphocyte arrest
coefficient (for example, Patel *et al.*, 2023,
doi:10.12688/f1000research.129923.2). This use requires a clear qualification:
the literature threshold ordinarily summarizes the proportion of time spent
below 2 µm/min; converting every such step into a hard run boundary is an
operational sensitivity analysis, not an independently validated segmentation
method. The 1- and 4-µm/min settings bracket the conventional value.

Pools preserve the **per-track nested structure** (needed for cell-balanced sampling, clustered bootstrap, auditability). Cell-balanced weights `w_i = 1/(n_tracks · n_{track(i)})` primary; step-weighted `1/n` as sensitivity. Full audit table + pool diagnostics (contributing tracks/runs, quantiles, per-track weight share, ECDF/CCDF, between-track and between-video heterogeneity, effective sample size).

**The segmentation sensitivity sweep is a first-class deliverable of this phase, not an appendix** — report μ̂, GoF p, and headline G across a grid of turning-angle and speed thresholds. Thresholds are *never* chosen by looking at the resulting μ. If conclusions are not stable across the grid, that is itself the finding.

**Phase 3 — Fitting.** `dist/mixture.py` (the model above), `families.py`
(alternatives), `mle.py`, `gof.py`, `compare.py`, `bootstrap.py`. One jointly
estimated exponent and crossover **per condition**, not per cell. Deliverables:
\(\hat\mu\) and \(\hat c\) with whole-track clustered-bootstrap intervals, KS
bootstrap p-value with both parameters refitted **reported alongside a power
curve** (a non-rejection at low power is how this literature goes wrong),
model comparison with both parameters refitted inside track-held-out grouped
CV, and the fixed-crossover \(\hat\mu\)-versus-\(c\) sensitivity curve.

**Phase 4 — Trajectory classes + replay + outcomes.** The control ladder, each rung removing exactly one structure:

| Class | Lengths | Directions | Order |
|---|---|---|---|
| `empirical` | real | real | fully preserved |
| `block_resampled` | real blocks | real within block | short-range only |
| `ordered_length_rotated` | real, in order | **per-vector random rotation** | length order kept |
| `pool_isotropic` | pool draw | isotropic | destroyed |
| `fitted_mu` | mixture at μ̂ | isotropic | destroyed |
| `mu2` | mixture at μ=2 | isotropic | destroyed |

Adjacent-rung differences are the interpretable quantities: *empirical vs
ordered-length* isolates directional correlation; *ordered-length vs pool* is
a composite replacement of the focal track's realized length sequence,
boundary runs, and fixed endpoint count by condition-level complete-run draws
with replacement; *pool vs fitted-μ̂* tests whether the parametric fit
reproduces real search behaviour; *fitted-μ̂ vs μ=2* tests whether μ=2 is
special. The requested per-vector rotation control appears twice — on the
ordered per-track sequence (`ordered_length_rotated`) and on pooled vectors
(`pool_isotropic`) — but these controls do not use the same length source.

Budget `B_i` = the source track's own cumulative relocation length; **overshoot discarded** by default (truncation manufactures a spike of short final steps that is larger for heavy-tailed classes, biasing toward the study's own hypothesis) with `truncate` as a labelled sensitivity. Never a `recycle` policy. Detection checked **only at relocation endpoints**, not at the initial position.

Outcomes: detection probability; **restricted mean detection distance** `E[min(τ,B)]`; conditional mean among detected, explicitly labelled and never used alone. Shape sensitivity `G = max_shape RMDD / min_shape RMDD` plus log form, and contrasts like `log RMDD(ball) − log RMDD(line)`. The independent unit is **track/video — not a Monte-Carlo replay**; all CIs route through the clustered bootstrap over 728 tracks (naive i.i.d. CIs would be far too narrow). Monte-Carlo error reported separately from biological uncertainty.

**Phase 5 — Figures, report, manifest.** All requested figures; `report.md`, `methods.md`, `limitations.md`, `analysis_manifest.json` (source hashes, config hash, git commit, versions, seeds, outputs).

`limitations.md` must state: sampling interval affects observed step lengths; segmentation choices affect run definitions; **the dataset contains no targets and no detection events**; an effective exponent does not prove an exact Lévy walk; deviations may reflect persistence, pauses, correlations or tissue structure; pools cannot exceed observed support; pooling can hide heterogeneity; **the z-slab makes empirical directions non-isotropic by construction**; **face-on area matching leaves a randomly-oriented disk effectively smaller than the ball**; and each condition-specific crossover is correlated with its fitted exponent and can be weakly identified in small samples.

**Phase 6 — Your C run at μ̂.** Emit `experiments/configs/detection_time_fitted_mu.conf` with `rangemu_LevyDistrib=μ̂`, the data-derived `l_max`, `list_shapes`, `surface_selector` and target sizes, plus the matching `.sh`. You run it with the existing binary; the pipeline ingests the resulting CSV for comparison. No C source is modified.

---

## Validation

Replacing the parity-suite apparatus, since no new C exists:

1. **Analytic single-step check** — the strongest test available. With one relocation from a uniform start, translation invariance on the torus gives `P(detect) = Vol(d-neighbourhood)/side³` in closed form. Assert the Monte-Carlo estimate matches within 4 binomial standard errors **for all three shapes**. This validates wrapping, membership, volume formulas and sampling simultaneously against a *known* answer — catching the "engine agrees with itself and both are wrong" failure that cross-checking cannot.
2. **Python-vs-C cross-check** — run the Python engine in **unbounded mode** reproducing the C experiment (C membership, d=1, line along y, face-on matching, matched μ and lmax) and compare mean detection time against your C run statistically. Validates the engine against trusted code with no new C.
3. **MLE recovery** — simulate ≥20,000 draws from the mixture at known μ and recover it within a statistically justified tolerance; stability at μ≈1.

Other tests: semicolon parsing (plus `SQUARE_GT..csv` as a real fixture, and `LTDB004_b` for the ID-mismatch path); splitting at missing frames; turning-angle correctness; segmentation on a path with known turns and pauses **including the tie-breaking convention**; first/last-run exclusion; inverse-CDF sampling; clustered bootstrap preserving whole tracks; isotropy of direction sampling; minimum-image wrapping at adversarial inputs (exact multiples of `side`, tiny negatives); ball/disk/line membership including boundary points; endpoint-only detection; budget matching and overshoot handling; cell-balanced giving equal per-track selection probability; step-weighted giving equal per-run probability; pool support containment; ordered-length preservation; block resampling (within-block vectors intact, one rotation per block); and a small end-to-end run on synthetic tracks.

---

## Verification

```bash
python -m ltdb_levy inspect    --config configs/primary.yaml   # format report + anomalies
python -m ltdb_levy preprocess --config configs/primary.yaml
python -m ltdb_levy fit        --config configs/primary.yaml   # mu_hat, CI, GoF, model comparison
python -m ltdb_levy replay     --config configs/primary.yaml --pilot
python -m ltdb_levy replay     --config configs/primary.yaml
python -m ltdb_levy report     --config configs/primary.yaml
pytest -q
```

Acceptance: `pytest` green; the analytic single-step check passes for all three shapes; MLE recovers known μ; the pilot returns detection probabilities away from 0 and 1 **before** the confirmatory `matched_areas` are frozen (and not changed afterwards); Python unbounded mode agrees statistically with your C run at μ̂.

---

## Risks

1. **Segmentation thresholds drive everything.** The relocation-length distribution *is* a function of the angle and speed thresholds — a power law can be created or destroyed by moving them. Mitigated by making the sensitivity sweep a Phase-2 deliverable.
2. **Crossover and exponent can be weakly identified jointly.** One crossover
   is fitted per condition, with 3.75 µm used only as the optimizer start.
   Whole-track joint bootstrap intervals, boundary flags, held-out fitting, and
   the fixed-crossover sensitivity curve expose the resulting uncertainty.
3. **Pool minimums may exclude much of the corpus.** With 15-track/100-run thresholds and 66 tracks under 10 points, expect several conditions to be reported ineligible — reported, never silently dropped.
4. **The mixture may simply not fit.** That is a legitimate result; the report says so plainly and switches to "effective exponent" language.
5. **Quasi-2D confinement** confounds 3-D isotropy comparisons — diagnosed and documented per the decision taken.
6. **Face-on matching is not orientation-fair** (caveat 1 above) — documented, so G is not over-read.
7. **Metadata verification is a human checkpoint** before any published claim.
8. **Never fabricate.** Any analysis that cannot be completed is reported as such rather than filled with plausible numbers.
