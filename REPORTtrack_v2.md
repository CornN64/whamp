# Branch Tracking Without Hand-Tuned Scales (`srcm/whamp_track_v2.m`)

AI-assisted work by Walter Puccio 2026, extending
[`REPORTtrack.md`](REPORTtrack.md) (which documents the original tracking
layer, the `+whampkit` extraction, and the measurement harness this report
reuses).

> **Naming.** The three solvers were renamed after this report was first
> written, so that the file you are meant to run has the plain name. Where
> [`REPORTtrack.md`](REPORTtrack.md) says `whamp_track.m` it now means
> `whamp_track_test.m`:
>
> | now | was | what it is |
> |---|---|---|
> | `whamp_track_test.m` | `whamp_track.m` | modular Version 1, calls `+whampkit`; kept as the reference build |
> | `whamp_track.m` | `whamp_track_standalone.m` | Version 1, self-contained |
> | `whamp_track_v2.m` | `whamp_track_standalone_v2.m` | Version 2, no hand-tuned scales |
>
> This report uses the new names throughout.

## Abstract

[`REPORTtrack.md`](REPORTtrack.md) ended with five suggested next steps. This
report implements four of them and initially declines one; the C back-port has
since been implemented in `srccopt64`, while the Fortran back-port remains
deferred.

The headline is this: **`whamp_track_v2` is the first of the four solvers in
this tree that reports no branch hops *and* no unjudgeable points.** Over the
same nine adversarial sweep lines, measured against the same step-refined
reference at one-eighth the sweep spacing, `whamp.m` records 13 branch hops,
Version 1 in both its builds records 0 hops but 19 points past where the
reference could follow them, and v2 records **0 hops and 0 points
past the end** — it reports exactly the branch that exists, and nothing else.
It does this at **2.83 s against 3.02 s** for the self-contained v1 on that
benchmark, and at **1.14× v1** on the large two-dimensional whistler sweep that
first exposed a real cost — the guarantee is now cheap, but §11 records that it
took three separate fixes to get there and that the first two of them did not
help the case that mattered.

The work was not a clean implementation. Three of the five items produced a
result that contradicted the obvious interpretation, and each one is the useful
part:

* **Sizing a continuity tolerance from a prediction is circular.** Doing so —
  the natural reading of item 1 — reproduced the baseline's 13 hops exactly.
  The tolerance absorbs precisely the error it was meant to detect.
* **`+whampkit/dispersion.m` reports `DP = 0` and `DZ = 0` on the axis planes**,
  for its own numerical reasons. Every sweep line starts on one of those planes,
  where the branch genuinely moves by a factor of 2.4. A zero derivative is the
  *absence* of a prediction, not a prediction of no movement, and treating it as
  the latter throws away the first point of every line.
* **The branch in the headline example ends in a fold at `kpar ≈ 0.0425`.**
  The reference stops there, and so does the exact implicit step. The extra
  coverage that made v1 look more thorough past that point is on a *different
  branch*, reached because the hand-tuned `relFloor` was generous enough to
  place a lucky seed. Removing `relFloor` does not lose coverage; it declines
  to claim coverage that was never there.

Item 2 (winding-number certification) removed the lower-bound caveat from the
hop count. Item 4 (two-dimensional sweeps) initially broke v2 badly — the two
sweep orderings disagreed by 76 — and fixing it required replacing the
visit-order history with a per-axis local rate map, after which the orderings
agree to 3·10⁻¹¹. Item 5 (a confidence-aware plotting helper) is delivered.
The Fortran portion of item 3 remains deferred; the all-double C engine now
has a v2 implementation in `srccopt64`.

`whamp.m` and both Version 1 builds are **unmodified in behaviour** — the only
edits to them were the renames in the table above. v2 reproduces v1's output
bit-for-bit on the regression line (`max |v2 − v1| = 0.000e+00`) while differing
from it wherever v1 was wrong.

## 1. Scope: what was done, and what was not

The five items from §10 of [`REPORTtrack.md`](REPORTtrack.md), with the outcome
of each:

| # | Suggested step | Status | Where |
|---|---|---|---|
| 1 | Replace the median step-scale and `relFloor` with the exact implicit-differentiation step from `DP`/`DZ`/`DX` | **Done**, but not as suggested — see §2 | §2 |
| 2 | Winding-number ground truth, to remove the lower-bound caveat | **Done** | §3 |
| 3 | Back-port to `src/whamp_engine.f90` and the C engines | **C done; Fortran deferred** | §7 |
| 4 | Test `varyKzFirst = 1` and genuine 2-D sweeps | **Done**; found and fixed a real defect | §5 |
| 5 | `flagSuspectBranch`-aware plotting helper | **Done** | §6 |

Two files were added, both untracked at the time of writing:

* [`srcm/whamp_track_v2.m`](srcm/whamp_track_v2.m) — 2314 lines, generated from
  a source-of-truth head file plus the inlined `+whampkit` dependency closure,
  in the same style as [`srcm/whamp_track.m`](srcm/whamp_track.m).
* [`srcm/whamp_track_plot.m`](srcm/whamp_track_plot.m) — 253 lines, item 5.

## 2. Item 1: removing the last hand-tuned constant

### 2.1 What v1 did, and why it was unsatisfactory

Version 1 decided whether an accepted root was still on the branch it was
following by comparing it against a scale built from two quantities:

* the **median of the branch's recent step sizes**, and
* **`relFloor = 0.25 · |x|`**, a floor on that scale.

The report's own complaint is correct: `relFloor` is a magic number with no
derivation. What was not appreciated until this work is that it was doing *two*
jobs, and removing it breaks both.

### 2.2 The implicit step

The dispersion relation `D(x; p, z) = 0` defines `x` implicitly as a function of
the wave numbers. Its exact first-order movement is

$$
\Delta x \;=\; -\,\frac{(p - p_0)\,D_P + (z - z_0)\,D_Z}{D_x}
$$

evaluated at the branch's own converged root. `+whampkit/dispersion.m` already
returns `DP`, `DZ` and `DX` for the group velocity, so this costs nothing. v2
uses it as the prediction:

```matlab
xPred = a.x - ((p - a.p) * DP + (z - a.z) * DZ) / DX;
```

`relFloor` is deleted. The only floor remaining is `sqrt(eps) * |x|`, which is
not a tuning knob but the statement that a step below the root's own round-off
cannot be resolved.

§10 of [`REPORTtrack.md`](REPORTtrack.md) measured v1's complaint on the flat
anisotropic branch at `h = 0.004`: the branch stepped ~0.004, `scale` floored at
0.247, and the resulting jump tolerance `6 × scale` was **1.48 — about 60× larger
than it should be**. Measuring v2 on that same line and spacing gives a median
`scale` of **0.0078**, i.e. a median jump tolerance of **0.047**, roughly 30×
tighter, against a largest genuine step of 0.34. The tolerance now tracks the
branch instead of flooring against it.

### 2.3 The circularity trap

The obvious way to finish item 1 was to size the *rejection* threshold from the
prediction's own recent error. This was implemented, measured, and **made things
dramatically worse: 13 hops, identical to the unfixed engine.**

The reason is structural and worth stating plainly, because it is the same trap
that §4.1 of [`REPORTtrack.md`](REPORTtrack.md) describes from the other
direction. A tolerance derived from the prediction is an allowance for the
prediction being wrong. If the prediction is wrong because the branch is about
to be lost, the tolerance widens by exactly enough to forgive losing it. The
auditor and the audited are the same number.

This is why v2 does **not** size the tolerance from the prediction. It sizes it
from the *step the branch actually took* (§2.5), which is an observation about
the branch and not a restatement of the thing being tested.

### 2.4 `dispersion.m` zeroes the derivatives on the axis planes

The second surprise. `+whampkit/dispersion.m` contains, in its own finite
differencing:

```matlab
if ZZ(1) ~= 0     % ... compute DZ
end
if PP(1) ~= 0     % ... compute DP
end
```

The guards avoid a `1/k` division in the differencing stencil. They are correct
for their purpose. But **every sweep line in this test suite starts at
`kpar = 0`**, and at `kpar = 0` the implicit step evaluates to `xPred = a.x`
exactly — a confident prediction that the branch has not moved at all, when in
fact it moves from `x = 0.408` to `x = 0.981` in a single step on the golden
model at `kperp = 0.005`.

A zero derivative is the **absence** of a prediction, not a prediction of zero
movement. v2's `local_predict` therefore returns a second value:

```matlab
function [xPred, informative] = local_predict(...)
```

and the reference falls back to the grid anchor where the derivatives are not
informative. Suppressing the prediction entirely at `kpar = 0` was also tried
and also failed — it let v2 accept the engine's hop to `x = 131.6` on line 1,
because with no prediction there is nothing for the continuity check to
disagree with. The anchor is the right fallback precisely because it is a
*weaker* claim, not a stronger one.

### 2.5 The yardstick: a rate, not a residual

The remaining question is what to compare against. The natural candidate — the
model residual `|root − prediction|` — is the wrong one, and
`diag_flags.m` shows why: on a completely clean `kperp` line with no hop
anywhere, that residual grows from 3·10⁻⁴ to 0.11 purely because of local
curvature. It is a property of how strongly `D` bends, not of whether the branch
is intact. Using it as the tolerance produced **7 false suspect flags on a line
v1 solves perfectly**, and taking a running maximum instead of a median made it
worse, not better.

A **step** is the right yardstick because a branch hop lands a branch-spacing
away, whereas curvature keeps the residual small compared with one step. v2
therefore measures the rate `|dx| / |dk|` at which the branch moves through the
grid, and sets the tolerance from the rate local to the cell being entered.

## 3. Item 2: certification by the argument principle

Every accepted root in v2 is certified with a winding-number count of `D` around
a small circle centred on it — an adaptive polygon starting at 8 samples and no
residual threshold anywhere (§11.2). Because `D` is meromorphic, the count is
`Z − P`:

| count | meaning |
|---|---|
| `+1` | a simple zero, and no pole inside — the reported root is real |
| `0` | neither zero nor pole inside — the reported point is **not a root at all** |
| `−1` | one pole and no zero — the point sits on a species gyroharmonic |

This works precisely where a residual test cannot. `|D|` ranges over roughly
thirty orders of magnitude across this domain, so no threshold on it is
meaningful; the winding number is a property of `D` alone and is scale-free.

On the regression line, **99 of 99 accepted roots certify as `+1`**, with zero
`flagNotCertified` and zero `flagCertUnknown`. On the anisotropic line of §8,
all six accepted roots certify `+1`.

The consequence for the headline claim matters: in
[`REPORTtrack.md`](REPORTtrack.md) the hop count was a **lower bound**, limited
by how far the step-refined reference itself converged. Certification is not a
lower bound — it is a direct statement about each individual root, independent
of any reference. The caveat is gone from the per-point claim, and what remains
bounded is only the *comparison* between solvers.

Three separate output flags keep "wrong branch" and "not a root" from collapsing
into one another: `flagSuspectBranch`, `flagNotCertified`, `flagCertUnknown`.

## 4. The fold, and why v1's extra coverage was not coverage

This is the finding that reframes the whole comparison, and it came from walking
the branch at 32× the sweep resolution (`h = 0.00025`) rather than from any
solver change.

On the anisotropic plasma at `kperp = 0.005`, the tracked branch does not
continue indefinitely in `kpar`. It **terminates in a fold at `kpar ≈ 0.0425`**,
where `dx/dkpar → 0` and `x` plateaus at `0.970741`. Past the fold there is no
root connected to this branch; the next root Newton finds is at `x = 1.382233`,
a different mode.

This has three consequences.

**The reference stops at the fold too.** The step-refined reference is itself a
continuation, so it dies at the same place. Points beyond it are not "harder
points the reference failed on" — they are points where refinement *cannot
judge at all*. Counting them as ordinary misses was the error §8.1 of
[`REPORTtrack.md`](REPORTtrack.md) identified and this work had not fully
corrected.

**v1's extra points are on a different branch.** `verify_fold.m` confirms that
v1 reports 4 points past the fold on this line. Since the branch does not exist
there, those 4 points are necessarily on another branch, reached because
`relFloor = 0.25·|x|` kept the tolerance wide enough that a scan-offset seed
happened to land in the right basin. v1's own continuity check cannot see this,
because the check is calibrated by a floor that permits the jump that causes it.

**`pastEnd` is the metric that exposes it.** The harness in §8 therefore reports
`pastEnd` — points a solver reports beyond where the reference gave out — as a
first-class number alongside hops. It is the only place a silent branch switch
can hide, and it is the number v2 takes to zero.

## 5. Item 4: two-dimensional sweeps and `varyKzFirst`

This item was listed as a test, and it found a real defect.

### 5.1 How v1's history fails in 2-D

Version 1 keeps its step history in **visit order**. On a single-`kperp`
`kpar` line that is exactly right. On a genuine 2-D sweep it is not: at the
start of a new line, the previously visited point sits at the *far end* of the
previous line, at a completely different wave number. The first point of every
line after the first is therefore judged against a step taken between two
unrelated points.

v2 replaces the ordered history with a map of accepted roots indexed by grid
cell, and takes each point's predecessor from its inner neighbour, falling back
to the same inner index on the previous outer line. On 1-D lines this reduces
exactly to v1's behaviour.

### 5.2 The failure, and the two bugs behind it

With grid anchoring alone, `varyKzFirst = 0` produced first-column roots

```
0.409   0.878   12.55   26.58   45.82   69.96
```

with 17 suspect flags, and the two sweep orderings disagreed by **7.59·10¹**.
Two separate causes:

**Pooling steps across axes.** The history mixed steps taken along `kpar` with
steps taken along `kperp`. On this model the branch moves roughly 70× faster in
`kpar` than in `kperp`, so the pooled scale was meaningless on one axis or the
other. Fixed by making the history a **rate per axis**, `|dx| / |dk|`, stored
separately for each axis. The `dk` in that ratio is
`abs(p − p₀) + abs(z − z₀)`, which needs no knowledge of which axis is varying
— a predecessor always differs from its successor in exactly one coordinate.
(Assuming the inner axis was always `kperp` fixed one ordering and broke the
other; this is why the symmetric form matters.)

**A windowed median is not a local rate.** After the axis fix, one point still
came out as `0.4887` where the true root is `0.9046`. `diag_cross.m` — a full
root scan at every grid point — proved that `0.48871` **is not a root at all**:
Newton started from that seed converges to `0.90457`. The tolerance was simply
too wide, because a median taken over a *visit-order window* was drawn from the
flat tail of the previous line, where `dx/dkpar ≈ 0.06`, while the point being
entered sat near `kpar = 0`, where `dx/dkpar ≈ 57`. A window in visit order is
not a neighbourhood in the grid.

The fix is a **local rate map**: `rateMap(ip, iz, axis)` stores the rate at the
grid cell each step *arrived at*, and `local_rateNear` reads the neighbouring
cells by expanding Manhattan rings, taking the maximum of the ties. The
tolerance is then set by the part of the branch the solver is entering, which is
the only thing that could be relevant.

### 5.3 Result

| measure | v1 | v2 |
|---|---|---|
| `varyKzFirst = 0`, first-column roots | 0.4092 0.8781 0.9531 0.9669 0.9700 0.9707 | identical |
| `varyKzFirst = 1`, first-column roots | 0.4092 0.8781 0.9531 0.9669 0.9700 0.9707 | identical |
| cross-ordering deviation (24 common points) | — | **3.19·10⁻¹¹** |
| 2-D sweep vs single-line runs | — | **exactly 0** |
| suspect flags, `varyKzFirst = 0` | 1 | 1 |
| suspect flags, `varyKzFirst = 1` | 4 | 4 |

Both orderings now agree with each other to round-off and with v1 exactly. The
visit-order history assumption in v1 is confirmed unsound in 2-D, and v2 does
not depend on it.

## 6. Item 5: `whamp_track_plot.m`

[`srcm/whamp_track_plot.m`](srcm/whamp_track_plot.m) draws a dispersion diagram
on which branch confidence is visible directly:

| mark | meaning |
|---|---|
| solid line + filled dot | found and certified |
| open circle | `flagNotCertified` — reported point is not a root |
| grey square | `flagCertUnknown` — certification not available |
| red ✗ | `flagSuspectBranch` — continuity check failed |
| dashed overlay | the implicit first-order prediction, against the accepted root |

Gaps are never bridged: marker plots carry `'LineStyle','none'`, so a missing
point is drawn as a gap rather than as a straight line across it.

Options: `Mode` (`'confidence'` / `'plain'`), `ShowPredicted`, `Frequency`,
`YScale`, `Transpose`, `MarkerSize`, `LineWidth`. The return value is a struct of
`gobjects` (`lines`, `suspect`, `notCertified`, `certUnknown`, `predicted`) for
later annotation.

Two orientation facts cost more time than the plotting itself and are worth
recording: **WHAMP_TRACK returns arrays as (kpar × kperp)**, so a dispersion
curve at fixed `kperp` is a *column* of the raw array, and `Transpose,true`
transposes. The helper's internal field accessor applies the same transpose as
`f`, or the confidence flags silently describe the wrong points. The TeX
interpreter also rejects `k_{\parallel}`; the axis label uses `k_{par}`.

Verified headlessly (no display): confidence mode gives 4 curves, 4 suspect
marks, 0 not-certified, 0 unknown, 4 prediction overlays; `Transpose` gives 6
curves; `Mode','plain'` gives 0 suspect marks; v1's output plots too, with
absent fields handled. A rendered PNG confirms the red ✗ marks and the dashed
prediction diverging exactly where first-order continuation fails.

## 7. Item 3: C back-port implemented; Fortran deferred

The v2 tracker is implemented for the all-double C engine in
[`srccopt64/whamp_track.c`](srccopt64/whamp_track.c), with its public interface
in [`srccopt64/whamp_track.h`](srccopt64/whamp_track.h). Callers configure the
existing `comin` globals as for `whamp_engine()`, optionally call
`track_options_default()` to reset options, then call
`whamp_track_v2(&diagnostics)`. It writes the existing `comoutput` arrays and
the additional prediction, continuity, and winding-number arrays declared in
[`srccopt64/comtrack.h`](srccopt64/comtrack.h). The independent
`whamp_track_test` target checks both 2-D sweep orders, winding certification,
and exact agreement with the legacy engine's first root. The C source was
compiled and this test was run with GCC.

The Fortran back-port to `src/whamp_engine.f90` remains deferred. The original
engine entry points are unchanged; callers opt into the tracker by calling its
separate API. This C implementation removes the MATLAB-only limitation for
`srccopt64`, but does not yet change the Fortran engine or the other C variants.

### 7.1 Numerical agreement with the MATLAB reference

Agreement is checked by [`srcm/whamp_track_v2_vs_c.m`](srcm/whamp_track_v2_vs_c.m),
which runs `whamp_track_v2` and the C tracker over twelve sweeps and compares
every output array. The C side is
[`srccopt64/whamp_track_dump.c`](srccopt64/whamp_track_dump.c)
(`mingw32-make whamp_track_dump`), which prints each case's arrays at `%.17g`.
The cases cover both example plasmas, one- and two-dimensional sweeps in both
line orders, negative and logarithmic axis steps, an isotropic model, a
three-species beam, and a start frequency above the branch.

| case | points | max rel | verdict |
|---|---|---|---|
| `line_kpar`, `line_kpar_pfirst` | 7 | 1.9e-12 | MATCH |
| `line_kperp` | 6 | 2.8e-08 | close |
| `grid2d_zfirst`, `grid2d_pfirst` | 20 | 5.4e-09 | close |
| `grid2d_negstep` | 20 | 6.1e-13 | MATCH |
| `grid2d_log` | 6 | 8.2e-08 | close |
| `ex2_grid`, `start_high` | 9, 7 | 0 | MATCH |
| `ex2_log` | 15 | 1.4e-16 | MATCH |
| `isotropic` | 5 | 6.8e-13 | MATCH |
| `three_species` | 20 | 1.7e-09 | close |

**Every flag matches exactly in every case** — `solutionFound`,
`tooHeavilyDamped`, `noConvergence`, `suspectBranch`, `notCertified`,
`certUnknown` — and so does `certWinding`. The two implementations therefore
make identical decisions throughout the sweep; only the last digits of the
roots differ, at 1e-14 to 1e-10, which is the Newton iteration's own stopping
noise (`FREQ_REL_TOL = 1e-6` per relative correction). The single-point engine
comparison, `whamp_vs_srccopt64`, gives a floor of 3.2e-14 for the same reason.

Two of the comparisons need a yardstick other than the array's own magnitude.
The Poynting flux is a difference of nearly cancelling products, and `SGP`/`SGZ`
are `imag(f)/VG` where both are zero to machine precision; there a per-cell
ratio reads a difference of 1e-17 as a relative error of one. Those quantities
are instead judged by whether the observed gap could have come from the gaps in
the inputs they are computed from — `|dS| <= |dE||B| + |E||dB| + 2*eps*|E||B|`,
and likewise for the quotient. The bound holds everywhere (largest ratio 0.4),
so the difference is inherited from the fields and not a logic divergence. The
six field components share one yardstick, the largest field magnitude in the
sweep, because on the `kperp = 0` line `Ez` and `Sx` vanish by symmetry.

Two divergences found by this comparison are fixed in `whamp_track.c`:
`tr_predecessor` now walks back at most 64 cells along the outer axis, as
`local_predecessor` does, and `tr_winding` returns the half-integer-resolved
winding on an early resolved contour instead of rounding it to an integer, as
`kit_windingNumber` does. Neither is reachable on the twelve cases, which is
why both were invisible to the flag comparison; they matter on grids wider than
65 cells and for contours whose phase lands on a half-integer.

## 8. Results

Nine adversarial sweep lines — three multi-species plasmas (golden, anisotropic,
beam-plasma) at `kperp` = 0.005, 0.01, 0.02 — solved at spacing `h = 0.008` and
judged against a step-refined reference at `h/8` with `refineFactor = 8`,
`fstart = 0.4`. The reference could judge **435 points** in total.

| solver | hops | in range | past end | time |
|---|---|---|---|---|
| `whamp` | **13** | 42 | 19 | 0.75 s |
| `whamp_track_test` (modular v1) | 0 | 56 | 19 | 7.74 s |
| `whamp_track` (self-contained v1) | 0 | 56 | 19 | 3.02 s |
| **`whamp_track_v2`** | **0** | **57** | **0** | **2.83 s** |

*`hops` = roots disagreeing with the refined reference within the range the
reference reached. `in range` = reference-converged points the solver also got.
`past end` = points reported beyond where the reference gave out — roots
refinement cannot check, and the only place a silent branch switch can hide.*

**v2 is the only solver in this tree with 0 hops and 0 past-end.** The three
track-based solvers all report 0 hops; only v2 refuses to report anything past
the fold, which is the correct answer, and it gains a point of genuine coverage
in exchange (57 against 56).

Per-line, the two numbers that move are lines 1 and 3 of the golden model —
where `whamp` hops once and reports 7 points past the end, v2 reports 0 and 0 —
and lines 4–6 of the anisotropic model, where `whamp` hops four times on each
line and reports 4 points past the end on each, and v2 reports 0 and 0 on all
three. The beam-plasma lines (7–9) are clean for every solver.

**Cost.** v2 is marginally *cheaper* than the self-contained v1 on these nine
lines (2.83 s against 3.02 s). It was not — an earlier version of this table
showed 3.65 s, 25% more than v1, and the first real-world sweep was reported as
about ten times slower than v1. That was a defect, not the price of the method;
§11 gives the three causes, the fixes, and the proofs that none of them changed
a reported root. These nine lines are small and converged, so they flatter v2 —
on the large two-dimensional sweep that actually prompted the complaint v2 is
1.14× v1, which is the number to quote.

## 9. Verification

| Check | Result |
|---|---|
| `whamp.m` modified? | **No.** Byte-identical to its pre-work state |
| Version 1 behaviour (`whamp_track_test.m`, `whamp_track.m`) | **Unchanged** — renamed only; both still produce identical output |
| v2 vs v1 on the regression line | `max \|v2 − v1\| = 0.000e+00`, 99 points, 0 suspect |
| v2 certification on the regression line | 99/99 certified `+1`, 0 not-certified, 0 unknown |
| v2 certification on the anisotropic line | 6/6 certified `+1` |
| `eval_final` re-run after every change | stable at 0 hops / 57 in range / 0 past end |
| v2 tolerance on the flat branch (`h = 0.004`) | median `scale` 0.0078, 11/11 certified, 0 suspect |
| `varyKzFirst = 0` vs `= 1`, v2 | agree to **3.19·10⁻¹¹** over 24 common points |
| 2-D sweep vs single-line runs, v2 | deviation **0.00e+00** |
| First-column roots, both orderings, v1 and v2 | identical to 4 decimals |
| Fold at `kpar ≈ 0.0425` | confirmed by an independent walk at `h = 0.00025` and by `verify_fold.m` |
| Plot modes | all exercised headlessly; object counts as in §6; PNG visually confirmed |
| `relFloor` | **deleted** from v2; only `sqrt(eps)·\|x\|` remains |
| Rate-map lookup vs original unbounded ring walk | `verify_ring.m`, 8000 randomised lookups, grids to 40×40, sparse-biased fill, both axes: **0 mismatches** |
| Order-1 vs order-2 `D` inside certification | 40 random states: **bit-identical in `D` and `IERR`**, 0 mismatches |
| Adaptive contour vs fixed-32 reference | `verify_wind.m` 267 cases: **0 mismatches**; `verify_certify.m` start-resolution sweep: 0 mismatches at the shipped start of 8 (§11.2) |
| `certify = false` vs `true` | `max\|v2 − v2_nocert\| = 0` — certification changes no reported root |
| Batched `D` vs scalar seam | `verify_batch.m`: 12 564 values of `D` and `IERR` over randomised single- and multi-species models, **bit-identical, 0 mismatches**; 240 windings vs the scalar contour, 0 mismatches |
| Whistler sweep, v2 after §11.4 | `smoke_v2` `max\|v2 − v1\| = 0`, 99/99 certified; `eval_final` stable at 0 hops / 57 in range / 0 past end; `eval_2d` cross-ordering 3.19·10⁻¹¹ |

## 10. Known limitations

**Four points are still flagged suspect under `varyKzFirst = 1`.** This is the
same count as v1, so v2 is not worse — but it is not better either, and the
difference is that v2 *shows* them instead of silently accepting them. Removing
them entirely would need the branch rate stored per grid cell **before** the
step is taken, rather than recorded from the step that has just succeeded. The
current scheme bootstraps from the previous line, which is one step behind.

**The rate is recorded after the step, not predicted before it.** This is the
same root cause as the point above, and it is the obvious remaining weakness. A
rate field built by a cheap preliminary pass over the grid would close it.

**The reference still cannot judge past the fold.** Certification removes the
lower-bound caveat for *individual roots*, but the *comparison* between solvers
still rests on the refined reference, and that reference is a continuation and
stops where the branch stops. `pastEnd` is the honest reporting of this, not a
solution to it. A fold-aware reference — one that follows the branch through the
turning point by continuing in `(x, k)` rather than `x(k)` — would be needed to
say anything about the region past `kpar ≈ 0.0425`.

**Certification is local.** A winding number of `+1` proves the reported point
is a simple root of `D`. It does not prove it is the *same* root as at the
neighbouring point — that is what continuity is for. The two tests are
complementary and neither subsumes the other.

**The `dispersion.m` derivative guards remain.** v2 works around them by
treating a zero derivative as no prediction. The cleaner fix is for
`+whampkit/dispersion.m` to report the derivative on the axis planes by a
one-sided stencil, or to return a flag alongside it. That is a change to shared
physics code and was deliberately not made here.

**Only three plasmas and three `kperp` values were tested.** Zero hops is not
claimed to hold generally.

**The 10× slowdown was not reproduced, only bounded.** The reported case was
eventually obtained and v2 measured against v1 on it: 2.0× before the work of
§11.4, 1.14× after. That is a real regression that has been largely closed, but
it is not a factor of ten, and no input in this tree produces a factor of ten.
If the original observation came from a build older than the §11.2 adaptive
contour, or from a sweep that hit the degenerate-axis hang of §11.5, that would
explain the discrepancy; neither can be confirmed from here.

**Certification is still a tenth of the runtime, and it is optional.** On the
whistler sweep the default costs 1.11× what `certify = false` costs. The default
stays `true` because the guarantee is the point of v2, but a user producing a
figure rather than a claim can remove that overhead by turning it off, and
`max|v2 − v2_nocert| = 0` says no reported root changes.

**A degenerate axis still hangs the sweep.** §11.5 documents a pre-existing
non-termination in the shared sweep logic for a `[v 0 v]` axis. It affects
`whamp.m` (which hangs with no message at all) as much as v2, and it is left
unfixed so that the reference stays byte-identical.

## 11. Performance: v2 was too slow in practice, and why

The first real-world use of v2 reported it running roughly **ten times slower**
than v1 on the same sweep. That is a defect, not a cost of the method.

To be exact about what was and was not reproduced: on the small, fully-converged
examples in this report v2 was only ever about **2×** v1, so the headline 10×
could not be reproduced there at all. Profiling a larger, sparser run found two
pathologies whose cost depends on the *shape* of the sweep rather than its
length, which is the right signature for a discrepancy of that size. Both were
fixed without changing a single reported root, and both fixes were then *proved*
equivalent to what they replaced.

Neither of them helped the real case. The reported sweep was eventually obtained
and measured (§11.4), and there the cost was not shape-dependent at all — it was
the plain price of evaluating `D` nine times per accepted root. That is the third
cause, and the one that mattered. An unrelated defect surfaced while hunting for
these and is recorded in §11.5.

### 11.1 Cause 1: an expanding-ring search that was only linear when the map was full

Item 4 replaced the visit-order history with a per-axis **local rate map**
(§5.3): each grid cell that converges records the rate `dω/dk` it achieved, and
a later point seeds itself from the rates in nearby cells. The lookup walked
outward ring by ring — radius 0, then 1, then 2 — until it found a cell holding
a rate for the requested axis.

That loop was bounded by `2·nMax`, i.e. by the *whole grid*, and it terminates
early only if it finds something. On a grid where an axis has **no recorded
rates at all** — which is exactly what happens at the leading edge of a sweep,
in a region the solver has not yet converged into, or on a wide declared range
where most points fail — every lookup exhausts every ring. The cost is then
`O(N²)` per lookup in the grid size, and it is paid at every point. Measured on
a 201×201 map, per lookup:

| cells holding rates | original walk | hybrid lookup |
|---|---|---|
| 100% | 0.525 ms | 0.500 ms |
| 10% | 0.268 ms | 0.118 ms |
| 1% | 0.036 ms | 0.0025 ms |
| 0.1% | 0.245 ms | 0.0020 ms |
| **0%** | **0.926 ms** | **0.157 ms** |

The pathological column is the *empty* one: 0.926 ms to discover that there is
nothing there, repeatedly, at every point. Note also that the cost is not
monotone in the fill — the 1%-filled map is 7× *cheaper* than the 0.1%-filled
one, because at 1% a hit is usually found within a few rings, whereas at 0.1%
most queries travel a long way or give up entirely. Cost tracks the distance to
the nearest hit, which is exactly the quantity a real sweep does not control.

**The fix** keeps a list of the cells that actually hold rates and makes the
lookup hybrid: ring out to a small cap of 8 (cheap, and the common case), and if
that finds nothing, fall back to a single vectorised minimum-distance scan over
the occupied cells. The fallback is a *complete* search, so the answer cannot
depend on where the cap sits — the cap is a fast path, not an approximation.

This was verified rather than argued: `verify_ring.m` runs 8000 randomised
lookups over grids up to 40×40, with fill deliberately biased sparse and both
axes independently occupied, against the original unbounded walk. **0
mismatches.**

### 11.2 Cause 2: the certification contour was fixed at 33 points

Item 2 certifies each accepted root by counting the winding of `D` around a
small closed contour (§3). The contour was a fixed 33 samples of order-2 `D`.
That is 24.6 ms per accepted root, and it is paid once per root that succeeds —
so it scales with how well the solver is doing, which is an awkward incentive.

Two changes:

* **Order 1 instead of order 2.** `dispersion.m` accumulates the `DX` block only
  when `order ≥ 2`, and certification needs the value of `D`, not its
  derivatives. Order 1 was checked against order 2 at 40 random states: **bit
  identical in both `D` and `IERR`, 0 mismatches**, for 1.43× on the
  certification call. (Order 0 is *not* usable — `E(:,1:0)` is empty, the
  species susceptibility is never accumulated, and `D` comes out as garbage such
  as `-7.78e+30` where the truth is `3.49e+27`. Order 1 is the minimum correct
  order, and it is what is used.)
* **An adaptive contour.** Start at 8 samples; double while the largest phase
  step around the contour exceeds `π/2`, up to three doublings. The ceiling is
  therefore 32 samples — the old fixed cost — and a contour enclosing one simple
  zero has near-linear phase, so **98% of certifications resolve on the first
  8-sample pass**. Mean evaluations per certification fell from 32.63 to 9.12.

The starting resolution matters and was chosen by measurement, not taste.
`verify_certify.m` compares against a non-refining fixed-32 reference (the
original shipped behaviour) over three plasmas × accepted-root, jittered and
deliberately-far probes × three contour radii:

| start | mismatches vs fixed-32 | mean evaluations |
|---|---|---|
| 4 | 2 | 5.4 |
| 6 | 2 | 6.8 |
| **8 (shipped)** | **0** | **9.12** |
| 12 | 0 | 12.9 |
| 16 | 0 | 16.6 |
| 32 | 0 | 32.63 |

Starts of 4 and 6 occasionally under-resolve and disagree with the reference, so
the default is 8. `verify_wind.m` repeats the comparison over 267 cases at the
shipped setting: **0 mismatches.**

### 11.3 Result on the tuning suite — and why it did not settle the matter

On the 404-point sweep used for tuning, v2 went from **4.17 s to 2.29 s**, and
the certification share of runtime from **58% to 26%**.

The headline nine-line benchmark of §8, re-run after the fixes, now reads **2.83 s
for v2 against 3.02 s for the self-contained v1** — v2 is faster than the solver
it used to be 25% slower than, with the same 0 hops / 57 in range / 0 past end as
before the fixes. On two-dimensional sweeps (`bench_real`) v2 runs at
**0.60–0.94× the wall time of v1** across grids from 286 to 2626 points — it is
now the faster solver there too, which matters because 2-D is where the rate map
lives.

That looked like the end of the story. It was not, and the reason matters more
than the numbers above. When the *actual* case that produced the complaint was
obtained — a 100 × 100 two-dimensional whistler sweep over an anisotropic
electron population — neither fix made any measurable difference. Both §11.1 and
§11.2 address costs that this sweep does not incur: the rate map fills quickly
so the ring search always hits inside a few rings, and the contour was already
resolving on its first 8-sample pass. The tuning suite is small, converged and
one-dimensional at heart; it could not see the real cost, which was simply the
price of certification itself.

### 11.4 Cause 3: certification is nine full evaluations of D per root

The user's case is `Whistlereanisotropy.m`: a single electron species,
`Te = 100 eV`, `ne = 1 cm⁻³`, `B = 50 nT`, `a = 5`, on a 100 × 100 grid of
`kperp` and `kpar` from 0.01 to 1. Measured before any of the work in this
section, with three solvers on the identical input:

| solver | wall time |
|---|---|
| `whamp_track` (v1) | 27.07 s |
| `whamp_track_v2` after §11.1–11.2 | 54.71 s |

So v2 was **2.0× v1** and **2.4× the reference engine** — still not the reported
10×, but a real and reproducible regression. A run with `certify` disabled
isolated the cause immediately:

| configuration | 30 × 30 grid |
|---|---|
| v1 | 4.46 s |
| v2, `certify = true` | 8.30 s |
| v2, `certify = false` | 3.69 s |

Certification was **56% of v2**, and v2 without it was *faster than v1*. The
tracking machinery is not the cost; the guarantee is.

Why it costs what it does was established by direct measurement rather than
inference. Over 876 real certified roots the contour **resolved at 8 samples
100% of the time** — never once refined — with a worst phase step of 0.8222
against the `π/2` threshold. So the adaptive loop of §11.2 is already optimal
here and there is nothing left to win there. The cost is irreducibly nine scalar
evaluations of `D` per root, each of which rebuilds the dispersion state and
runs the eight-pole susceptibility sum from scratch.

**The fix was to evaluate the contour in one vectorised pass.** A contour is
eight or sixteen frequencies at a *single* shared `(p, z)`, so the state build
and the pole loop can be shared across all of them. `whampkit.evalDbatch` is
that vector mirror of `evalD(...,1)`, and `windingNumber` now builds the whole
vertex row and calls it once. Two things were needed to make it actually pay:

* **Freezing converged Taylor samples by masking rather than by indexing.** The
  first batched version carried a `live` mask and updated `rc(1,live)` and
  friends. On a nine-element row that is eight logical-indexed reads and writes
  per pass, and it made the batch *slower than the scalar loop it replaced* —
  0.97×, measured over all 876 real roots. Writing the freeze as a zeroed term
  with a substituted harmless denominator is bit-identical (adding zero leaves a
  partial sum untouched, and the scalar `break` adds nothing further either) and
  costs four elementwise products. That alone took the batch from 0.97× to 2.52×.
* **Stacking all eight poles into one call.** `rBatch` is elementwise and every
  pole shares the same `ala`, so the eight arguments can be laid out pole-major
  in one 72-element row and summed in a single call instead of eight 9-element
  ones. Grouping cannot change the result, because each sample's series is
  advanced independently including its own convergence test. This brought the
  single-certification microbenchmark to **9.3×** and the real-root measurement
  to **2.96×**.

`pole_constants` was also memoised with a `persistent` pair — same values, built
once instead of per call — which removes a fixed allocation from every
susceptibility evaluation in the toolkit.

The whole batched path is **bit-identical** to the scalar seam it replaces, and
this is asserted rather than assumed: `verify_batch.m` compares 12 564 values of
`D` and `IERR` across randomised single- and multi-species models, drifting and
loss-cone species included, and 240 windings against the scalar contour. **0
mismatches.**

The result on the user's own case, three interleaved repetitions each so that
machine drift is spread rather than charged to whichever solver runs last:

| solver | best | median |
|---|---|---|
| `whamp` | 18.03 s | 18.20 s |
| `whamp_track` (v1) | 21.42 s | 21.52 s |
| v2, `certify = false` | 21.97 s | 22.27 s |
| **v2, `certify = true`** | **24.49 s** | **24.63 s** |

v2 went from **2.0× v1 to 1.14× v1**, and from 2.4× to **1.36×** the reference
engine. Certification fell from 56% of runtime to about **10%**, and from 4.62 s
to 1.08 s in the profiler. The guarantee is now cheap enough to leave on by
default.

The `certify` option remains the blunt instrument if it is ever wanted: switching
it off removes the pass entirely, and `max|v2 − v2_nocert| = 0` confirms it
changes no reported root. `certSamples` sets the contour's *starting*
resolution, so `struct('certSamples', 32)` restores the original fixed-32
behaviour exactly.

An earlier, non-interleaved run of that option sweep appeared to show 16 samples
being *cheaper* than 8, which would have implied the adaptive loop was somehow
costing more than the evaluations it saved. It was not: re-measured with the
configurations interleaved over four repetitions, the ordering is monotone as it
should be — best 2.65 s at 8 samples, 2.71 s at 16, 3.02 s at 32 — and the roots
are identical throughout. The apparent inversion was machine drift charged to
whichever configuration happened to run first. The shipped default of 8 is
correct, and `bench_csamp.m` is the check.

**The general lesson survived being wrong.** §11.1 and §11.2 were real fixes to
real pathologies, but they were found by profiling the wrong workload, and the
report nearly claimed the case closed on the strength of them. Both were loops
whose average cost depends on how much data has been collected — genuinely
pathological at a sparse frontier, and genuinely irrelevant to a dense whistler
sweep where the cost was the raw number of `D` evaluations. The methodological
point stands but is weaker than it looked: a benchmark on a small converged
example cannot see a shape-dependent bug, *and* fixing the shape-dependent bugs
still does not license a claim about total cost. The workload that prompted the
complaint has to be the one that is measured.

### 11.5 A separate, pre-existing defect found on the way

Chasing the 10× led to sweeps that did not finish at all, and one of those is a
**genuine bug that predates v2 and is shared by `whamp.m`, both Version 1 builds
and v2.** It is recorded here rather than fixed, because fixing it changes the
behaviour of the reference engine, which §9 has kept byte-identical throughout.

`whampkit.axis_spec` accepts a degenerate axis `[v 0 v]` — start equal to end
with a zero step — and documents that this is legal:

```matlab
whampkit.axis_spec([0.005 0 0.005], 'kperp')   % -> [0.0050 0.0050 0.0000]
```

`fill_axis` handles it correctly and returns a single-point axis. But
`advance_sweep` walks the outer axis by *repeated addition of the step*, and with
a step of zero the position never changes, so the termination test
`PLG >= PM1 && PLG <= PM2` stays true forever:

```
start  PLG=0.0050  ZLG=0.0000  done=0
step 1  PLG=0.0050  ZLG=0.0000  done=0
step 2  PLG=0.0050  ZLG=0.0000  done=0
   ...  done never becomes true
```

The sweep therefore never sets `done`. `whamp.m` has no backstop at all, so it
loops forever with no message. The three track solvers carry a
`guard > 200000` backstop, so after tens of minutes they raise `Runaway` — "The
sweep did not terminate" — which at least stops, but reports a symptom rather
than the cause. A single `kperp` value written as `[0.005 0 0.005]` instead of
the scalar `0.005` is enough to trigger it, which is an easy thing to type.
`verify_zerosweep.m` reproduces it in three lines, without running a solve.

The natural fix is for `advance_sweep` to treat a zero step as a one-point axis
and set `done` after the first visit, matching what `fill_axis` already assumes.
That is a one-line change to shared code and it is deliberately **not** made
here: it alters `whamp.m`, and every claim in §9 rests on `whamp.m` being
untouched. It should be a separate change with the engine's own tests.

## 12. Reproducing

The harnesses live in the session workspace alongside the v2 head file and the
assembler; they are not part of the tracked tree.

```
% Build v2 from its head file plus the inlined +whampkit closure.
% -KitExtra must be QUOTED or PowerShell splits it on the comma.
.\assemble.ps1 -Main v2_head.m -MainDir <session files dir> `
    -Out whamp_track_v2.m `
    -KitExtra "evalD,evalDbatch,windingNumber" `
    -FuncName whamp_track_v2

smoke_v2        % v2 vs v1 bit-identity on the regression line
eval_final      % the headline table: hops / in range / past end
eval_2d         % item 4: varyKzFirst 0 and 1, cross-ordering, 2-D vs 1-D
diag_cross      % full root scan per grid point (proved 0.48871 is not a root)
diag_shape      % fine walk that found the fold
verify_fold     % how many points each solver reports past the fold
diag_tol        % tolerance actually used, and certification statistics

verify_ring     % 8000 randomised rate lookups vs the original unbounded walk
verify_wind     % adaptive contour vs fixed-32, 267 cases
verify_certify  % certification agreement by contour start resolution (table in 11.2)
bench_speed     % three-way timing on a 404-point sweep, with certify on/off
bench_real      % v1 vs v2 on two-dimensional sweeps of varying size
bench_scale     % 1-D scaling, 1x100 to 1x1600
bench_ring      % rate-map lookup cost against map fill fraction
bench_big       % large sparse grids, new v2 vs the original ring walk
bench_sparse2   % wide 1-D sweeps, new v2 vs the original ring walk
verify_zerosweep  % reproduces the non-terminating [v 0 v] axis of 11.5

verify_batch      % evalDbatch vs the scalar seam: 12564 values, 0 mismatches
bench_real_contour % contour cost at 876 REAL certified roots, scalar vs batched
bench_stable      % interleaved 100x100 whistler timing: whamp / v1 / v2 (table 11.4)
bench_w_cert      % 30x30 whistler, certify on / off / certSamples 4 and 16
bench_csamp       % certSamples 8/16/32 interleaved (guards the drift artifact)
prof_v2           % MATLAB profiler breakdown of v2 on the whistler case
```

`bench_real_contour` is the one to trust over any microbenchmark: it measures the
contour at the roots a real sweep actually produces, and it is what caught the
batched path being *slower* than the scalar loop it replaced.

## 13. Conclusion

**The last hand-tuned constant is gone, and removing it changed the result
rather than merely tidying it.** `relFloor = 0.25·|x|` was not a harmless
allowance. It was wide enough to let a scan-offset seed land in a neighbouring
basin past a fold, and v1's continuity check could not object because the floor
was what calibrated the check. Replacing it with the exact implicit step did not
lose coverage — it declined a claim that was never earned, and the honest
coverage is one point *higher* than before.

**Four of the five suggested steps are done, and the one that was listed as a
test found the real bug.** Two-dimensional sweeps broke v2 badly — orderings
disagreeing by 76 — and the cause was that a window in visit order is not a
neighbourhood in the grid. The per-axis local rate map that replaced it brings
the two orderings into agreement at 3·10⁻¹¹ and leaves the one-dimensional
results bit-identical. The visit-order history assumption in v1 is unsound in
2-D; it survives only because v1 was tested in 1-D.

**The guarantee is no longer expensive, but it was never free.** The first
real-world run was reported at roughly ten times v1's wall time. Three causes
were found. Two were shape-dependent loops — an expanding-ring search over the
rate map costing `O(N²)` per lookup whenever an axis held no rates, and a
certification contour fixed at 33 samples of order-2 `D` regardless of whether 9
sufficed — and both were replaced with versions *proved* equivalent to what they
shipped, 8000 randomised lookups and 267 contour cases, 0 mismatches. Neither
helped the reported case at all. The cost there was the raw price of
certification: nine scalar evaluations of `D` per accepted root, 56% of runtime.
Evaluating the contour as one vectorised pass, with converged Taylor samples
frozen by masking rather than indexed and all eight poles stacked into a single
call, took v2 from 2.0× v1 to **1.14× v1** and certification from 56% of runtime
to about 10%, bit-identically over 12 564 compared values. Hunting for these also
surfaced a genuine pre-existing bug, shared by `whamp.m` and both Version 1
builds, in which a degenerate `[v 0 v]` axis makes the sweep never terminate; it
is documented in §11.5 and deliberately not fixed, because fixing it would change
the reference engine.

**Three findings are more general than this solver.** Sizing a continuity
tolerance from a prediction is circular — the tolerance forgives exactly the
error it is meant to catch, and this was measured at 13 hops, identical to the
unfixed baseline. A derivative reported as zero for the internal convenience of
the code that computes it is not a prediction of zero, and every sweep line in
this suite starts where that happens. And the model residual `|root −
prediction|` measures local curvature, not branch integrity; using it as a
threshold flags clean lines. All three produce plausible output until the
measurement is checked against a direct run.

**What is claimed:** over nine sweep lines spanning three multi-species plasmas
and three perpendicular wavenumbers, `whamp_track_v2` records **0
branch hops and 0 points beyond the reference's reach**, where `whamp.m` records
13 and 19 and Version 1 records 0 and 19; every accepted root it reports is
independently certified as a simple root of `D` by the argument principle; it
agrees with itself across sweep orderings to 3·10⁻¹¹; and it contains no
hand-tuned scale. `whamp.m` is byte-identical to its pre-work state, and both
Version 1 builds are unchanged in behaviour.

**What is not claimed:** that zero hops holds generally — three plasmas and
three `kperp` values are not a proof. That the four suspect points under
`varyKzFirst = 1` are unavoidable; they are the visible residue of a rate that is
recorded after the step instead of predicted before it, and closing that is the
obvious next move. That certification establishes branch *identity* — it
establishes that each reported point is a root, and continuity does the rest.
And that anything is fixed for users who run the Fortran engine or the other C
variants: this back-port covers only the all-double `srccopt64` tree, and the
Fortran engine still carries the original defect.

**Suggested next steps**, in rough order of value:

1. Predict the branch rate per grid cell in a cheap preliminary pass, rather
   than recording it from the step that has already succeeded. This should
   remove the four remaining suspect points and makes the tolerance genuinely
   predictive.
2. Let `+whampkit/dispersion.m` report `DP` and `DZ` on the axis planes by a
   one-sided stencil, or return an explicit "not available" flag, so that
   callers do not have to infer it from a zero.
3. Make the reference fold-aware — continue in `(x, k)` through a turning point
   — so that the region past `kpar ≈ 0.0425` becomes judgeable and `pastEnd`
   can shrink instead of merely being reported.
4. Port v2 to the Fortran engine, using the MATLAB implementation and the
   all-double C tracker as references. The other C variants can be considered
   separately if they need the same tracking API.
5. Sweep the option space (`trustFactor`, `jumpFactor`, `certRadius`,
   `certSamples`) to confirm the defaults are not load-bearing, since the point
   of this work was to remove exactly that kind of dependence.
