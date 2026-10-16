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
report implements four of them, declines one, and in doing so overturns part of
the headline result it was building on.

The headline is this: **`whamp_track_v2` is the first of the four solvers in
this tree that reports no branch hops *and* no unjudgeable points.** Over the
same nine adversarial sweep lines, measured against the same step-refined
reference at one-eighth the sweep spacing, `whamp.m` records 13 branch hops,
Version 1 in both its builds records 0 hops but 19 points past where the
reference could follow them, and v2 records **0 hops and 0 points
past the end** — it reports exactly the branch that exists, and nothing else.
The cost is 3.7 s against 3.0 s for the same work.

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
Item 3 (back-port to the Fortran and C engines) was **declined by decision, not
omission**: v2 stays a MATLAB-only reference.

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
| 3 | Back-port to `src/whamp_engine.f90` and the C engines | **Declined** (owner decision) | §7 |
| 4 | Test `varyKzFirst = 1` and genuine 2-D sweeps | **Done**; found and fixed a real defect | §5 |
| 5 | `flagSuspectBranch`-aware plotting helper | **Done** | §6 |

Two files were added, both untracked at the time of writing:

* [`srcm/whamp_track_v2.m`](srcm/whamp_track_v2.m) — 1995 lines, generated from
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
a small circle centred on it — 33 sample points, no residual threshold anywhere.
Because `D` is meromorphic, the count is `Z − P`:

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

## 7. Item 3: the Fortran and C back-port was declined

The back-port to `src/whamp_engine.f90` (`root_finding`, lines 182–241) and
`srccopt64/whamp_engine.c` (`ddDX` seed at line 355) is **not done, by
decision**. v2 remains a MATLAB-only reference implementation.

This is worth stating explicitly rather than leaving as a silent omission,
because it means the defect described in [`REPORTtrack.md`](REPORTtrack.md) —
a single carried seed with no branch-identity test — is still live in the
engines that most users actually run. The MATLAB reference now demonstrates the
fix works without a hand-tuned constant; porting it is a separate, well-specified
piece of work with an existence proof behind it.

## 8. Results

Nine adversarial sweep lines — three multi-species plasmas (golden, anisotropic,
beam-plasma) at `kperp` = 0.005, 0.01, 0.02 — solved at spacing `h = 0.008` and
judged against a step-refined reference at `h/8` with `refineFactor = 8`,
`fstart = 0.4`. The reference could judge **435 points** in total.

| solver | hops | in range | past end | time |
|---|---|---|---|---|
| `whamp` | **13** | 42 | 19 | 0.72 s |
| `whamp_track_test` (modular v1) | 0 | 56 | 19 | 7.72 s |
| `whamp_track` (self-contained v1) | 0 | 56 | 19 | 2.92 s |
| **`whamp_track_v2`** | **0** | **57** | **0** | 3.65 s |

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

**Cost.** v2 costs about 25% more than the self-contained v1 (3.65 s against
2.92 s for the nine lines). The difference is the certification pass: a 33-point
contour evaluation of `D` per accepted root. It buys a per-point guarantee
rather than a statistical one, which seems the right trade, and `certify` can be
switched off.

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

## 11. Reproducing

The harnesses live in the session workspace alongside the v2 head file and the
assembler; they are not part of the tracked tree.

```
% Build v2 from its head file plus the inlined +whampkit closure.
% -KitExtra must be QUOTED or PowerShell splits it on the comma.
.\assemble.ps1 -Main v2_head.m -MainDir <session files dir> `
    -Out whamp_track_v2.m `
    -KitExtra "evalD,windingNumber" `
    -FuncName whamp_track_v2

smoke_v2        % v2 vs v1 bit-identity on the regression line
eval_final      % the headline table: hops / in range / past end
eval_2d         % item 4: varyKzFirst 0 and 1, cross-ordering, 2-D vs 1-D
diag_cross      % full root scan per grid point (proved 0.48871 is not a root)
diag_shape      % fine walk that found the fold
verify_fold     % how many points each solver reports past the fold
diag_tol        % tolerance actually used, and certification statistics
```

## 12. Conclusion

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
And that anything is fixed for users who run the Fortran or C engine, because
the back-port was declined and those engines still carry the original defect.

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
4. Revisit the Fortran and C back-port with v2 as the reference. The MATLAB
   version now shows the fix needs no tuned constants, which removes the main
   objection to carrying it in a compiled engine.
5. Sweep the option space (`trustFactor`, `jumpFactor`, `certRadius`,
   `certSamples`) to confirm the defaults are not load-bearing, since the point
   of this work was to remove exactly that kind of dependence.
