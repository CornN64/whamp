# WHAMP: Translation, Porting and Branch Tracking — Consolidated Report

AI-assisted work by Walter Puccio, 2026.

This report merges and supersedes five earlier ones — which have been deleted.
It records six phases of work on the WHAMP plasma-wave code and the evidence for
each.

## Scope

| Phase | Deliverable | Verified by |
|---|---|---|
| 1 | `src/` Fortran 90 → `srcc/` C, byte-for-byte | 67 differential comparisons vs a gfortran oracle |
| 2 | `srccopt/` (tidied, byte-exact) and `srcc64/` (all-`double` experiment) | 30-session golden suites + oracle replay |
| 3 | `srcm/whamp.m` pure-MATLAB port, refactored | 30-case `%.20g` regression + conformance to `srcc64` |
| 4 | `whamp_track` v1 — branch tracking with a trust region | nine adversarial sweep lines vs a step-refined reference |
| 5 | `whamp_track` — tracking with no hand-tuned scales | same harness, plus winding-number certification |
| 6 | Tracker + `mexwrap` MEX gateway back-ported to C | `whamp_track_test`, `whamp_track_dump`, MATLAB suite |

Toolchain throughout: `gcc` / `gfortran` 6.3.0 (`C:\MinGW\bin`), MATLAB R2020b,
Windows 10/11 x64.

The method is the same in every phase: **anchor to the real reference program
before editing anything**, dump at full precision rather than comparing with a
tolerance, change one thing at a time, and re-diff after every change.

---

## 1. Fortran → C translation (`srcc/`)

WHAMP's interface is a line-oriented REPL whose exact column layout, separators,
blank records and help texts are part of its observable behaviour, so "the same
numbers" is not an acceptable acceptance criterion. Three rules made exact
comparison possible:

1. **One-to-one file mapping** — each `src/<name>.f90` becomes `srcc/<name>.c`,
   same routines in the same order; Fortran `MODULE`s become
   `comin`/`comcout`/`comoutput` `.h`/`.c` pairs.
2. **A byte-exact gfortran format emulator** — instead of trusting `printf`,
   `srcc/whamp.h` implements gfortran 6.3's field-generation rules for every
   edit descriptor WHAMP uses: `fE()` (`1P…Ew.d` with gfortran's exponent
   adjustment, round-half-even, `********` overflow), `fF()`, `fint()`, and
   `ld_real/ld_int/ld_logical/ld_complex()` for list-directed records including
   separator rules and `(re,im)` complex form. `F32()` reproduces the
   single/double distinction.
3. **Faithful control flow** — named `CYCLE`, labelled `GOTO`, `ENTRY`, `SAVE`
   locals, implied-DO lists and trailing-`/` `FORMAT` semantics translated
   explicitly rather than approximated.

Verification is differential: a reference binary is built with `gfortran -O2`
from the untouched `src/`, both binaries are driven with identical arguments,
working directory and piped stdin, and raw stdout bytes are compared. Because
WHAMP loops forever on end-of-file (identically in both), the harness reads a
bounded prefix and kills the process.

**Result: 64 of 67 comparisons byte-identical** — the `doc/test.txt` acceptance
session, all 19 output-format letters on two plasma models, `-debug` and
`-maxiterations` paths, parameter overrides, grid edge cases, and the callable
engine test.

### Bugs differential testing found (invisible to inspection)

* **List-directed separator after a colon-terminated character item.** Fortran
  suppresses the separator; the C port printed an extra space.
* **Trailing `/` in `FORMAT` emits an extra empty record.** Two places (`FORMAT
  131` in `typin`, `FORMAT 105` in `output`) were missing a newline.
* **`CYCLE output_format_loop` restarts the *outer* loop**, but had been
  rendered as `continue` on the inner one — so after printing output help the C
  program kept scanning the input string instead of re-prompting.

Plus format-routing corrections: the `x` (phase) case needs `F8.2` editing (which
prints ` NaN`, not MinGW's `1.#R`), and `p`/`z` need `F12.7`.

### Known divergences (all ≤ 1 ULP, or undefined behaviour in the original)

* **`Ex2` with `-debug`** — a 1-ULP seed in a case the solver never converges is
  amplified to ~2×10⁻⁶ in the *debug trace* only; all final and physical output
  is byte-identical. Isolated with two probes proving libm and C99 `_Complex`
  arithmetic are bit-identical between the two front ends; the residual comes
  from expression-level association, and is insensitive to `-O0…-O3` and
  `-ffp-contract`.
* **`R_negf`** — re-solve at `f = -0.1` diverges at iteration 4 where the
  condition number is ~4.5×10⁶; final `X` differs by 8×10⁻¹⁶, twelve orders
  inside the solver's own tolerance.
* **Engine tests ET6–ET8** — `whamp_engine.f90` never initialises grid entries
  the solver does not visit, so the Fortran prints denormal garbage that changes
  between runs *of the same binary*. The C uses `calloc`, so those entries are
  deterministically `0.0`.

Not translated, deliberately: `src/wrfi.f90` (dead — no `CALL WRFI` anywhere)
and `src/mexwhamp.F` (MATLAB gateway; superseded in phase 6).

---

## 2. Optimised variants (`srccopt/`, `srcc64/`)

* **`srccopt/`** — dead code removed, the Fortran edit-descriptor emulation
  extracted into `fedit.c/.h`, `whamp.h` slimmed, argument parsing and record
  reading rewritten with named helpers, and `ALLOC1`/`ALLOC2` collapsing 21
  hand-written `free()` calls on error paths into one unwinding sequence.
  **Byte-for-byte identical to `srcc/` on all 30 golden sessions.**
* **`srcc64/`** — the all-`double` experiment: every `F32(literal)` became a
  plain literal, `typin.c`'s four `float` variables (the only single-precision
  variables in the whole Fortran source) became `double`, and constants were
  collected in `constants.h`.

Both build warning-free at `-Wall -Wextra -O2` and at `-Wall -Wextra -Wpedantic`.
Under the broadest set gcc 6.3 implements, diagnostics fall from 238 (`srcc/`)
to 82 (`srccopt/`) to 72 (`srcc64/`); the residue is `-Wfloat-equal` (48,
inherent — Fortran compares reals with `==`) and `-Wdouble-promotion` (26,
deliberate). Silencing them would mean changing behaviour.

Against the authoritative oracle, `srcc/` and `srccopt/` differ on 3 of 25
sessions (the pre-existing ULP residuals above); `srcc64/` differs on 18,
which is the float→double experiment.

### 2.1 Single precision is load-bearing — a sweep silently loses its endpoint

For the sweep `p0,.05,.01z.0022f.1`:

| Build | Points emitted | Last p |
|---|---|---|
| `srcc/`, `srccopt/` | 5 | 0.0400000 |
| `srcc64/` | **6** | **0.0500000** |

The original **drops the requested endpoint**. The loop advances while
`PLG <= PM2`; with float rounding the accumulator reaches
`0.049999998882412910` against a bound of `0.049999997019767761` and fails the
test, while in double it lands exactly on the bound and succeeds. A sweep from 0
to 0.05 in steps of 0.01 *should* have six points, so `srcc64` is arguably
more correct — and further from the reference. **Widening `float` to `double`
changes program behaviour, not just its last digit.**

### 2.2 The Fortran source disagrees with itself about π and about *c*

`av.f90` uses a 14-digit π; `whamp_engine.f90`, `rasy.f90`, `rint.f90` and
`output.f90` use a 15-digit one. `av.f90` uses `299.79` while `energy.f90` uses
`299.792458`. `energy.c` scales the wave magnetic field by `1/c` and `av.c` then
forms `|c·B|²/|E|²`, so the two uses must agree; the disagreement left that ratio
wrong by 1.6×10⁻⁵.

The π difference is provably output-safe (the 14-digit value only reached
`POYN[]`, which the original computes and discards). The *c* difference is not:
`enden=` is `Re(ERG·conj(ERG)/(ERG+conj(ERG)))·2 − ENFIELD`, two nearly equal
terms, so a 1.6×10⁻⁵ change in `c²` becomes a **34 % change** in the printed
value. Neither printed value is trustworthy — the original already destroys 4–5
significant digits there. The exact SI `c` is used because it is *correct*, not
because it is *stable*; `enden=` should not be quoted to more than one
significant digit.

### 2.3 A real mistranslation: complex `COS`/`SIN`/`EXP` lost their imaginary part

`-Wconversion` reported five `discards imaginary component` warnings pointing at
three statements present identically in all trees: `COT = cos(PI*Y)/sin(PI*Y)`
in `rasy.c` and `rint.c`, and `Pv = exp(AL*C - Y*Xloc)` in `rint.c`. In every
case the argument is complex, so Fortran's elemental `COS`/`SIN`/`EXP` return a
complex value. C99 specifies type-generic macros that select the complex
overload; **MinGW-w64 gcc 6.3.0 does not implement that**, so the argument is
silently converted to `double`.

It survived the original byte-for-byte verification, a 25-session golden suite,
and review, because `RYLA` routes to `RASY`/`RINT` only for strongly damped
waves and every case in the battery used the weakly damped `p0z.0022f.1`.
Against the oracle on `heavy_damping` (`p2z.5f1.5`):

| Build | Result |
|---|---|
| Fortran oracle | no convergence, `I= 51 IRK= 1` |
| before the fix | `Local minima!`, converges at `I= 7` — a spurious root |
| after the fix | **byte-identical to the oracle** |

Fixed with explicit `ccos`/`csin`/`cexp` in all trees. `#include <tgmath.h>` was
also built and is byte-identical, but was rejected because it would silently
retype unrelated future calls. Every trig/exp/log/sqrt/pow call in all 18
translation units was then audited against its argument's declared type: these
three were the only ones receiving a complex argument.

### 2.4 Coverage gap

No golden session exercises `enden=` / `enfl_p=` / `enfl_z=`, and the
`all_outputs` case prints only help text because its output string is rejected
before use. The energy-output path in `av.c` is **untested in every tree**.

### 2.5 Recommendation

Use **`srccopt/`** where output must match the original WHAMP. Treat
**`srcc64/`** as an experiment — clean and self-consistent, but not closer to
reference behaviour. **Enable `-Wconversion` on any port to a new toolchain**: it
was the only warning class that caught a real defect, and `-Wall -Wextra` was
clean both before and after the fix.

---

## 3. Pure MATLAB port (`srcm/whamp.m`)

Refactored from 938 to 1129 lines (19 functions → 34); the increase is
documentation and extracted helpers, not new behaviour. `solve_point` was split
into `solve_point`/`dispersion_state`/`newton_root`/`evaluate_at` with named
tolerances; the 14 file-scope predictor variables of the C engine became named
fields of a `sweep` struct; the triplicated `DC` expression was collapsed;
residue tables and Gauss–Legendre abscissae were hoisted into helpers.

**Verification.** `whamp_dump.m` runs 30 cases and writes 21 fields at `%.20g`
(2616 lines of golden text), covering both sweep orders, negative and
logarithmic axes, three-species and loss-cone plasmas, ion drift, heavy damping,
non-convergence, `kperp = 0` and `kpar = 0`, three magnetic fields, and all 7
expected error identifiers. After 9 refactor batches the dump is **byte-identical
to the baseline**, and it stayed so after each individual batch.

Crucially the harness was anchored *outward* first: the pre-refactor code was
compared against the real `srcc64` transcripts and **9 of 9 comparable
sessions matched**, the closest at **3.2×10⁻¹⁴** relative error. The refactor
therefore preserves agreement with the reference program, not merely
self-agreement. `checkcode` went from 7 warnings to 0.

The conformance scripts caught a mistake *in the test*: a ~3×10⁻⁴ disagreement
that looked like a numerical defect was a units error — `Models/Ex1` specifies a
1 eV electron temperature, not `whamp.m`'s built-in 100 eV.

### Deliberate divergences from the C reference

All predate the refactor, are commented at the point of divergence, and are cases
where the C behaviour is a memory-lifetime artifact rather than intended physics:

1. **No cross-point state leakage.** The C engine keeps `D`, `DX`, `SG` in file
   scope, so a failed point continues with values left by an unrelated earlier
   point. `whamp.m` resets per point.
2. **Heavily-damped points return immediately** (status 2) rather than carrying
   on with a Newton step derived from stale globals.
3. **Grid indices are clamped** to `[1, count]`; the C writes out of bounds for
   negative steps.

---

## 4. Branch tracking, version 1 (`whamp_track`)

### The symptom and its cause

`whamp.m` sweeps `kpar` carrying one complex predictor. On a drift-unstable
plasma the reported frequency leaps between distant modes and then rides the new
branch — e.g. reporting 14.04, 28.10, 47.33, 71.42 where the true branch is flat
at ≈0.97. Nothing flags it: **every one of those values is a genuine,
well-converged root.**

The cause is not a numerical defect. Multiple roots genuinely coexist at a
single `(k⊥, k∥)`; `D(ω)` is meromorphic with one pole per species gyroharmonic;
so a Newton step started near a pole lands in whichever basin it falls into.
`roots()` needs a polynomial and `fzero` is real-only. **The fix must be a policy
layer on top of the existing iteration, not a better iteration.**

### The design

Physics was mechanically extracted into `srcm/+whampkit/` (36 files) so every
solver evaluates a bit-identical dispersion function; `whamp.m` was not modified
and its regression still passes byte-for-byte. Three mechanisms:

1. **Trust-region Newton** — retries confined to a disc centred on the *previous
   root on this branch*, radius `trustFactor ×` the branch's recent step size.
2. **Continuity check** — a converged root that moved more than `jumpFactor ×`
   the branch step is rejected and the point retried.
3. **Multi-start re-anchoring** — engine extrapolation, then previous root, then
   scaled offsets; first candidate that converges *and* passes continuity wins.

> **The decision that mattered: every scale is measured on the *branch*, never on
> the extrapolated seed.** A bad seed is bad by roughly the size of the jump it
> is about to cause, so any radius expressed in terms of the seed is circular —
> it admits precisely the failure it is meant to catch. This single change took
> the solver from ineffective to zero hops.

The **first seed is tried with `radius = Inf`** — it *is* `whamp`'s solve. On a
well-behaved line the machinery never engages.

### Bugs found

* **Trust-region Newton converging against the wall.** When the trust budget was
  exhausted `limitStep` returned `correction = 0`, which passes the step-size
  convergence test trivially — so the iteration "converged" while pinned against
  the disc, recording residuals of order 10⁸–10¹⁰ as roots. Recognised by a
  "root" appearing at *exactly* `trustFactor × scale`. Fixed by returning an
  `exhausted` flag and requiring `~wasClamped`. **Any trust region must be able
  to report that it bit.**
* **A bug in the measurement.** `judge()` skipped any point the reference had not
  converged to; on hard lines the reference dies early, so every baseline hop
  beyond that was invisible and the baseline appeared to make zero hops. It now
  reports `checked` and `skipped` alongside every count.
* **Diagnostics accumulated only on success**, describing the easy points and
  hiding the ones where the seed list was worked hardest.
* **"Too heavily damped" was flattened into "did not converge"**, so
  `flagTooHeavilyDamped` never fired. Fixed by asking the engine's own seed for
  its status.

### Result

Nine sweep lines (three plasmas × `kperp` ∈ {0.005, 0.01, 0.02}), `h = 0.008`,
reference at `h/8`:

| Solver | Branch hops | Cost per line |
|---|---|---|
| `whamp` | **13** | ~0.1 s |
| **`whamp_track`** | **0** | 0.5–1.9 s |
| `whamp_global` | 7 | ~80–95 s |

`whamp_track` reproduces the reference to **10⁻¹²–10⁻¹¹** where the baseline is
off by 10¹–10². On the three well-behaved beam lines the two are indistinguishable.

`whamp_global` (enumerate all roots per point, match by optimal assignment) is
structurally more appealing and worse on both quality measures at ~800× the cost:
in cold multi-species plasmas the susceptibility expansion fails over most of the
complex plane, so the enumeration cannot sample the box reliably. Retained,
marked **not recommended**, so the next person does not spend a week
rediscovering why.

**"Missing points" is not a quality metric** — a coverage count that rewards a
solver for jumping onto a neighbouring branch is measuring the wrong thing.
`flagSuspectBranch` is the honest version of the same information.

---

## 5. Tracking without hand-tuned scales (`whamp_track`)

> **Naming.** `whamp_track_test.m` is modular v1 (calls `+whampkit`);
> `whamp_track.m` is self-contained v1; `whamp_track.m` supersedes both.

**Headline: v2 is the first solver in this tree that reports no branch hops *and*
no unjudgeable points**, at 2.83 s against 3.02 s for self-contained v1 on the
nine-line benchmark and **1.14× v1** on the large 2-D sweep that first exposed a
real cost.

| solver | hops | in range | past end | time |
|---|---|---|---|---|
| `whamp` | 13 | 42 | 19 | 0.75 s |
| `whamp_track_test` (modular v1) | 0 | 56 | 19 | 7.74 s |
| `whamp_track` (self-contained v1) | 0 | 56 | 19 | 3.02 s |
| **`whamp_track`** | **0** | **57** | **0** | **2.83 s** |

### 5.1 Removing `relFloor` — and the trap in doing it

`relFloor = 0.25·|x|` was v1's least principled number: on a flat branch the
reference steps ~0.004 while the scale floored at 0.247, giving a jump tolerance
of 1.48 — about 60× too large. v2 replaces it with the exact implicit step

$$\Delta x = -\frac{(p-p_0)D_P + (z-z_0)D_Z}{D_x}$$

which `dispersion.m` already returns for the group velocity. The only remaining
floor is `sqrt(eps)·|x|`, which is not a knob but the statement that a step below
the root's own round-off cannot be resolved. Median tolerance on the flat branch
falls from 1.48 to 0.047, roughly 30× tighter.

Three findings came out of this, each more useful than the change itself:

* **Sizing a continuity tolerance from a prediction is circular.** Doing the
  obvious thing — sizing the *rejection* threshold from the prediction's own
  recent error — reproduced the baseline's **13 hops exactly**. A tolerance
  derived from the prediction is an allowance for the prediction being wrong; if
  it is wrong because the branch is about to be lost, the tolerance widens by
  exactly enough to forgive losing it. v2 sizes the tolerance from the *step the
  branch actually took*, which is an observation, not a restatement.
* **`dispersion.m` reports `DP = 0` and `DZ = 0` on the axis planes**, for its own
  `1/k` reasons — and every sweep line in the suite starts at `kpar = 0`, where
  the branch genuinely moves by a factor of 2.4. A zero derivative is the
  *absence* of a prediction, not a prediction of no movement. `local_predict`
  therefore returns an `informative` flag and falls back to the grid anchor.
  Suppressing the prediction entirely also failed — it let v2 accept the engine's
  hop, because with no prediction there is nothing for continuity to disagree
  with.
* **The model residual `|root − prediction|` is the wrong yardstick.** On a
  completely clean line it grows from 3·10⁻⁴ to 0.11 purely from local curvature,
  producing 7 false suspect flags. A *step* is the right yardstick, because a hop
  lands a branch-spacing away while curvature stays small compared with one step.
  v2 measures the rate `|dx|/|dk|`.

### 5.2 Certification by the argument principle

Every accepted root is certified by a winding-number count of `D` around a small
circle — an adaptive polygon starting at 8 samples, with **no residual threshold
anywhere**. Since `D` is meromorphic the count is `Z − P`: `+1` is a simple zero
and the root is real; `0` means the reported point is **not a root at all**; `−1`
means it sits on a species gyroharmonic. This works precisely where a residual
test cannot, because `|D|` ranges over ~30 orders of magnitude while the winding
number is scale-free. On the regression line **99 of 99 accepted roots certify as
`+1`**. Certification is not a lower bound — it is a direct statement about each
root, independent of any reference, which removes the caveat from the per-point
claim. Three separate flags keep "wrong branch", "not a root" and "unknown" apart.

### 5.3 The fold, and why v1's extra coverage was not coverage

Walking the branch at 32× resolution showed that on the anisotropic plasma at
`kperp = 0.005` the tracked branch **terminates in a fold at `kpar ≈ 0.0425`**,
where `dx/dkpar → 0` and `x` plateaus at 0.970741. Past it there is no root
connected to this branch.

The step-refined reference is itself a continuation and dies at the same place, so
points beyond it are not "hard points the reference failed on" — they are points
refinement *cannot judge at all*. v1 reports 4 points past the fold, which are
necessarily on another branch, reached because `relFloor` kept the tolerance wide
enough for a lucky scan-offset seed; v1's own continuity check cannot object
because the floor calibrates the check. `pastEnd` is therefore reported as a
first-class metric, and it is the number v2 takes to zero. **Removing `relFloor`
did not lose coverage; it declined a claim that was never earned — and gained one
point of genuine coverage (57 vs 56).**

### 5.4 Two-dimensional sweeps: the item listed as a test found the real bug

v1 keeps its step history in **visit order**. On a genuine 2-D sweep the
previously visited point at the start of a new line sits at the far end of the
previous line, at a completely different wave number — so the first point of
every line after the first is judged against a step between unrelated points.
With grid anchoring alone, `varyKzFirst = 0` gave first-column roots
`0.409 0.878 12.55 26.58 45.82 69.96` with 17 suspect flags, and the two sweep
orderings disagreed by **7.6×10¹**. Two causes:

* **Pooling steps across axes.** The branch moves ~70× faster in `kpar` than in
  `kperp`, so a pooled scale was meaningless on one axis. Fixed by a **rate per
  axis**, with `dk = |p − p₀| + |z − z₀|` — a symmetric form needing no knowledge
  of which axis varies (assuming the inner axis was always `kperp` fixed one
  ordering and broke the other).
* **A windowed median is not a local rate.** One point came out as `0.4887` where
  the true root is `0.9046`. A full root scan at every grid point proved
  `0.48871` **is not a root at all**. The tolerance was too wide because a median
  over a *visit-order window* was drawn from the flat tail of the previous line
  (`dx/dkpar ≈ 0.06`) while the point being entered sat near `kpar = 0`
  (`dx/dkpar ≈ 57`). **A window in visit order is not a neighbourhood in the
  grid.** Fixed by a local rate map `rateMap(ip, iz, axis)` with Manhattan-ring
  lookup.

Result: both orderings agree to **3.19×10⁻¹¹** over 24 common points, 2-D sweeps
match single-line runs **exactly**, and 1-D results are bit-identical to v1.

### 5.5 `whamp_track_plot.m`

Draws a dispersion diagram on which branch confidence is visible: filled dot =
found and certified, open circle = not certified, grey square = certification
unavailable, red ✗ = suspect branch, dashed overlay = the implicit prediction
against the accepted root. Gaps are never bridged (`'LineStyle','none'`). Two
orientation facts cost more time than the plotting: **the tracker returns arrays
as (kpar × kperp)**, so a fixed-`kperp` curve is a *column*, and the helper's
field accessor must apply the same transpose as `f` or the flags silently
describe the wrong points.

### 5.6 Performance: v2 was reported ~10× slower, and why

Three causes. **The first two were real fixes to real pathologies found by
profiling the wrong workload, and neither helped the reported case at all.**

1. **An expanding-ring search that was only linear when the map was full.** The
   rate-map lookup walked ring by ring bounded by the whole grid, so on an axis
   with *no* recorded rates — the leading edge of a sweep — every lookup
   exhausted every ring: `O(N²)` per lookup. On a 201×201 map, 0.926 ms per
   lookup at 0% fill against 0.157 ms for the hybrid replacement (ring out to a
   cap of 8, then one vectorised minimum-distance scan over occupied cells — a
   *complete* search, so the cap is a fast path, not an approximation). Cost
   tracks distance to the nearest hit, which a real sweep does not control.
   Verified with 8000 randomised lookups against the original walk: 0 mismatches.
2. **The certification contour was fixed at 33 samples of order-2 `D`.** Order 1
   suffices (checked bit-identical in `D` and `IERR` at 40 random states; order 0
   is *not* usable — the susceptibility is never accumulated and `D` comes out as
   garbage). An adaptive contour starting at 8 samples and doubling while the
   largest phase step exceeds π/2 resolves 98% of certifications on the first
   pass, cutting mean evaluations from 32.63 to 9.12. The start of 8 was chosen by
   measurement: starts of 4 and 6 occasionally under-resolve against a fixed-32
   reference. 267 cases at the shipped setting: 0 mismatches.
3. **The actual cost was certification itself: nine full evaluations of `D` per
   root.** On the user's case — a 100×100 2-D whistler sweep — v2 was 2.0× v1,
   and with `certify = false` it was *faster than v1*. Certification was 56% of
   runtime. Over 876 real certified roots the contour resolved at 8 samples 100%
   of the time, so there was nothing left in the adaptive loop: the cost was
   irreducibly nine scalar `D` evaluations, each rebuilding the state and running
   the eight-pole sum from scratch.

   The fix evaluates the whole contour in **one vectorised pass**, since a contour
   is eight or sixteen frequencies at a single shared `(p, z)`. Two things were
   needed to make it pay: freezing converged Taylor samples **by masking rather
   than indexing** (the logical-indexed version was *slower than the scalar loop
   it replaced*, 0.97×; writing the freeze as a zeroed term with a substituted
   harmless denominator is bit-identical and cost four elementwise products,
   taking it to 2.52×), and **stacking all eight poles into one call** (pole-major
   in one 72-element row; grouping cannot change the result because each sample's
   series advances independently). Result on the user's case, interleaved:

   | solver | best |
   |---|---|
   | `whamp` | 18.03 s |
   | `whamp_track` (v1) | 21.42 s |
   | v2, `certify = false` | 21.97 s |
   | **v2, `certify = true`** | **24.49 s** |

   v2 went from 2.0× v1 to **1.14× v1**, and certification from 56% of runtime to
   ~10%, **bit-identically** over 12 564 compared values of `D` and `IERR`.

An earlier non-interleaved run appeared to show 16 contour samples being *cheaper*
than 8, which would have implied the adaptive loop cost more than it saved. It was
machine drift charged to whichever configuration ran first; interleaved over four
repetitions the ordering is monotone. **The methodological lesson: a benchmark on
a small converged example cannot see a shape-dependent bug — *and* fixing the
shape-dependent bugs still does not license a claim about total cost. The workload
that prompted the complaint has to be the one that is measured.**

### 5.7 A pre-existing non-termination bug, deliberately not fixed

`whampkit.axis_spec` accepts a degenerate axis `[v 0 v]` (start = end, step 0) and
`fill_axis` correctly returns a single-point axis. But `advance_sweep` walks the
outer axis by repeated addition of the step, so with a zero step the position
never changes and the termination test stays true forever. `whamp.m` loops forever
with no message; the three track solvers carry a `guard > 200000` backstop and
raise `Runaway` after tens of minutes. Writing a single `kperp` as
`[0.005 0 0.005]` instead of `0.005` is enough to trigger it. The fix is one line
in shared code and is **deliberately not made**, because every claim in this
report rests on `whamp.m` being byte-identical; it should be a separate change
with the engine's own tests.

---

## 6. C back-port: tracker and `mexwrap` MEX gateway

`src/mexwhamp.F` had never been ported. It is now reimplemented in C as
`mexwrap.c`, exposing exactly the struct interface the MATLAB code uses:

```matlab
Output = mexwrap(PlasmaModel, InputParameters)          % optionally (…, Options)
```

`PlasmaModel` needs `B` and `Species`; `InputParameters` needs `fstart`, `kperp`,
`kpar`, `useLog`, with `varyKzFirst` and `maxIterations` optional. The output
struct carries `kperp`, `kpar`, the 16 result arrays (`f`, `Ex`…`u`), and the
flags `flagSolutionFound`, `flagTooHeavilyDamped`, `flagNoConvergence`,
`flagSuspectBranch`, `certWinding`, `flagNotCertified`, `flagCertUnknown`, plus
`predicted` and `diagnostics`. `matlab/mexwrap.m` is a shim that locates the
built MEX.

The tracker itself is `whamp_track.c` / `whamp_track.h` with shared declarations
in `comtrack.h`; callers configure the existing `comin` globals as for
`whamp_engine()`, optionally call `track_options_default()`, then call
`whamp_track()`. It writes the existing `comoutput` arrays. The original engine
entry points are unchanged — callers opt in.

**Ported to all four C trees** — `srcc64`, `srcc64`, `srccopt`, `srcc` — for
parity. Two subtleties had to be handled:

* **`PI` collision.** `srcc` and `srccopt` declare file-scope `static const double
  PI` in five files, so a shared `PI` macro would make those declarations
  ill-formed. Their `constants.h` omits `PI` and the tracker uses a file-local
  `TR_PI`.
* **The tracker must match *its own* engine.** The bit-for-bit trees emulate
  gfortran with `F32()`, so their engine uses single-rounded `1836.1`, `1e-6` and
  `1e-5`. With plain-double constants the tracker and engine computed *different
  dispersion functions* and disagreed about the root. In `srcc/constants.h` and
  `srccopt/constants.h` these are therefore defined via `((double)(float)...)`;
  the all-double trees keep them unrounded, matching `whamp_track.m`.
  `DEBYE_KM` (12405) and `ELECTRON_REST_ENERGY_KEV` (511) are exactly
  representable in single and are written plain in every tree.

**Verification.** `whamp_track_test` checks both 2-D sweep orders, winding
certification, and exact agreement with the legacy engine's first root; it passes
in all four trees. `whamp_track_dump` prints each case's arrays at `%.17g` for
comparison against MATLAB by `srcm/whamp_track_vs_c.m`, whose **14 cases all
MATCH (worst relative deviation 2.834×10⁻¹²) with every flag matching exactly**.
Dump output is byte-identical between `srcc64` and `srcc64`. `srcc` and
`srccopt` differ on 9 of 14 cases — a pre-existing difference between the two
bit-for-bit trees, not a tracker issue: `srccopt` applies `F32()` at two sites
`srcc` does not, one of them a routing threshold in `xsi.c` that selects the
integration path (and they use different input-editing helpers, `ctof.c` versus
`fedit.c`). Each tree's tracker agrees exactly with its own engine, which is what
`whamp_track_test` asserts. All four MEX builds return identical results on a
common problem.

Two divergences this comparison found were fixed in `whamp_track.c`:
`tr_predecessor` now walks back at most 64 cells along the outer axis as
`local_predecessor` does, and `tr_winding` returns the half-integer-resolved
winding on an early resolved contour instead of rounding it. Neither is reachable
on the 14 cases, which is why both were invisible to the flag comparison; they
matter on grids wider than 65 cells and for contours whose phase lands on a
half-integer.

### 6.1 Which engine a result came from

Four trees build the same gateway, and they are *not* numerically
interchangeable, so "which one ran" has to be answerable after the fact rather
than guessed from a path. Three pieces:

* **Stamp.** Each Makefile defines `VARIANT` and passes
  `-DWHAMP_VARIANT_NAME=$(VARIANT)` when building the MEX. The name is
  stringified inside `mexwrap.c` rather than quoted on the command line,
  because `-D` quoting differs between shells. A build without the flag still
  compiles and reports `unknown_srcm` — never a wrong answer.
* **Report.** Every result carries `Output.engine`, and
  `mexwrap('engine')` returns the name without running anything.
* **Refuse.** `matlab/mexwrap.m` reads `WHAMP_ENGINE` (default `srcc64`),
  loads *only* that directory's MEX, and asks the loaded binary to confirm.
  A mismatch, an unstamped binary, or a missing build is an error naming what
  was asked for, what was found, and where. The earlier shim walked a
  preference list `{srcc64, srccopt, srcc64, srcc}` and took the first hit,
  which silently substituted a different engine whenever one build went
  missing — that is what this replaces.

Self-containment was checked rather than assumed: `nm -u` on all four MEX
files shows zero undefined engine symbols, so each links only the
`libwhamp.a` built beside it. The risk was never cross-tree linking; it was
selection.

Switching variants mid-session needs `clear mexwrap`, since MATLAB keeps a
MEX loaded until it is cleared; the shim reports that as an error rather than
pretending to switch.

The Poynting flux and `SGP`/`SGZ` need a yardstick other than the array's own
magnitude — they are differences of nearly cancelling products, where a
per-cell ratio reads a difference of 1e-17 as a relative error of one. They are
instead judged by whether the gap could have come from the inputs they are
computed from (`|dS| ≤ |dE||B| + |E||dB| + 2ε|E||B|`); the bound holds
everywhere, so the difference is inherited from the fields, not a logic divergence.

---

## 7. Findings that generalise

1. **Rounding behaviour *is* behaviour.** Replacing `float` with `double` is not
   a cleanup — it is a change of program, justified by evidence like any other.
2. **`-Wall -Wextra` was clean before and after the only real bug found in the
   whole exercise.** `-Wconversion` caught it. Treat warnings as part of
   verification, not housekeeping — especially on a toolchain that does not
   implement C99 type-generic `math.h`.
3. **Differential testing finds what inspection misses** — three control-flow and
   formatting defects in phase 1, and the complex-trig defect only because no
   test in the battery reached the affected code path.
4. **Any trust region must be able to report that it bit.** A cap that silently
   shrinks a step to zero turns a step-size convergence test into a test of the
   cap, and reports the wall of its own disc as a root.
5. **Report the denominator alongside the count.** A harness that quietly skips
   unjudgeable points makes the baseline look perfect on the hardest line.
6. **Sizing a continuity tolerance from a prediction is circular** — the tolerance
   forgives exactly the error it is meant to catch. Measure against the branch.
7. **A derivative reported as zero for the convenience of the code computing it
   is not a prediction of zero.**
8. **A window in visit order is not a neighbourhood in the grid.**
9. **Coverage counts are not quality metrics** when a solver can earn them by
   jumping branches.
10. **Benchmark the workload that prompted the complaint.** Two genuine
    pathologies were found by profiling the wrong one, and fixing both changed
    nothing about the reported cost.

---

## 8. Open items

1. Predict the branch rate per grid cell in a cheap preliminary pass rather than
   recording it from the step that already succeeded — should remove the four
   remaining suspect points under `varyKzFirst = 1`.
2. Let `+whampkit/dispersion.m` report `DP`/`DZ` on the axis planes by a one-sided
   stencil, or return an explicit "not available" flag.
3. Make the reference fold-aware (continue in `(x, k)` through a turning point) so
   the region past `kpar ≈ 0.0425` becomes judgeable.
4. Fix the degenerate `[v 0 v]` axis non-termination — as a separate change with
   the engine's own tests, since it alters `whamp.m`.
5. Close the `av.c` energy-output coverage gap (`enden=`, `enfl_p=`, `enfl_z=`)
   before touching that file again.
6. Decide the sweep-endpoint question of §2.1 explicitly rather than by rounding
   accident.
7. Port the tracker to the Fortran engine, which still carries the original
   single-seed defect.
8. Wire `whamp_verify` and the golden suites into CI so a MATLAB edit is gated the
   same way a C edit is.
9. Consider whether four near-duplicate C trees (`srcc`, `srcc64`, `srccopt`,
   `srcc64`) are worth their long-term maintenance cost now that the
   diagnostic value of the experiment has been extracted.

---

## 9. Layout and reproduction

```
whamp/
├── src/            ORIGINAL Fortran 90 (read-only reference, 19 files)
├── srcc/           C translation, bit-for-bit (F32 emulation, ctof.c)
├── srcc64/         all-double + tidied; the validated tracker/MEX reference
├── srcm/           pure MATLAB: whamp.m, whamp_track{,_test,_v2}.m, +whampkit/
├── matlab/         MATLAB front end, mexwrap.m shim
├── fortran_ref/    gfortran oracle + differential harness
├── tools_golden/   30-session golden suites (golden/, golden64/)
├── Models/         plasma model files (Ex1, Ex2, …)
└── doc/            documentation, incl. test.txt (acceptance spec)
```

Each C tree builds with `mingw32-make` (`make` is not on PATH) and exposes
`all`, `whamp_track_test`, `whamp_track_dump` and `mexwrap`. `make clean` uses
`rm` and does not work on Windows; use `git clean -fdX`.

```powershell
# C: build a tree and run its tracker self-test
cd srcc64 ; C:\MinGW\bin\mingw32-make.exe all whamp_track_test whamp_track_dump
.\whamp_track_test.exe            # exit 0 on success

# C MEX (must be run from inside the tree; make stamps the variant name)
C:\MinGW\bin\mingw32-make.exe mexwrap

# Tracker vs the MATLAB reference (14 cases; driver path is srcc64)
& matlab.exe -batch "cd srcm; addpath(pwd); whamp_track_vs_c"

# MATLAB regression and conformance
& matlab.exe -batch "addpath('srcm'); whamp_verify('baseline.txt')"
& matlab.exe -batch "addpath('srcm'); whamp_vs_repl_golden; whamp_vs_srcc64"

# Golden suites (exit 1 on any difference)
powershell -ExecutionPolicy Bypass -File tools_golden\run_cases.ps1 `
    -Exe srccopt\whamp.exe -Cases tools_golden\cases.json `
    -OutDir tools_golden\actual -CompareDir tools_golden\golden

# Branch-tracking benchmark (nine sweep lines)
& matlab.exe -batch "addpath('srcm'); compareExample"
```

The oracle must be built with `-O2`; `src/Makefile` omits it and yields a
different binary. The tracking harnesses (`verify_*`, `bench_*`, `eval_*`,
`diag_*`) live in the session workspace alongside the v2 head file and the
assembler, and are not part of the tracked tree.

---

## 10. What is claimed, and what is not

**Claimed.** `srccopt/` reproduces `srcc/` byte-for-byte on 30 sessions and the
documented acceptance test exactly. `srcm/whamp.m` produces byte-identical output
to its pre-refactor state across 30 cases and 21 fields at full precision, and
that state matches `srcc64` on all 9 comparable golden sessions. Over nine
sweep lines spanning three multi-species plasmas and three perpendicular
wavenumbers, `whamp_track` records **0 branch hops and 0 points beyond the
reference's reach** where `whamp.m` records 13 and 19; every accepted root it
reports is independently certified as a simple root of `D`; it agrees with itself
across sweep orderings to 3·10⁻¹¹; and it contains no hand-tuned scale. The
tracker and `mexwrap` build and self-test in all four C trees, and the tracker
matches the MATLAB reference on 14 cases with every flag identical. `whamp.m` is
byte-identical to its pre-work state and both v1 builds are unchanged in
behaviour.

**Not claimed.** That zero hops holds generally — three plasmas and three `kperp`
values are not a proof, and the *comparison* between solvers still rests on a
reference that stops at the fold. That `srcc64` matches the original (it does
not, in two specific documented ways). That the golden suites are complete — §2.3
and §2.4 are direct evidence that they are not. That the residual oracle
differences are understood beyond being ULP-level and pre-existing. That the
MATLAB port is a faithful reference for the Fortran original — it is a faithful
reference for `srcc64`. And that anything is fixed for users of the Fortran
engine, which still carries the original single-seed defect.
