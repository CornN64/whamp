# WHAMP in Pure MATLAB (`srcm/whamp.m`): Refactor and Conformance Report

AI-assisted work by Walter Puccio 2026, extending
[`REPORT64.md`](REPORT64.md) (which documents the `src/` → `srcc/` → `srccopt64/`
C translations and their byte-for-byte golden verification).

## Abstract

This report documents the refactoring of [`srcm/whamp.m`](srcm/whamp.m), a pure
MATLAB port of the WHAMP plasma-wave dispersion solver, and the evidence that
the refactor changed no behaviour.

The method follows the C work: build a bit-exact regression harness first,
confirm the harness is anchored to the real reference program rather than to
the MATLAB code's own output, then refactor in small batches and re-diff after
every one.

* **30 cases** covering every code path dump **21 output fields at 20
  significant digits** — 2616 lines of golden text.
* The pre-refactor code was checked against the actual `srccopt64` transcripts
  in [`tools_golden/golden64/`](tools_golden): **9 of 9 comparable sessions
  match**, the closest at **3.2 × 10⁻¹⁴ relative error**.
* After 9 refactor batches the dump is **byte-identical** to the baseline.
* Static analysis (`checkcode`) improves from **7 warnings to 0**.

Three divergences from the C reference remain. All are deliberate, documented at
the point of divergence in the source, and unchanged by this work (§4).

---

## 1. What was refactored

[`srcm/whamp.m`](srcm/whamp.m) grew from 938 to 1129 lines (29.8 → 40.4 KB).
The increase is documentation and extracted helpers, not new behaviour: 195
comment lines, and 19 functions became 34.

| Area | Change |
|---|---|
| Input validation | Extracted `normalise_species`, `is_finite_real_scalar`; removed dead branches |
| Axis handling | Documented `axis_spec` / `fill_axis` / `sweep_origin`; dropped an unreachable `count < 1` test |
| `prepare_plasma` | Documented REN/T/ST normalisation, the electron unit-mass convention, loss-cone collapse, relativistic CV |
| `solve_point` | Split into `solve_point` / `dispersion_state` / `newton_root` / `evaluate_at`; magic numbers named (`freqRelTol`, `freqAbsTol`, `minimaTol`) |
| `dispersion` | Triplicated `DC` expression collapsed into `dcTerm(E,col)`; locals `av/bv/cv` → `Av/Bv/Cv` to match `difu.c` and avoid shadowing `plasma.cv` |
| `susceptibility` | Eight-pole residue tables hoisted into `pole_constants`; unused `ib`/`b`/`delta` arguments dropped; the four derivative stages labelled |
| `wave_energy` | Documented the cofactor construction (Melrose, *Plasma Astrophysics* Vol I p.43) and the 1994-04-13 L-mode normalisation |
| R-function | Documented `r_function` / `r_taylor` / `r_asymptotic` / `r_integral`; Gauss–Legendre abscissae hoisted into `gauss_legendre_16` |
| Main function | Split into `species_arrays`, `sweep_state`, `remember_point`, `advance_sweep`, `next_frequency`, `record_point`. The 14 file-scope predictor variables of the C engine are now named fields of a `sweep` struct |

## 2. Test method

Two independent layers, both run headless under MATLAB R2020b.

### 2.1 Bit-exact regression harness

| File | Role |
|---|---|
| [`srcm/whamp_dump.m`](srcm/whamp_dump.m) | Runs 30 cases, writes 21 fields per case at `%.20g` |
| [`srcm/whamp_verify.m`](srcm/whamp_verify.m) | Re-dumps and byte-diffs against `baseline.txt`, printing the first 25 differing lines |
| [`srcm/baseline.txt`](srcm/baseline.txt) | 2616-line golden dump taken before any edit |

The case set exercises: single points; both sweep orders (`varyKzFirst` 0 and
1); negative steps; linear and logarithmic axes; three-species and loss-cone
plasmas; ion drift; heavy damping; non-convergence; `maxIterations` of 1 and 3;
`kperp = 0` and `kpar = 0` (the guarded `DP`/`DZ` branches); magnetic fields of
10, 100 and 300 nT; partial species defaults; and 7 error paths. Every one of
the 7 expected error identifiers (`InvalidMaxIterations`,
`InvalidStartFrequency`, `InvalidLogOption`, `InvalidWaveNumberStep`,
`InvalidSpecies`, `InvalidMagneticField`, `TooManySpecies`) is reproduced.

Complex values are dumped as real and imaginary parts at full precision, so a
single flipped low-order bit in any output is a diff.

### 2.2 Reference conformance

The harness alone would only prove the code still equals *itself*. Two further
scripts compare against the real `srccopt64` output:

* [`srcm/whamp_vs_srccopt64.m`](srcm/whamp_vs_srccopt64.m) — parses the `fOUT=`
  complex vectors from the engine sessions in
  [`tools_golden/golden64/`](tools_golden).
* [`srcm/whamp_vs_repl_golden.m`](srcm/whamp_vs_repl_golden.m) — parses the
  `kp kz fre fim` result rows from the REPL transcripts (7 significant digits
  in the golden text, so the tolerance there is looser by construction).

Model parameters are translated into `whamp.m` units with
`B = gyro / 0.0279928`. `Models/Ex1` uses a **1 eV** electron temperature and
`Models/Ex2` an electron-first species order; neither matches `whamp.m`'s
built-in defaults, and getting this mapping wrong produces a ~3 × 10⁻⁴ mismatch
that looks like a code bug but is not.

## 3. Results

### 3.1 Conformance to `srccopt64`

| Case | Points | Max relative | Max absolute | Verdict |
|---|---:|---:|---:|---|
| `engine_doc` | 1 | 3.228e-14 | 3.435e-15 | **MATCH** |
| `engine_sweep` | 5 | 3.228e-14 | 3.435e-15 | **MATCH** |
| `engine_log` | 15 | 7.833e-15 | 3.970e-15 | **MATCH** |
| `engine_change` | 1 | 3.228e-14 | 3.435e-15 | **MATCH** |
| `doc_pzf` | 1 | 6.563e-10 | — | **MATCH** |
| `sweep_p` | 6 | 5.821e-08 | — | **MATCH** |
| `sweep_neg_step` | 5 | 5.817e-08 | — | **MATCH** |
| `log_scale` | 15 | 1.474e-08 | — | **MATCH** |
| `no_convergence` | 1 | 1.028e-09 | — | **MATCH** |

The engine rows compare against full-precision output and agree to roughly
double-precision round-off. The REPL rows are limited by the 7-significant-digit
golden text.

### 3.2 Regression after refactoring

```
>> whamp_verify('baseline.txt')
IDENTICAL to baseline.txt
```

All 2616 lines match byte-for-byte after the final batch. The same result held
after each of the 9 intermediate batches, which is what made the incremental
approach safe.

### 3.3 Static analysis

| | Functions | `checkcode` warnings |
|---|---:|---:|
| Before | 19 | 7 |
| After | 34 | **0** |

All seven warnings were removed, and the regression stayed byte-identical while
they were.

**Four dead increments** in the sweep loop:

```
186: The value assigned to variable 'ip' might be unused.
193: The value assigned to variable 'iz' might be unused.
202: The value assigned to variable 'iz' might be unused.
209: The value assigned to variable 'ip' might be unused.
```

These are `ip = ip + 1` / `iz = iz + 1` in the four branches that advance the
inner and outer axes. They are dead because the top of the loop recomputes both
from the axis values:

```matlab
while true
    ip = grid_index(sweep.PLG, kpSpec, numel(kp));   % overwrites
    iz = grid_index(sweep.ZLG, kzSpec, numel(kz));
    ...
    ip = ip + 1;      % never read again before the next iteration
end
```

The increments are leftovers from an earlier sequential-counter formulation that
`grid_index` superseded. `grid_index` derives the row from the walked axis value
instead, and additionally clamps to `[1, count]` so a rounding slip cannot index
out of range. A probe replaying the walk over linear, negative-step and
logarithmic axes found `grid_index(next) == grid_index(current) + 1` at every
step, so the increment was **redundant rather than wrong** — `grid_index` is
strictly safer, but the two agree everywhere the sweep actually goes. **Fix:
delete them**, which the extraction of `record_point` and `advance_sweep` did
naturally.

**Three unused arguments** to `susceptibility`, all vestiges of the C signature
`CHI(XSI, J, IB, KOL, IERR)`:

```
645: Input argument 'ib' might be unused, although a later one is used.
645: Input argument 'b' might be unused.
645: Input argument 'delta' might be unused.
```

* `ib` — in C, `CHI` reads `A = AA[J][IB]` and so genuinely needs the column
  index. In the port the caller (`dispersion`) already loops over `ib` and passes
  the selected value in as `a`, making the parameter redundant.
* `b` and `delta` — `CHI` has no such arguments. They were added to the port
  speculatively and never referenced; their effect lives entirely in the caller's
  `factor * q` weighting.

**Fix: drop them**, giving `susceptibility(state, j, order, a)`, with a comment
recording where the C `IB` went.

Both classes look more dangerous than they are, which is the argument for having
the harness before touching them: deleting an increment inside a sweep loop is
exactly the kind of edit a tolerance-based test would wave through, and the only
convincing evidence that it is inert is that all 2616 dump lines still match.

## 4. Deliberate divergences from the C reference

These are **not** introduced by this refactor; they predate it, are commented at
the point of divergence in [`srcm/whamp.m`](srcm/whamp.m), and were re-confirmed
rather than "fixed".

1. **No cross-point state leakage.** The C engine keeps `D`, `DX` and `SG` in
   file scope, so a failed point continues with values left by an unrelated
   earlier point, and NaN frequencies leak between sweep points. `whamp.m`
   resets per point. This is why `engine_ex2` is excluded from the conformance
   table — its golden output depends on the leakage.
2. **Heavily-damped points return immediately.** The C sets `IERR` and carries
   on with a Newton step derived from stale globals; `whamp.m` returns status 2.
3. **Grid indices are clamped.** `grid_index` bounds the index to `[1, count]`;
   the C writes out of bounds for negative steps.

Each is a case where the C behaviour is a memory-lifetime artifact rather than
intended physics. Matching it would require reproducing undefined behaviour, and
the first point of every sweep agrees either way.

## 5. Reproducing

```powershell
cd srcm
& "C:\Program Files\Polyspace\R2020b\bin\matlab.exe" -batch "addpath(pwd); whamp_verify('baseline.txt')"
& "C:\Program Files\Polyspace\R2020b\bin\matlab.exe" -batch "addpath(pwd); whamp_vs_repl_golden; whamp_vs_srccopt64"
```

## Conclusions

**The refactor is behaviour-preserving, and this is verified rather than
asserted.** The strongest evidence is not that the output is unchanged — that
only follows from the harness being correct — but that the harness was anchored
to the external `srccopt64` reference *before* any edit was made, at
3 × 10⁻¹⁴ relative agreement. The refactor therefore preserves agreement with
the reference program, not merely self-agreement.

**Three things made the incremental approach worth the overhead.** Dumping at 20
significant digits rather than comparing with a tolerance turns verification into
a byte diff that cannot hide a regression behind a generous epsilon. Batching one
function group per run keeps every failure attributable to a known edit — it
caught a genuine semantic error during the `solve_point` split, where a "tidied"
re-evaluation was substituted for the C engine's deliberate step-back, which
restores `X` but intentionally keeps the stale `D` for the local-minima test. And
running the harness after each batch, rather than at the end, meant that bug was
found while its cause was still the only thing that had changed.

**The conformance scripts caught a mistake in the test, not in the code.** The
first comparison attempt showed a ~3 × 10⁻⁴ disagreement that looked like a
numerical defect. It was a units error: `Models/Ex1` specifies a 1 eV electron
temperature, not `whamp.m`'s built-in 100 eV. A refactor verified only against
itself would have shipped that blind spot.

**What is claimed:** `srcm/whamp.m` produces byte-identical output to its
pre-refactor state across 30 cases and 21 fields at full precision; that
pre-refactor state matches `srccopt64` on all 9 comparable golden sessions; and
`checkcode` is clean. **What is not claimed:** that the golden set covers all of
WHAMP (it does not exercise relativistic beams, all ten species slots, or the
`enden=` field-energy output flagged in §3.3 of [`REPORT64.md`](REPORT64.md)),
that the three documented divergences in §4 are the right engineering call for
every consumer, or that the MATLAB port is a faithful reference for the Fortran
original — it is a faithful reference for `srccopt64`.

**Suggested next steps**, in rough order of value:

1. Add golden cases for the uncovered paths above, particularly multi-species
   and relativistic configurations, before further changes to
   `susceptibility`.
2. Decide explicitly whether the §4 divergences should be configurable — a
   `strictReference` option would let `whamp.m` reproduce `engine_ex2` for
   differential testing without giving up deterministic behaviour by default.
3. Wire `whamp_verify` into CI alongside the existing
   [`tools_golden/run_cases.ps1`](tools_golden/run_cases.ps1), so a MATLAB edit
   is gated the same way a C edit already is.
