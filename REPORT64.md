# WHAMP C Translation: Optimised Variants (`srccopt`, `srccopt64`)

AI-assisted work by Walter Puccio 2026, extending
[`REPORT.md`](REPORT.md) (which documents the original `src/` → `srcc/`
translation and its byte-for-byte verification).

## Abstract

This report documents two derived C trees built from the byte-for-byte verified
translation of the WHAMP plasma-wave code in [`srcc/`](srcc), and the
differential-testing evidence for each. The method is unchanged from
[`REPORT.md`](REPORT.md): build with a fixed toolchain, replay 30 scripted
sessions through the REPL and the callable engine, and compare captured stdout
byte-for-byte against a golden baseline.

Two trees were produced.

* [`srccopt/`](srccopt) — a tidied and optimised C rewrite that reproduces
  `srcc/` **byte-for-byte** on all 30 golden sessions.
* [`srccopt64/`](srccopt64) — an experimental **all-`double`** variant built on
  top of `srccopt/`, in which every Fortran single-precision entity is widened
  to C `double` and all constants are collected in
  [`srccopt64/constants.h`](srccopt64/constants.h). It is verified against its
  own baseline, `tools_golden/golden64/`, on the same 30 sessions.

Both trees build warning-free with `gcc -std=c11 -Wall -Wextra -O2`, and both
pass their respective golden suites (25 REPL sessions + 5 engine sessions
each). The `doc/test.txt` acceptance session is reproduced exactly by
`srccopt/`; `srccopt64/` reproduces it to one digit in the eighth decimal
place, which is the intended and documented consequence of removing single
precision.

The all-`double` experiment produced four findings that are more valuable than
the refactor itself:

1. **Single precision is load-bearing in the original.** It is not incidental
   rounding; it changes the number of grid points a wavenumber sweep emits
   (§3.1).
2. **The Fortran source disagrees with itself about π and about the speed of
   light**, and the disagreement is observable (§3.2, §3.3).
3. **The `enden=` field-energy output is numerically meaningless** in the
   original code, because of catastrophic cancellation (§3.3).
4. **The translation silently discarded the imaginary part of three complex
   trigonometric expressions** — a genuine mistranslation inherited from
   `srcc/`, found only by compiling with `-Wconversion`, and now corrected in
   all three trees (§3.5).

The conclusion is therefore narrow but firm: **`srccopt/` is the tree to use**,
because matching the reference program is what makes the translation verifiable,
and `srccopt64/` — though cleaner and arguably closer to mathematical intent —
is not closer to reference behaviour. The constants extraction is the part worth
keeping in either tree.

---

## 1. What was done

### 1.1 `srccopt/` — tidy-up and optimisation, byte-exact

Derived from `srcc/` with these changes:

| Change | Where | Effect |
|---|---|---|
| Deleted dead code | `ctof.c` | Never called from any translation unit. |
| Extracted Fortran output editing | new `fedit.h` / `fedit.c` | The gfortran edit-descriptor emulation is now one reusable module instead of being inlined at each call site. |
| Slimmed the central header | `whamp.h` | Was a catch-all; now declares only what is shared. |
| Rewrote argument parsing | `whamp.c` | Helpers `get_padded_arg`, `write_padded_arg`, `adjustl_trim` and a named `FORTRAN_ARGLEN` replace open-coded padding arithmetic. |
| Stride-aware record reader | `input.c` | `read_reals()` replaces a `goto`-driven reader. |
| Allocation macros + axis helpers | `whamp_engine.c` | `ALLOC1`/`ALLOC2`/`axis_size`/`fill_axis` collapse 21 hand-written `free()` calls on the error paths into one unwinding sequence. |

`srcc/` is left untouched so that it remains the reference for differential
comparison.

### 1.2 `srccopt64/` — the all-`double` experiment

`srcc/` (and therefore `srccopt/`) reproduces gfortran exactly, which requires
emulating Fortran's default-`real` entities. That emulation is the `F32()`
macro (`#define F32(x) ((double)(float)(x))`) plus a handful of genuine C
`float` variables in `typin.c`.

`srccopt64/` removes the emulation:

* every `F32(literal)` became the plain literal;
* `typin.c`'s `float DEC`, `DEK`, `TV[3]` and `inputNumber` became `double`
  (these are the only single-precision variables in the whole Fortran source —
  `src/typin.f90` lines 12 and 24);
* the `(float)` casts and the now-redundant `(double)` casts were dropped;
* the `F32` macro was deleted from `whamp.h`.

### 1.3 Constants extraction (`srccopt64/constants.h`)

All constants and shared macros now live in
[`srccopt64/constants.h`](srccopt64/constants.h) (167 lines, 20 public
defines), included through `whamp.h` so no `.c` file needed a new include.

| Group | Defines |
|---|---|
| Mathematical | `PI` |
| Physical | `C_LIGHT`, `C_MS`, `EPSILON_0`, `PROTON_MASS_RATIO`, `DEBYE_KM`, `ELECTRON_REST_ENERGY_KEV` |
| Algorithmic | `FREQ_REL_TOL`, `FREQ_ABS_TOL`, `MINIMA_TOL`, `RTAY_SCALE`, `ASYMP_TOL`, `OVERFLOW_CEILING`, `SWEEP_END`, `POYNTING_COEF` |
| Shared macros | `XS(I,K)` (was duplicated in `av.c`, `difu.c`, `xsi.c`), `RC11`/`RC21`/`RC12`/`RC22` (was duplicated in `rint.c`, `rtay.c`, `rasy.c`) |

Five duplicate `static const double PI` definitions were removed. `1022.` is
now written `2. * ELECTRON_REST_ENERGY_KEV`, and `C_MS` is derived as
`C_LIGHT * 1.0e6`; both rewrites were checked to be bit-identical to the
literals they replace.

Truly file-local defines were deliberately left in place: `FORTRAN_ARGLEN`
(`whamp.c`), `NFIELD` (`input.c`), `FCOMPAT_BUF` (`fedit.h`), and the
`F2D`/`ALLOC1`/`ALLOC2` macros already in `comoutput.h`.

The `Makefile` header dependency was widened from `whamp.h` alone to
`whamp.h constants.h fedit.h`; a touch of `constants.h` was confirmed to
recompile all 18 objects and relink both executables.

---

## 2. Validation

### 2.1 Harness

[`tools_golden/`](tools_golden) drives a build with the same arguments, working
directory and piped stdin, captures raw stdout, and compares byte-for-byte
against a golden directory. `cases.json` holds 30 sessions: 25 for the REPL
(`whamp.exe`) and 5 for the callable engine (`whamp_engine_test.exe`). The
sessions cover the `doc/test.txt` acceptance run, all 19 output-format letters,
bad-input recovery, linear and log wavenumber sweeps, negative steps,
non-convergence, heavy damping, three-species plasmas, mid-session plasma
changes, `-debug`, `-maxiterations`, `-help`, and the no-model path.

### 2.2 Results

| Build | Suite | Result |
|---|---|---|
| `srccopt/whamp.exe` | `tools_golden/golden/` | **25 / 25 byte-identical** |
| `srccopt/whamp_engine_test.exe` | `tools_golden/golden/` | **5 / 5 byte-identical** |
| `srccopt64/whamp.exe` | `tools_golden/golden64/` | **25 / 25 byte-identical** |
| `srccopt64/whamp_engine_test.exe` | `tools_golden/golden64/` | **5 / 5 byte-identical** |

All three trees compile with **zero warnings** under
`gcc -std=c11 -Wall -Wextra -O2` (gcc 6.3.0, `C:\MinGW\bin`), and also under
`-Wall -Wextra -Wpedantic`; see §2.6 for the strictest audit.

Two caveats on the word *byte-identical*. First, these figures are against the
golden directories, which are regression baselines rather than the reference
program — agreement with the actual Fortran oracle is measured separately in
§2.5. Second, both golden directories were regenerated for one case
(`heavy_damping`) after the mistranslation described in §3.5 was fixed; the
pre-fix baseline is reproduced by reverting that one commit.

### 2.3 `doc/test.txt`

The acceptance session in [`doc/test.txt`](doc/test.txt) is
`whamp -file Models/Ex1` followed by `p0z.0022f.1` / `pzf`, whose documented
result is

```
0.0000000 0.0022000 1.0640799E-01 3.34E-11
```

| Build | Output | Match |
|---|---|---|
| `srcc/` (gfortran reference) | `1.0640799E-01  3.34E-11` | yes |
| `srccopt/` | `1.0640799E-01  3.34E-11` | **yes, byte-exact** |
| `srccopt64/` | `1.0640800E-01  3.34E-11` | one digit in the 8th decimal |

That single digit is the deliberate, expected consequence of removing single
precision, and it is the only deviation from the documented behaviour.

### 2.4 `srccopt64` versus the byte-exact golden

`golden64/` is a fresh baseline, so it proves `srccopt64` is *self*-consistent.
To measure what the experiment actually cost, `golden64/` was compared against
`golden/`:

| Category | Count | Meaning |
|---|---|---|
| Identical | 9 | No printed value changed. |
| Digit-level | 17 | Differ only in the last printed digit; max relative difference **1.8 × 10⁻⁷**. |
| Debug-trace only | 3 | `debug`, `plus_minus_debug`, `unknown_opt`. Intermediate Newton iterates differ by up to 1.4 × 10⁻³; the converged root differs only at the 1.8 × 10⁻⁷ level. |
| Structural | 1 | `sweep_p` — see finding 1 below. |

The 1.4 × 10⁻³ figure deserves comment because it looks alarming and is not.
It appears in `D=` values of order 10¹³ on Newton iteration lines. Both variants
stop iterating at the same tolerance while `|D|` is still ≈ 7.7 × 10¹³, so the
intermediate traces legitimately differ; tightening the tolerance to 10⁻¹⁶ in
**both** trees makes them converge to the same root
(`1.0640799E-01  3.70E-15`). The difference is an early-stopping artefact, not a
precision advantage in either direction.

### 2.5 Comparison against the Fortran oracle

The tables in §2.2 prove that each tree matches **its own golden directory**,
which is a regression baseline, not the reference program. To measure agreement
with the original, all three C trees were replayed against
`fortran_ref/whamp_ref.exe` — the authoritative oracle, built with
`gfortran -O2` following the recipe in [`REPORT.md`](REPORT.md) §8. (A scratch
build using `src/Makefile`, which omits `-O2`, is a *different binary*; both
disagree with the golden on the same cases, so the conclusion is unaffected, but
future work should use `fortran_ref/whamp_ref.exe`.)

| Tree | Sessions differing from the oracle | Nature |
|---|---|---|
| `srcc/` | 3 of 25 — `all_outputs`, `log_scale`, `no_convergence` | pre-existing, accepted |
| `srccopt/` | 3 of 25 — same three | pre-existing, accepted |
| `srccopt64/` | 18 of 25 | the float→double experiment (§2.4) |

The three shared differences are the single-ULP perturbations in
non-convergent or ill-conditioned iterations that [`REPORT.md`](REPORT.md) §9
already records as accepted: `log_scale` differs in the last printed digit
(`-4.21E-12` versus `-4.20E-12`), `no_convergence` shows `5.13E-16` versus
`5.59E-15`, and `all_outputs` carries extra plasma-header lines. They are not
caused by anything in this work, and the fix in §3.5 leaves them unchanged.

Before §3.5 was fixed, `heavy_damping` was a **fourth** difference and was
*structural*, not a rounding artefact. After the fix it is not a difference at
all: all three trees now reproduce the oracle's `heavy_damping` output
byte-for-byte, and `srccopt64/`'s count against the oracle fell from 19 to 18.
The remaining 18 are the float→double experiment and are unrelated to §3.5 —
the same 18 were measured before the fix, minus `heavy_damping`.

### 2.6 Strict-warning audit

The three trees were compiled with progressively stricter warning sets using
`tools_golden/strict_check.ps1`. Seven flags in the broadest set are not
implemented by gcc 6.3.0 and were dropped (`-Wimplicit-fallthrough`,
`-Wduplicated-branches`, `-Wdangling-else`, `-Wformat-overflow`,
`-Wformat-truncation`, `-Walloc-zero`, `-Wstringop-overflow`); the probe used to
determine this is `tools_golden/flag_probe.ps1`.

| Warning set | `srcc/` | `srccopt/` | `srccopt64/` |
|---|---|---|---|
| `-Wall -Wextra -O2` (the build flags) | 0 | 0 | 0 |
| `-Wall -Wextra -Wpedantic` | **0** | **0** | **0** |
| broad set + `-Wconversion -Wsign-conversion -Wfloat-equal …` | 238 | **82** | **72** |

The strictest set breaks down as follows.

| Category | `srcc/` | `srccopt/` | `srccopt64/` | Assessment |
|---|---|---|---|---|
| `-Wfloat-conversion` | 112 | 10 | 7 | mostly removed by the tidy-up |
| `-Wfloat-equal` | 100 | 48 | 48 | inherent — Fortran compares reals with `==` |
| `-Wdouble-promotion` | 26 | 26 | 19 | inherent — deliberate `float`/`double` mixing |
| `-Wswitch-default` | 2 | 2 | 2 | deliberate: Fortran `SELECT` had no default |
| `-Wshadow`, `-Wmissing-prototypes` | 4 | 0 | 0 | fixed in the tidy-up |

The residual warnings in `srccopt/` are therefore not sloppiness but the
arithmetic the translation is required to perform: emulating Fortran
single-precision rounding (`typin.c`), multiplying complex literals by `I`
(`xsi.c`), and Fortran's habit of testing reals for exact equality. Silencing
them would mean changing behaviour, which is the opposite of what a
byte-exactness-preserving refactor is for. `-Wstrict-overflow=5` additionally
emits six optimistic integer-range warnings in `typin.c` and was excluded as
noise.

---

## 3. Findings

### 3.1 Single precision is load-bearing — a sweep silently loses its endpoint

This is the most significant result. For the sweep `p0,.05,.01z.0022f.1`:

| Build | Points emitted | Last p |
|---|---|---|
| `srcc/` (gfortran) | 5 | 0.0400000 |
| `srccopt/` | 5 | 0.0400000 |
| `srccopt64/` | **6** | **0.0500000** |

The original Fortran **drops the requested endpoint p = 0.05**. The mechanism
was confirmed by instrumenting the range test in both trees:

```
srccopt    PM2 = 0.049999997019767761   (0.05 rounded to float)
srccopt64  PM2 = 0.050000000000000003   (0.05 as double)
```

The loop keeps advancing while `PLG <= PM2`, where `PLG` is the running
accumulator. The instrumented values explain everything:

| | bound `PM2` | `PLG` after 5 steps of `PM3` | `PLG <= PM2` |
|---|---|---|---|
| `srccopt` | `0.049999997019767761` | `0.049999998882412910` | **false** — endpoint rejected |
| `srccopt64` | `0.050000000000000003` | `0.050000000000000003` | **true** (exactly equal) — endpoint kept |

In `srccopt` the float-rounded step `0.0099999997764825821` accumulates to a
value that slightly *overshoots* the float-rounded bound `0.049999997019767761`,
so the sixth point fails the range test and the sweep stops at p = 0.04. In
`srccopt64` the double accumulation lands exactly on the double bound, `<=`
succeeds, and p = 0.05 is emitted.

Note which answer is *mathematically* intended: a sweep from 0 to 0.05 in steps
of 0.01 has six points, so `srccopt64` is arguably the more correct of the two —
the original Fortran silently loses a requested grid point to single-precision
rounding. That is precisely the dilemma this experiment exposes: `srccopt64` is
closer to the user's intent but **further from the reference program**, and
matching the reference is what makes the translation verifiable.

Either way, endpoint inclusion is decided by rounding noise rather than by
design. The important conclusion is that **widening `float` to `double` changes
program behaviour, not just its last digit** — which is why `srccopt/` keeps the
single-precision emulation and why `srccopt64/` is labelled an experiment rather
than an improvement.

### 3.2 The Fortran source disagrees with itself about π

| File | Value | Digits |
|---|---|---|
| `src/av.f90:36` | `3.1415926535897` | 14 |
| `src/whamp_engine.f90:22`, `rasy.f90:8`, `rint.f90:51`, `output.f90:44` | `3.14159265358979` | 15 |

Per instruction, `constants.h` keeps the longer value everywhere. This is
provably output-safe: the 14-digit π only ever reached the `POYN[]` Poynting
estimate in `av.c`, which the original computes and then discards (`(void)POYN`).

### 3.3 …and about the speed of light — and why `enden=` is meaningless

| File | Value |
|---|---|
| `src/av.f90:46,54,165` | `299.79`, `2.9979e8` |
| `src/energy.f90:66` | `299.792458` |

`energy.c` scales the wave magnetic field by `V = 1/C_LIGHT`, and `av.c` then
forms the field-energy ratio as `|C_LIGHT·B|²/|E|²`. Those two uses of `c` must
agree; the Fortran's disagreement left that ratio wrong by 1.6 × 10⁻⁵. Using one
value fixes it, so `constants.h` keeps the exact SI value `299.792458`.

**But the visible effect is large, for a bad reason.** `av.c` computes

```c
ENDEN = creal(ERG * conjg(ERG) / (ERG + conjg(ERG))) * 2.0 - ENFIELD;
```

where `ENFIELD = 1 + |B|²/|E|² · c²`. The two terms are nearly equal, so `ENDEN`
retains only the leading digits of their difference. Measured on `p0z.0022f.1`
with output `pzyu`:

| | term | `ENFIELD` | printed `enden=` |
|---|---|---|---|
| `c = 299.79` | 320862.74203051405 | 320847.19842531346 | `0.155E+02` |
| `c = 299.792458` | 320862.74203051405 | 320852.45972948492 | `0.103E+02` |

A 1.6 × 10⁻⁵ relative change in `c²` becomes a **34 % change** in the printed
value. This is textbook catastrophic cancellation: the 8th significant digit of
`ENFIELD` becomes the 1st significant digit of `ENDEN`. Neither printed value is
trustworthy — the original Fortran already destroys 4–5 significant digits there.
The exact `c` is used because it is *correct*, not because it is *stable*.
Fixing it properly means reformulating the energy partition, which is out of
scope.

### 3.4 Coverage gap found during validation

No golden session exercises the `enden=` / `enfl_p=` / `enfl_z=` output, and the
`all_outputs` case (which requests every output letter) prints only the help
text, because its output string is rejected before use. The energy-output code
path in `av.c` is therefore **untested in both trees**. Any future change to
`av.c` should add a case that captures those fields at full precision.

### 3.5 A real mistranslation: complex `COS`/`SIN`/`EXP` lost their imaginary part

Raising the warning level to `-Wconversion` produced five diagnostics of the
form

```
warning: conversion to 'double' from 'cd {aka _Complex double}'
         discards imaginary component [-Wconversion]
```

pointing at three statements that exist identically in `srcc/`, `srccopt/` and
`srccopt64/`:

| File | Statement | Fortran original |
|---|---|---|
| `rasy.c` | `COT = cos(PI * Y) / sin(PI * Y);` | `rasy.f90:11` `COT = COS(PI*Y)/SIN(PI*Y)` |
| `rint.c` | `COT = cos(PI * Y) / sin(PI * Y);` | `rint.f90:65` `COT = COS(PI*Y)/SIN(PI*Y)` |
| `rint.c` | `Pv  = exp(AL * C - Y * Xloc);` | `rint.f90:96` `P = EXP(AL*C - Y*X)` |

In every case the argument is complex (`Y` is declared `complex(kind=d2p)` in
Fortran and `cd` in C), so Fortran's elemental `COS`/`SIN`/`EXP` return a
complex value. The C code intended the same thing.

**Why it happened.** C99 specifies that `cos`, `sin` and `exp` are type-generic
macros that select the complex overload when given a complex argument. MinGW-w64
gcc 6.3.0 does not implement that: `<math.h>` declares only the `double`
versions, so a `double complex` argument is silently converted to `double` and
the imaginary part is thrown away. Verified directly on this toolchain for
`z = 0.3 + 0.4i`:

```
cos(z)   ->  0.72654252800536245 + 0i
ccos(z)  ->  0.14581924981527741 - 0.94020085815400634i
```

The author of the translation knew the rule — `rint.c` already uses `cexp`
correctly two lines above the defective `exp` — so this is an oversight in three
places, not a systematic misunderstanding.

**Why it was never caught.** `RYLA` routes to `RASY`/`RINT` only when
`3*(AL-10) > AY && AY*AY < 15*AL`, i.e. for strongly damped waves. The original
verification battery (`fortran_ref/battery2.ps1`) uses `p0z.0022f.1` for
essentially every case, which is weakly damped and never reaches those
routines.

**Evidence that it changes results.** Against the authoritative oracle
`fortran_ref/whamp_ref.exe` (built `gfortran -O2` per [`REPORT.md`](REPORT.md)
§8), the `heavy_damping` session (`p2z.5f1.5` / `pzf`) behaved differently:

| Build | Result |
|---|---|
| Fortran oracle | no convergence, `I= 51 IRK= 1`, `X= 0.14E+05 -0.13E+03` |
| `srcc` / `srccopt` before the fix | `Local minima!`, converges at `I= 7 IRK= 8`, `X= 0.23E+01 -0.40E+00` |
| `srcc` / `srccopt` after the fix | **byte-identical to the oracle** |

The unfixed code reported a spurious converged root in a region the reference
program never reaches.

**The fix.** Replace the three statements with the explicit complex functions
`ccos`, `csin` and `cexp`. An alternative — adding `#include <tgmath.h>` — was
also built and produces byte-identical output on all 30 sessions; the explicit
`c`-prefixed forms were chosen instead because `tgmath` would also silently
retype unrelated calls in future edits (e.g. `sqrt(float)` → `sqrtf`).

The fix was applied to **all three trees**, including the tracked `srcc/`, since
the defect originated there. Consequently `tools_golden/golden/` and
`tools_golden/golden64/` were regenerated for the single affected case
(`heavy_damping`). The pre-fix baseline output is the middle row of the table
above; no other case changed in either directory.

**Completeness of the audit.** Every `cos`, `sin`, `tan`, `exp`, `log`, `log10`,
`sqrt`, `sinh`, `cosh`, `tanh`, `pow`, `atan` and `fabs` call in all 18
translation units was checked against the declared type of its argument. The
three statements above are the only ones receiving a complex argument; the
remaining calls are genuinely real (for instance `energy.c` `Q = sqrt(Q)`
operates on a `double` produced by `creal()`). After the fix,
`-Wconversion ... discards imaginary component` reports **zero** occurrences in
all three trees.

---

## 4. Recommendation

* Use **`srccopt/`** whenever output must match the original WHAMP or the
  published documentation — it is byte-exact with `srcc/`, including
  `doc/test.txt`, and (with §3.5 applied to both) agrees with the Fortran oracle
  on every golden session except the three pre-existing ULP-level residuals.
* **`srcc/` itself was changed.** The complex-trig defect originated there, so
  the fix was applied to the tracked reference tree as well. This is the one
  behavioural change in a tracked file made by this work; it is a bug fix, not a
  refactor, and it is the only case where `srcc/` output differs from its
  pre-fix state.
* Treat **`srccopt64/`** as an experiment. It is clean, self-consistent and
  warning-free, but it is *not* closer to the reference behaviour: it emits a
  sweep point the original drops, and it disagrees with `doc/test.txt` in the
  eighth decimal. Its value is that building it exposed findings 3.1–3.4, and
  that raising the warning level for the audit exposed 3.5.
* **Enable `-Wconversion` on any future port of this code to a new toolchain.**
  The defect in §3.5 was invisible at `-Wall -Wextra` and is invisible on any
  platform where `cos` *is* type-generic. It is a portability trap that will
  reappear the moment someone rewrites this translation by hand.
* The constants extraction is the part worth keeping regardless of which tree
  is used going forward. Back-porting it to `srccopt/` is **not** a
  find-and-replace: `srccopt/` must retain `F32()` semantics and the original
  two π values and two `c` values to stay byte-exact, so the header there would
  need distinct names (e.g. `PI_AV` versus `PI`) and a different set of
  documented divergences.

---

## 5. Reproducing

```powershell
# Build (MinGW gcc 6.3.0 on PATH)
cd srcc      ; mingw32-make all      # tracked reference tree
cd ..\srccopt   ; mingw32-make all
cd ..\srccopt64 ; mingw32-make all

# Golden suites (exit code 1 on any difference)
powershell -ExecutionPolicy Bypass -File tools_golden\run_cases.ps1 `
    -Exe srccopt\whamp.exe              -Cases tools_golden\cases.json `
    -OutDir tools_golden\actual  -CompareDir tools_golden\golden
powershell -ExecutionPolicy Bypass -File tools_golden\run_cases.ps1 `
    -Exe srccopt\whamp_engine_test.exe  -Cases tools_golden\cases.json `
    -OutDir tools_golden\actual  -CompareDir tools_golden\golden
powershell -ExecutionPolicy Bypass -File tools_golden\run_cases.ps1 `
    -Exe srccopt64\whamp.exe            -Cases tools_golden\cases.json `
    -OutDir tools_golden\actual64 -CompareDir tools_golden\golden64
powershell -ExecutionPolicy Bypass -File tools_golden\run_cases.ps1 `
    -Exe srccopt64\whamp_engine_test.exe -Cases tools_golden\cases.json `
    -OutDir tools_golden\actual64 -CompareDir tools_golden\golden64
```

WHAMP's REPL never exits on end-of-input, so every scripted session must end
with `s` at the prompt; `cases.json` already does this.

The `clean` targets in all three `Makefile`s invoke `rm -f` and therefore do not
work on Windows; delete `*.o` and `*.exe` manually.

To reproduce the oracle comparison of §2.5, capture the reference output first
and then diff each tree against it:

```powershell
# The oracle must be built with -O2 (see REPORT.md section 8); the
# src/Makefile default omits it and yields a different binary.
powershell -ExecutionPolicy Bypass -File tools_golden\run_cases.ps1 `
    -Exe fortran_ref\whamp_ref.exe -Cases tools_golden\cases.json `
    -OutDir $env:TEMP\out_ref
```

The strict-warning audit of §2.6 is driven by two scripts kept alongside the
harness:

```powershell
# Which -W flags this gcc actually implements
powershell -ExecutionPolicy Bypass -File tools_golden\flag_probe.ps1

# Compile every .c in a tree with an explicit flag set and tally diagnostics
powershell -ExecutionPolicy Bypass -File tools_golden\strict_check.ps1 `
    -Src srccopt -Flags '-Wall -Wextra -Wpedantic'
```

---

## 6. Size

| Tree | `.c` files | `.c` lines | `.h` files | `.h` lines |
|---|---|---|---|---|
| `srcc/` (reference) | 18 | 3 000 | 4 | 386 |
| `srccopt/` | 18 | 3 096 | 5 | 264 |
| `srccopt64/` | 18 | 3 066 | 6 | 438 |

Counts include the four explanatory comment lines added with the §3.5 fix in
each tree.

`srccopt/` has fewer header lines than `srcc/` because `whamp.h` was slimmed
and the edit-descriptor code moved into a `.c` file. `srccopt64/` adds
`constants.h` (167 lines, mostly documentation of the findings above).

Toolchain for all measurements: `gcc` / `gfortran` 6.3.0 (`C:\MinGW\bin`),
Windows 10/11 x64.

---

## 7. Conclusion

Two things were achieved, and one assumption was disproved.

**The tidy-up is safe.** `srccopt/` is a genuine refactor: dead code removed,
the Fortran edit-descriptor emulation isolated in its own module, the central
header slimmed, allocation and argument-parsing boilerplate collapsed into named
helpers. None of it changed behaviour. All 30 golden sessions are byte-identical
to `srcc/`, `doc/test.txt` reproduces exactly, and the build is warning-free. It
can be dropped in wherever `srcc/` is used today, including behind the MATLAB
front end, with no change to any downstream consumer.

**The precision experiment is not.** Widening every `float` to `double` was
expected to be a purely numerical change — the same program, with better-rounded
last digits. It is not. Finding 3.1 shows that single precision participates in
control flow: the original's float-rounded sweep bound silently rejects a
requested grid point, and `srccopt64/` accepts it. That is a behavioural
difference, invisible to any test that does not specifically count sweep points,
and it is the reason `srccopt64/` is presented as an experiment rather than as an
improvement. Higher precision and higher fidelity to the reference program turned
out to be in tension, and for a translation whose entire value proposition is
fidelity, fidelity wins.

**What the experiment was nonetheless worth.** It was the only way to find out
how much of the original's behaviour is rounding rather than arithmetic, and it
produced three durable results: the source disagrees with itself about π and
about `c` (§3.2, §3.3); the `enden=` energy output is destroyed by catastrophic
cancellation and should not be quoted to more than one significant digit, if at
all (§3.3); and the energy code path is not covered by any test in either tree
(§3.4). These are properties of the *original* WHAMP, uncovered here for the
first time, and they stand independently of which tree is adopted.

**The most valuable result was an accident.** None of the above was the point of
raising the warning level, but compiling with `-Wconversion` to audit strictness
turned up §3.5 — a genuine mistranslation in which three complex
trigonometric expressions silently lost their imaginary part on a toolchain that
does not implement type-generic `math.h`. It had survived the original
byte-for-byte verification, survived the construction of a 25-session golden
suite, and survived review of the translation, because no test in the suite
reached the affected code path and the default warning level could not see it.
The code reported a confidently-wrong converged root where the reference program
reports non-convergence. Fixing it moved all three trees *closer* to the Fortran
oracle, which is the one direction a fidelity-preserving refactor is not supposed
to move.

That is the argument for treating warnings as part of the verification, not as
housekeeping: `-Wall -Wextra` was clean in all three trees before the fix and is
clean after it, and it was blind to the only real bug found in the whole exercise.

**What is claimed, and what is not.** Claimed: all three trees compile
warning-free under the stated toolchain at `-Wall -Wextra -O2` and at
`-Wall -Wextra -Wpedantic`; each reproduces its own golden baseline exactly;
`srccopt/` reproduces `srcc/` byte-for-byte; the documented acceptance test
passes; and all three now agree with the Fortran oracle on `heavy_damping`, the
case affected by §3.5. Not claimed: that `srccopt64/` matches the original (it
does not, in two specific and documented ways), that the golden suite is complete
(it is not — §3.4, and §3.5 is direct evidence of that), that the three residual
oracle differences in `srcc/` and `srccopt/` are understood beyond being
ULP-level and pre-existing, or that the numerical findings apply to WHAMP builds
other than this Fortran-90-derived one.

**Suggested next steps**, in rough order of value:

1. Close the coverage gap of §3.4 by adding a golden case that captures
   `enden=`, `enfl_p=` and `enfl_z=` at full precision. This is a prerequisite
   for touching `av.c` again.
2. Decide the sweep-endpoint question of §3.1 explicitly rather than by
   rounding accident — either document that WHAMP drops a requested endpoint, or
   fix the range test (e.g. compare against `PM[2] - PM3/2`) in a deliberate,
   tested change to both trees.
3. If `constants.h` is back-ported to `srccopt/`, give the divergent constants
   distinct names (`PI_AV`, `C_LIGHT_AV`) so the two original spellings remain
   visible and byte-exactness is preserved.
4. Add `-Wconversion` to the regular build, or at least to a scheduled check,
   and treat `discards imaginary component` as an error. It is the only warning
   class that caught a real defect here, and `heavy_damping` now guards the
   regression but only if someone runs it.
5. Leave `srccopt64/` in the tree as a labelled experiment, or delete it. Its
   diagnostic value has already been extracted into this report; keeping a
   second, subtly-different copy of the whole program around is a long-term
   maintenance cost.

The broader lesson is worth stating plainly: in legacy numerical code, rounding
behaviour *is* behaviour. A refactor that preserves results must preserve
precision, and an "obvious modernisation" such as replacing `float` with
`double` is not a cleanup — it is a change of program, to be justified by
evidence in the same way as any other.
