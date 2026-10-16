# Branch Tracking for the Pure MATLAB WHAMP Solver (`srcm/whamp_track.m`)

AI-assisted work by Walter Puccio 2026, extending
[`REPORTmatlab.md`](REPORTmatlab.md) (which documents the `srcm/whamp.m`
refactor and its byte-for-byte regression harness).

> **Note on names.** The solver this report calls `whamp_track.m` is now
> `srcm/whamp_track_test.m` (the modular build, kept as the reference); the
> plain name `srcm/whamp_track.m` now belongs to its self-contained build, and
> [`whamp_track_v2.m`](srcm/whamp_track_v2.m) supersedes both. See the rename
> table in [`REPORTtrack_v2.md`](REPORTtrack_v2.md). Everything below is
> unchanged and describes the Version 1 work as it was done.

## Abstract

WHAMP's root finder loses the branch it is tracking: the frequency it reports
at successive points of a `k` sweep can jump onto a different, unrelated mode.
This report establishes the cause, documents two new solvers that fix it, and
measures the result.

The cause is not a numerical defect. Multiple roots genuinely coexist at a
single `(k⊥, k∥)`, and `whamp.m` carries one complex seed from point to point
with no test that the root it lands on is the one it was following. Because
`D(ω)` is meromorphic — one pole per species gyroharmonic — a Newton step
started near a pole lands in whichever basin it falls into. The jump is correct
Newton behaviour; what is missing is branch identity.

Two solvers were added, both driving physics extracted into a shared
`+whampkit` package so that every version evaluates a bit-identical dispersion
function. `whamp.m` itself is **unmodified**, and its 30-case `%.20g` regression
still passes byte-for-byte.

* **`whamp_track.m`** — local continuation with a trust region, a continuity
  check, and multi-start re-anchoring. Over nine adversarial sweep lines
  measured against the baseline on a step eight times finer:
  **13 branch hops → 0**, at ~5–10× the baseline cost.
* **`whamp_global.m`** — enumerate all roots per point, match them across
  points by optimal assignment. Structurally more appealing, and **worse on
  both quality measures at ~800× the cost**. Kept, documented as not
  recommended, as evidence that the global route is the wrong trade-off here.

Four bugs were found on the way, and the most instructive were in the
measurement and the diagnostics rather than in the solver: a trust-region Newton
that **falsely converged against the wall of its own disc**, manufacturing roots
with `|D|` of order 10⁸–10¹⁰ (§5); and a comparison harness that **silently
skipped points the reference had not reached**, which understated the baseline's
hop count (§6). Two further fidelity gaps in `whamp_track` — diagnostics that
excluded exactly the points that worked hardest, and a lost "too heavily damped"
flag — are in §7.

---

## 1. The symptom

`whamp.m` sweeps `kpar` at fixed `kperp`, carrying a predictor from the previous
point. On a drift-unstable plasma the reported frequency does not move
continuously along a mode — it leaps between distant values, and then rides the
new branch for the rest of the line.

Concretely, from `compareExample.m` line 5 (the anisotropic plasma at
`kperp = 0.01`), comparing the coarse run against the same code on a step eight
times finer:

| `kpar` | reference branch | `whamp.m` reports |
|---|---|---|
| 0.016 | 0.9565857 | 14.03655 |
| 0.024 | 0.9677542 | 28.1007 |
| 0.032 | 0.970421 | 47.32616 |
| 0.040 | 0.9709409 | 71.41516 |

The reference branch is nearly flat; the coarse run has jumped onto a steeply
rising mode at the very first step and never returns. Nothing in the output flags
this — every one of those values is a genuine, well-converged root of the
dispersion relation.

## 2. Diagnosis

Three facts, established before any solver was written.

**All reported points are genuine roots.** An earlier probe suggested `whamp.m`
was reporting non-roots. That was wrong — an artifact of evaluating candidates at
the wrong `z`. A direct check confirmed every point the baseline reports is a
zero of `D`. The sole failure mode is branch selection.

**Multiple roots genuinely coexist.** At a given `(k⊥, k∥, z)` the dispersion
determinant has many zeros. There is no such thing as "the" root, so a solver
that does not decide which branch it wants will have one chosen for it by the
starting value.

**`D(ω)` is meromorphic**, with a pole at each species gyroharmonic. This rules
out the obvious tools: `roots()` needs a polynomial, and `fzero` is real-valued
only. It also explains the mechanism — an unconstrained Newton step started near
a gyroharmonic pole is thrown arbitrarily far.

The conclusion is that the fix must be a *policy* layered on top of the existing
iteration, not a better iteration. `whamp.m`'s Newton is fine.

## 3. Approach: extract the physics, add policy on top

The physics core of `whamp.m` was mechanically extracted into a package,
`srcm/+whampkit/` (36 files, 1369 lines), so that the baseline and both new
solvers evaluate exactly the same dispersion function, susceptibility, wave
energy and sweep predictor. Nothing in `whamp.m` was changed — it remains the
reference, and its existing regression harness
([`srcm/whamp_dump.m`](srcm/whamp_dump.m), [`srcm/baseline.txt`](srcm/baseline.txt))
is the guard that the extraction was truly behaviour-preserving.

| File | Lines | Role |
|---|---|---|
| [`srcm/whamp_track.m`](srcm/whamp_track.m) | 384 | Version 1 — local continuation (recommended) |
| [`srcm/whamp_global.m`](srcm/whamp_global.m) | 507 | Version 2 — global enumeration (experimental) |
| [`srcm/compareSolvers.m`](srcm/compareSolvers.m) | 145 | Measurement harness |
| [`srcm/compareExample.m`](srcm/compareExample.m) | 36 | Nine-line driver that reproduces §8 |
| [`srcm/+whampkit/`](srcm/+whampkit/) | 1369 | Shared physics, 36 files |

## 4. `whamp_track.m` — the recommended solver

Drop-in for `whamp`: same signature, same `Output` struct, plus
`Output.flagSuspectBranch` and `Output.diagnostics`. Three mechanisms:

1. **Trust-region Newton.** Every retry is confined to a disc centred on the
   *previous root on this branch*, with radius `trustFactor ×` the branch's own
   recent step size. A Newton step can no longer be thrown into a distant basin.
2. **Continuity check.** A converged root is compared against the previous root
   on the same sweep line. If it moved more than `jumpFactor ×` the branch step,
   it is *rejected* rather than recorded, and the point is retried.
3. **Multi-start re-anchoring.** Seeds are tried in priority order: the engine's
   own extrapolation, then the previous root, then offsets of the previous root
   scaled by the branch step. The first candidate that converges **and** passes
   continuity wins.

### 4.1 The design decision that matters

> Every scale is measured on the **branch** — the roots the sweep has already
> accepted — never on the extrapolated seed.

The seed is produced by extrapolating the previous point. When it is bad, it is
bad by roughly the size of the jump it is about to cause. Any trust radius or
rejection threshold expressed in terms of the seed is therefore *circular*: it
admits precisely the failure it is meant to catch. The branch's own recent step
size is independent of the seed, so it still says something useful after the seed
has already gone wrong.

This was the single change that took the solver from ineffective to zero hops.

### 4.2 Why it costs nothing when there is no problem

The **first seed is tried with `radius = Inf`** — it *is* `whamp`'s solve, from
the same extrapolated seed, unconstrained. On a well-behaved sweep line the
solver lands on the same root in the same number of iterations and the machinery
never engages. Lines 7–9 of the test set are identical to the baseline to
~10⁻¹⁰, and cost 0.6–0.8 s against the baseline's 0.1 s (the residual overhead is
the continuity bookkeeping, not extra Newton work).

### 4.3 History must span the whole sweep

The tracking history is kept in **visit order across the entire sweep**,
mirroring the engine's own `XO`/`PO`/`ZO` predictor, which also survives line
boundaries. Resetting it at each line makes the continuity check vacuous:
`advance_sweep` treats every `kpar` point as a new outer line when there is a
single `kperp`, so a per-line history is empty at every point and nothing is ever
rejected. This is why the first version of `whamp_track` appeared to do nothing
at all.

### 4.4 Options

| Option | Meaning | Default |
|---|---|---|
| `.trustFactor` | Newton disc radius, in branch steps | 3 |
| `.jumpFactor` | Reject a root further than this from the branch | 6 |
| `.relFloor` | Floor on the branch step, as a fraction of the previous root, so a flat line does not reject round-off | 0.25 |
| `.scanOffsets` | Last-resort seed offsets, in branch steps | `[0.5 -0.5 1 -1 2 -2]*0.5` |
| `.useScan` | Try the scan offsets when nothing else works | `true` |
| `.maxIterations` | Newton cap | from `InputParameters` |

`Output.flagSuspectBranch` is a per-point 0/1 array in the same orientation as
`Output.f`: a converged root that the continuity check did not like. Plot it to
see where a sweep is least trustworthy.

## 5. A real bug: trust-region Newton converging against the wall

Found in [`srcm/+whampkit/newtonRoot.m`](srcm/+whampkit/newtonRoot.m) while
investigating a single residual hop.

When the trust budget was exhausted, `limitStep` returned `correction = 0`. The
caller's convergence test was

```matlab
abs(correction) <= freqRelTol * abs(X) || abs(correction) < freqAbsTol
```

A zero correction passes trivially. The iteration therefore **"converged" while
pinned against the wall of the trust disc, with no regard for whether `D` was
small** — residuals of order 10⁸–10¹⁰ were being recorded as roots. The symptom
in the wild was a "root" appearing at *exactly* `trustFactor × scale` from the
previous root, which is how it was recognised.

The trust region was manufacturing fake roots.

Two guards, both required:

* `limitStep` now returns a third output `exhausted`; when set, the main loop
  gives up (`rootFound = false`) and lets the caller try another seed.
* Convergence additionally requires `~wasClamped`, since a truncated
  correction's magnitude says nothing about `|D|`.

The general lesson is worth recording: **any trust region must be able to report
that it bit.** A cap that silently shrinks a step to zero turns a convergence
test on the step size into a test of the cap.

The consequence is that `whamp_track` now *refuses* points the baseline used to
record. That is correct behaviour — see §8.1.

## 6. A bug in the measurement

The harness initially reported that the baseline made **zero** hops on the golden
plasma at `kperp = 0.02`. Running the baseline directly showed it hopping
identically (348.9, 544.7, 777.1, 1038.4 …).

`judge()` in [`srcm/compareSolvers.m`](srcm/compareSolvers.m) matched each coarse
point against the reference and `continue`d past any point the reference had not
converged to. On hard lines the reference dies early (21–23 of 81 points), so
every baseline hop beyond that point was invisible. A hop count computed this way
is not a count — it is a count over an undisclosed subset.

`judge()` now returns `checked` and `skipped`, and every output line prints them:

```
  whamp         hops 1  missing 1  worst dev 1.31e+02    0.2s  (3 checked, 7 unjudgeable)
```

**All hop counts in §8 are lower bounds**, and the baseline's is the one that
benefits from the slack — the true baseline error rate is worse than 13.

## 7. Two fidelity gaps in `whamp_track`

Both were found by reading the diagnostics against a direct run rather than by
any test failing, and both made the solver look less informative than it was.

### 7.1 Diagnostics excluded the points that worked hardest

`attemptsPerPoint`, `rejected` and `clampedSteps` were accumulated inside the
`if status == 1` branch — so they counted only points that **succeeded**. The
points where the seed list was worked hardest are precisely the ones that failed,
and they were being left out.

On the golden plasma at `kperp = 0.01` the per-point attempt array read
`1 1 2 1 1` (7 entries, over successes only); after moving the accumulation
outside the success branch it reads `1 1 2 8 8 8 7 8 1 1` over all 10 visited
points. The old array described the easy points and hid the interesting ones.

### 7.2 "Too heavily damped" was flattened into "did not converge"

`whamp.m` distinguishes status 2 (the susceptibility reports the point as too
heavily damped — a *physical* statement) from status 0 (the iteration gave up — a
*numerical* one), and exposes the first as `Output.flagTooHeavilyDamped`. The
multi-start loop collapsed both to 0, so the flag never fired.

The fix asks the **engine's own seed** for its status, since that is exactly the
solve `whamp` would have performed had it been on this branch. If that seed was
refused on damping grounds, the point is damped regardless of what the fallback
seeds managed. On the golden plasma at `kperp = 0.01` the flags now read:

| point | `kpar` | `whamp` | `whamp_track` |
|---|---|---|---|
| 1 | 0.000 | found 0.408227 | found 0.408227 |
| 2 | 0.008 | found 0.985457 | found 0.985457 |
| 3 | 0.016 | **damped** | found 1.000741 (matches reference) |
| 4 | 0.024 | found 296.363 | damped |
| 5 | 0.032 | found 493.420 | damped |
| 6 | 0.040 | found 729.924 | damped |
| 7 | 0.048 | found 996.551 | found 1.892077 (rescued by fallback) |
| 8 | 0.056 | found 1284.148 | found 1284.148, **flagged suspect** |
| 9 | 0.064 | found 1584.345 | found 1584.345 |
| 10 | 0.072 | found 1889.886 | found 1889.886 |
| 11 | 0.080 | — | — |

## 8. Results

Nine sweep lines: three plasmas × `kperp ∈ {0.005, 0.01, 0.02}`, coarse step
`h = 0.008`, reference = baseline at `h/8`. Reproduce with
[`srcm/compareExample.m`](srcm/compareExample.m).

| Solver | Branch hops | Missing points | Cost per line |
|---|---|---|---|
| `whamp` | **13** | 11 | ~0.1 s |
| **`whamp_track`** | **0** | 13 | 0.5–1.9 s |
| `whamp_global` | 7 | 86 | ~80–95 s |

Per line (`whamp` → `whamp_track`):

| # | Plasma | `kperp` | hops | checked | `whamp` worst dev | `whamp_track` worst dev |
|---|---|---|---|---|---|---|
| 1 | golden | 0.005 | 1 → 0 | 3 | 1.3e+02 | 5.3e-12 |
| 2 | golden | 0.01 | 0 → 0 | 2 | 3.2e-12 | 3.2e-12 |
| 3 | golden | 0.02 | 0 → 0 | 2 | 3.3e-11 | 3.3e-11 |
| 4 | aniso | 0.005 | 4 → 0 | 6 | 6.9e+01 | 2.3e-12 |
| 5 | aniso | 0.01 | 4 → 0 | 6 | 7.0e+01 | 1.6e-12 |
| 6 | aniso | 0.02 | 4 → 0 | 6 | 7.6e+01 | 3.2e-11 |
| 7 | beam | 0.005 | 0 → 0 | 10 | 9.8e-10 | 9.8e-10 |
| 8 | beam | 0.01 | 0 → 0 | 10 | 6.8e-10 | 6.8e-10 |
| 9 | beam | 0.02 | 0 → 0 | 10 | 2.6e-10 | 2.6e-10 |

Plasmas (all `Species` entries `{m, n, t, a, vd}`):

* **golden** — `B = 100`; He⁺ `{4, 1, 0.05, 5, 1}`, O⁺ `{16, 0.5, 0.01, 1, 0}`, e⁻ `{0, 1, 0.001, 1, 0}`. Strong anisotropy plus drift; this is the case that exposes the bug.
* **aniso** — `B = 30`; H `{1, 0.6, 0.05, 2, 0}`, He⁺ `{4, 0.4, 0.02, 3, 0}`, e⁻ `{0, 1, 0.01, 1, 0}`.
* **beam** — `B = 100`; H `{1, 0.7, 0.05, 1, 2.5}`, O⁺ `{16, 0.3, 0.05, 1, 0}`, e⁻ `{0, 1, 0.001, 1, 0}`. Drifting beam, well-behaved branch.

Every hop the baseline makes is on a line where the reference could judge it, and
`whamp_track` reproduces the reference to **10⁻¹²–10⁻¹¹** where the baseline is
off by **10¹–10²** — fourteen orders of magnitude. On the three well-behaved beam
lines the two are indistinguishable.

### 8.1 The coverage difference: "missing" is not a quality metric

`whamp_track` reports 13 missing points against the baseline's 11, all of the
difference on line 2. That number is misleading in both directions, and it is
worth being precise about why.

Line 2 is the golden plasma at `kperp = 0.01`, detailed in §7.2. The reference
itself converges only to `kpar = 0.022` — 23 of 81 fine points — so points 4–11
cannot be judged at all. What the two solvers actually do there:

* `whamp` reports 9 points, of which **the last 7 are on a branch that had
  already left the reference** at point 3 and rose to 1889.886.
* `whamp_track` reports 7 points, of which the first 3 match the reference
  exactly, point 7 (1.892077) is a plausible continuation of the same slowly
  rising branch, and only point 8 (1284.148) is a jump — and that one **is
  flagged**.

So `whamp_track` has *fewer* points and *more* of them are on the right branch.
Counting raw coverage rewards the solver that jumped. **Coverage counts are not
comparable between two runs that have diverged onto different branches** — from
the divergence point on they are measuring different things.

The trade is deliberate. A hole, or a flagged point, is more honest than an
unflagged wrong point; prefer completeness over purity and plot
`flagSuspectBranch`, but do not read a larger "missing" as a worse solver.

### 8.2 Why the global approach lost

`whamp_global.m` is the structurally attractive answer: at each point enumerate
*all* roots of `D` in a frequency box, order them with `esort`, and match across
points by optimal assignment. Branch identity becomes a property of the root set
rather than an accident of the starting value.

It loses for a reason that is structural rather than a matter of tuning. In cold
multi-species plasmas the susceptibility expansion fails over most of the complex
plane, so the enumeration cannot sample the box reliably. The root set it
produces is incomplete and step-size dependent, and a global matcher built on an
incomplete root set cannot beat a local continuation that never needed the set in
the first place. It is 7 hops / 86 missing at ~800× the baseline cost.

It is retained, clearly marked **experimental — not recommended**, with these
measurements in its header, both as a documented alternative and so that the next
person does not spend a week rediscovering why.

## 9. Verification

| Check | Result |
|---|---|
| `whamp.m`, `matlab/whamp.m` modified | **No** — `git diff` empty against the pre-work commit |
| 30-case `%.20g` regression (`whamp_dump`) | **Byte-identical** to `baseline.txt` (60492/60492 chars), verified repeatedly including after each of the §7 fixes |
| `checkcode` over all 40 new/changed files | **0 errors** |
| `whamp_track` default-plasma smoke test | Matches baseline `f(1) = 0.1063736104`, `flagSuspectBranch(1) = 0` |
| `whamp_global` after the `solvePoint` signature change | Works (new `center` argument is optional) |
| Nine-line comparison after the §7 fixes | Unchanged: 13 hops → 0 |

The regression harness is the load-bearing evidence that extracting the physics
into `+whampkit` changed nothing: the baseline is the code users run, and it
still produces the same 2616 lines at 20 significant digits.

## 10. Known limitations

1. **`relFloor = 0.25` is the least principled number in the design.** It is too
   loose on flat branches: the reference steps ~0.004 while `scale` floors at
   0.247, giving a jump tolerance of 1.48 — about 60× larger than it should be.
   The solution struct carries `DP`, `DZ` and `DX`, so the exact
   implicit-differentiation step `|((p − P₀)·DP + (z − Z₀)·DZ)/DX|` is available
   and would replace the median heuristic and its floor entirely. This is the
   highest-value follow-up.
2. **The reference is itself limited.** On the two hardest lines it converges for
   only 21–23 of 81 fine points, so most of those lines are unjudgeable and the
   hop counts are lower bounds. A better reference — the winding-number counter
   in [`+whampkit/windingNumber.m`](srcm/+whampkit/windingNumber.m) is the natural
   candidate, being scale-free — would let more points be scored.
3. **Cost is 5–10× the baseline** on hard lines, from the extra seed attempts.
   Acceptable for offline dispersion work; not automatically for a contour sweep
   with thousands of lines.
4. **`IERR` (too heavily damped) aborts the solve** at both `solvePoint` entry and
   inside `newtonRoot`, identically to the baseline. In diagnostics a residual of
   `Inf` means `IERR`, not a large residual — a source of confusion during this
   work.
5. **Not back-ported.** `src/whamp_engine.f90` and the C versions share the same
   single-seed structure and have the same defect.
6. **Untested configurations**: `varyKzFirst = 1`, genuine 2-D sweeps,
   relativistic beams, and all ten species slots.

## 11. Reproducing

```matlab
addpath('srcm');            % -batch does not put the CWD on the path
compareExample;             % the nine-line table from section 8
```

Or for a single line, with the reference at a finer step:

```matlab
L = struct('plasma', myPlasma, 'kperp', 0.01, 'z1', 0.08);
opts = struct('h', 0.008, 'refineFactor', 8, 'fstart', 0.4, ...
              'solvers', {{@whamp, @whamp_track}}, 'verbose', true);
compareSolvers(L, opts);
```

Note the double cell around `@whamp`: `struct('solvers', {a, b})` builds a 1×2
struct array, and `opts.solvers` then expands to two comma-separated arguments.

To inspect what the machinery did at each point, use
`Output.diagnostics.perPoint` — one struct per visited point, with the seed tried,
its status, the resulting root and the Newton displacement. Note that `tried` is a
cell array nested one level inside the packed struct.

And to confirm the baseline is untouched:

```matlab
whamp_dump('now.txt');      % then byte-diff against srcm/baseline.txt
```

## 12. Conclusion

**Root jumping in WHAMP is not a numerical failure and cannot be fixed by a
better root finder.** Multiple roots coexist at every point of the sweep, the
dispersion function is meromorphic with a gyroharmonic pole at every species
harmonic, and `whamp.m`'s single carried seed simply does not decide which branch
it is on. Newton does the rest. The fix is a policy layer, not a new iteration.

**`whamp_track.m` eliminates the problem on every line tested: 13 branch hops
reduced to 0**, with agreement to 10⁻¹² against a refined reference where the
baseline is off by 10¹–10², at 5–10× the baseline cost, and with *zero* change to
`whamp.m` — whose bit-for-bit regression still passes. On well-behaved lines it is
indistinguishable from the baseline, because the first seed it tries is the
baseline's own solve, unconstrained.

**The one idea that made it work was measuring everything against the branch
rather than against the seed.** The earlier, ineffective version sized its trust
region from the extrapolated seed — which is circular, because a bad seed is bad
by exactly the size of the jump it causes. Switching to the branch's own
historical step size took the solver from no effect to zero hops.

**Three of the four bugs were not in the solver's core.** A trust region that
silently reduces a step to zero turns a step-size convergence test into a test of
the cap, and reports the wall of its own disc as a root. A comparison harness that
quietly skips unjudgeable points makes the baseline look perfect on the hardest
line. Diagnostics that accumulate only on success describe the easy cases and hide
the interesting ones. All three produce plausible output and pass every test you
would think to write; each was caught only by checking the measurement against a
direct run, and by reporting the denominator alongside the count.

**The "missing points" metric should not be trusted on its own.** `whamp_track`
reports 13 missing where the baseline reports 11, and on the one line where they
differ the baseline's extra points are all on a branch it hopped onto and never
left. A coverage count that rewards a solver for jumping is not a quality metric;
`flagSuspectBranch` is the honest version of the same information.

**What is claimed:** over nine sweep lines spanning three multi-species plasmas
and three perpendicular wavenumbers, `whamp_track` records no branch hops where
`whamp.m` records thirteen, while `whamp.m` remains byte-identical to its
pre-work state.

**What is not claimed:** that zero hops holds generally — the hop counts are lower
bounds, limited by how far the refined reference itself converges, and only
`varyKzFirst = 0` and one-dimensional sweeps were tested. That the coverage counts
are comparable between solvers that have diverged onto different branches (they
are not). That `relFloor` is right; it is a heuristic, and it is the obvious thing
to replace. And that the Fortran and C engines are fixed, because they are not.

**Suggested next steps**, in rough order of value:

1. Replace the median step-scale and `relFloor` with the exact
   implicit-differentiation step from `DP`/`DZ`/`DX`, removing the last
   hand-tuned constant.
2. Strengthen the reference so more points are judgeable — a winding-number based
   ground truth would remove the lower-bound caveat from the headline result.
3. Back-port the tracking to `src/whamp_engine.f90` and the C engines, which carry
   the same defect.
4. Test `varyKzFirst = 1` and genuine two-dimensional sweeps, where the
   visit-order history assumption needs re-examining.
5. Add a `flagSuspectBranch`-aware plotting helper so users can see branch
   confidence directly on a dispersion diagram.
6. Consider whether `whamp_global.m` should stay in the tree at all, now that its
   header records why it should not be used.
