function Output = whamp_track_test(PlasmaModel, InputParameters, Options)
%WHAMP_TRACK_TEST Branch tracking, modular reference build (Version 1).
%   This is the original modular Version 1: it calls the +whampkit package
%   rather than inlining it, and is kept so that whamp_track.m can be checked
%   against it bit-for-bit.  For normal use call WHAMP_TRACK (self-contained
%   Version 1) or WHAMP_TRACK_V2 (no hand-tuned scales).
%   OUTPUT = WHAMP_TRACK_TEST(PlasmaModel, InputParameters) walks the same
%   two-axis k sweep as WHAMP, but adds three mechanisms that stop the
%   reported root from hopping between coexisting branches:
%
%     1. TRUST-REGION NEWTON.  Every retry is confined to a disc centred on
%        the *previous root on this branch*, sized from the branch's own
%        recent step size.  Because D(w) is meromorphic -- one pole per
%        species gyroharmonic -- an unconstrained Newton step started near a
%        pole can land in a distant branch's basin.  The cap makes that
%        impossible.  The engine's own extrapolated seed is still tried
%        first and unconstrained, so a well-behaved line reproduces WHAMP
%        exactly and pays nothing.
%
%     2. CONTINUITY CHECK.  A converged root is compared against the
%        previous root on the same sweep line, measured on a robust scale
%        built from that line's own recent step sizes.  A root that moves
%        implausibly far is not silently recorded; it is rejected and the
%        solver retries from other seeds.
%
%     3. MULTI-START RE-ANCHORING.  When a seed fails or is rejected, the
%        point is retried from a short priority list of alternative seeds
%        (previous root, then offsets of it scaled by the branch step).  The
%        first candidate that converges *and* passes continuity wins.
%
%   The tracking history is kept in *visit order across the whole sweep*,
%   mirroring the engine's own XO/PO/ZO predictor, which also survives line
%   boundaries.  Resetting it at each line would leave the continuity check
%   with no previous root and make it vacuous.
%
%   OUTPUT = WHAMP_TRACK_TEST(...) with no arguments uses the default plasma.
%
%   Options (struct, all optional):
%   Every scale below is measured on the *branch* -- the roots the sweep has
%   already accepted -- never on the extrapolated seed.  That distinction is
%   the whole design: the seed is produced by extrapolating the previous
%   point, so when it is bad it is bad by about the size of the jump it is
%   about to cause, and any threshold expressed in terms of it is circular.
%
%     .trustFactor    Newton disc radius = trustFactor * branch step.
%                                                          Default 3
%     .jumpFactor     Reject a root > jumpFactor * branch step from
%                     the previous root.                   Default 6
%     .relFloor       Floor on the branch step, as a fraction of the
%                     previous root, so a flat line does not reject
%                     round-off.                           Default 0.25
%     .scanOffsets    Offsets (multiples of the branch step) applied to the
%                     previous root as last-resort seeds.
%                                                Default [0.5 -0.5 1 -1 2 -2]*0.5
%     .useScan        Try the scan offsets when nothing else works.  Default true
%     .maxIterations  Newton iteration cap.  Default from InputParameters
%
%   OUTPUT.diagnostics reports what the machinery did, so this solver can be
%   scored against WHAMP_GLOBAL:
%     .attemptsPerPoint   seeds tried per point
%     .clampedSteps       Newton steps cut by the trust region
%     .rejected           roots rejected by the continuity check
%     .rescued            points saved by a fallback seed
%     .flagged            points recorded despite failing continuity
%     .perPoint           array of structs, one per visited point
%
%   Output.flagSuspectBranch is the per-point 0/1 array (same orientation as
%   Output.f).  A one means the point is a converged root of the dispersion
%   relation, but one the continuity check did not like -- either it moved too
%   far for this branch, or it is all that was left after every seed was
%   exhausted.  Plot it to see where the sweep is least trustworthy.
%
%   Known limitation
%   ----------------
%   Refusing to record a discontinuous root can cost coverage.  Over the nine
%   lines in compareExample.m this solver reports 13 points where WHAMP
%   reports 11, and all of the difference is on one line: there the baseline
%   jumps onto a distant branch at kpar=0.016 and then rides it happily for
%   the remaining points, while this solver refuses the jump, keeps trying,
%   and gives out sooner.  That is the intended trade -- a hole is more
%   honest than a wrong point -- but it means the coverage counts are not
%   directly comparable once the two have diverged onto different branches,
%   and it is why Output.flagSuspectBranch exists: prefer a flagged point to
%   a hole if completeness matters more than purity.
%
%   The physics is shared with WHAMP through the whampkit package, so every
%   version evaluates exactly the same dispersion function.  WHAMP itself is
%   left untouched and its bit-for-bit regression baseline still applies.

d = fileparts(mfilename('fullpath'));
addpath(fileparts(d));

if nargin < 1 || isempty(PlasmaModel)
    PlasmaModel = whampkit.complete_plasma(whampkit.default_plasma());
else
    PlasmaModel = whampkit.complete_plasma(PlasmaModel);
end
if nargin < 2 || isempty(InputParameters)
    InputParameters = whampkit.complete_parameters(whampkit.default_parameters());
else
    InputParameters = whampkit.complete_parameters(InputParameters);
end
if nargin < 3 || isempty(Options)
    Options = struct();
end

opt = local_options(Options, InputParameters);

[mass, density, temperature, anisotropy, delta, drift] = ...
    whampkit.species_arrays(PlasmaModel.Species);
plasma = whampkit.prepare_plasma(mass, density, temperature, anisotropy, ...
    delta, drift, 0.0279928 * PlasmaModel.B);

kpSpec = whampkit.axis_spec(InputParameters.kperp, 'kperp');
kzSpec = whampkit.axis_spec(InputParameters.kpar, 'kpar');
kp = whampkit.fill_axis(kpSpec);
kz = whampkit.fill_axis(kzSpec);
kpPhysical = whampkit.physical_axis_value(kp, InputParameters.useLog);
kzPhysical = whampkit.physical_axis_value(kz, InputParameters.useLog);

shape = [numel(kp), numel(kz)];
values = whampkit.emptyValues(shape);
flagFound = zeros(shape, 'int32');
flagDamped = zeros(shape, 'int32');
flagNoConvergence = zeros(shape, 'int32');
flagSuspect = zeros(shape, 'int32');

cycleZFirst = InputParameters.varyKzFirst ~= 0;
sweep = whampkit.sweep_state(kpSpec, kzSpec, InputParameters.useLog);
x = InputParameters.fstart;

% Accepted roots in the order the sweep visited them.  This mirrors what the
% engine's own XO/PO/ZO predictor does -- it too survives line boundaries --
% and is what makes a single-kperp kpar sweep trackable: the engine has no
% history of its own there, because every kpar point is a new outer line.
history = zeros(1, 0);

diag = local_diag();
% Cell arrays rather than struct arrays: struct([]) is a 0x0 with no fields,
% so growing it by assigning a filled struct fails.
perPoint = {};

guard = 0;
while true
    guard = guard + 1;
    if guard > 200000
        error('whamp_track_test:Runaway', 'The sweep did not terminate.');
    end

    ip = whampkit.grid_index(sweep.PLG, kpSpec, numel(kp));
    iz = whampkit.grid_index(sweep.ZLG, kzSpec, numel(kz));

    [solution, status, pointDiag] = local_solve_point( ...
        plasma, sweep.p, sweep.z, x, opt, history);

    if status == 1
        values = whampkit.record_point(values, ip, iz, solution);
        flagFound(ip, iz) = 1;
        if pointDiag.flagged
            flagSuspect(ip, iz) = 1;
            diag.flagged = diag.flagged + 1;
        end
        if pointDiag.rescued
            diag.rescued = diag.rescued + 1;
        end
        history(end+1) = solution.X; %#ok<AGROW>
        sweep = whampkit.remember_point(sweep, solution);
    elseif status == 2
        flagDamped(ip, iz) = 1;
    else
        flagNoConvergence(ip, iz) = 1;
    end

    % Effort counters accumulate on *every* visited point, not only the ones
    % that succeeded.  The points that failed are precisely where the seed
    % list was worked hardest, so counting only successes makes the
    % diagnostics report the easy cases and hide the interesting ones.
    diag.rejected = diag.rejected + pointDiag.rejected;
    diag.clampedSteps = diag.clampedSteps + pointDiag.clamped;
    diag.attemptsPerPoint(end+1) = pointDiag.attempts; %#ok<AGROW>

    perPoint{end+1} = local_packPoint(ip, iz, sweep.p, sweep.z, x, ...
        status, pointDiag); %#ok<AGROW>

    sweep = whampkit.advance_sweep(sweep, kpSpec, kzSpec, cycleZFirst, status ~= 1);
    if sweep.done
        break
    end
    x = whampkit.next_frequency(sweep);
end

Output = local_assemble(InputParameters, PlasmaModel, ...
    kpPhysical, kzPhysical, values, flagFound, flagDamped, ...
    flagNoConvergence, flagSuspect);
diag.converged = double(nnz(flagFound));
diag.pointsVisited = numel(perPoint);
if ~isempty(diag.attemptsPerPoint)
    diag.meanAttempts = mean(diag.attemptsPerPoint);
else
    diag.meanAttempts = NaN;
end
diag.perPoint = perPoint;
Output.diagnostics = diag;
end

% ---------------------------------------------------------------- Version 1

function opt = local_options(Options, InputParameters)
opt = struct();
opt.trustFactor = local_get(Options, 'trustFactor', 3);
opt.jumpFactor = local_get(Options, 'jumpFactor', 6);
opt.relFloor = local_get(Options, 'relFloor', 0.25);
opt.scanOffsets = local_get(Options, 'scanOffsets', ...
    [0.5 -0.5 1 -1 2 -2] * 0.5);
opt.useScan = local_get(Options, 'useScan', true);
opt.maxIterations = local_get(Options, 'maxIterations', ...
    InputParameters.maxIterations);
end

function v = local_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function diag = local_diag()
diag = struct('attemptsPerPoint',zeros(1,0),'clampedSteps',0, ...
    'rejected',0,'rescued',0,'flagged',0);
end

function [solution, status, pointDiag] = local_solve_point( ...
        plasma, p, z, xseed, opt, history)
%SOLVE_POINT Try a priority list of seeds, keep the first credible root.
%
% The first seed is always the engine's own extrapolation, so on a
% well-behaved line this costs the same as WHAMP.  Fallbacks only run when
% that seed fails to converge or produces a root the continuity test rejects.
%
% Everything is measured against the *branch*, not against the seed.  That is
% the whole point: the seed is built by extrapolating the previous point, so
% when it is wrong it is wrong by an amount of the same order as the jump it
% is about to cause.  Deriving the trust radius or the scan offsets from it
% makes the test circular -- it lets through precisely the failure it is
% meant to catch.  The branch's own recent step size is independent of the
% seed, so it still says something when the seed has already gone bad.

% Note the nested cell for 'tried': struct('f',{}) builds a 0x0 struct
% array, which would leave every field empty.
pointDiag = struct('attempts',0,'clamped',0,'rejected',0, ...
    'rescued',false,'flagged',false,'seedUsed',xseed,'tried',{{}});

xPrev = [];
if ~isempty(history)
    xPrev = history(end);
end
scale = local_stepScale(history, xPrev, opt);

% The engine's own seed is tried unconstrained, so a well-behaved line costs
% exactly what WHAMP costs and lands on exactly the same root.  The trust
% disc only applies to the branch-centred retries, where it is centred on the
% branch and wide enough for a few normal steps.
radius = Inf;
if ~isempty(xPrev)
    radius = opt.trustFactor * scale;
end

seeds = local_seedList(xseed, xPrev, scale, opt);

solution = [];
status = 0;
% The status of the engine's own seed, kept so the point can be reported with
% WHAMP's own semantics.  "Too heavily damped" is a physical statement about
% the plasma and "did not converge" is a numerical one about the iteration;
% the multi-start loop would otherwise flatten the first into the second and
% Output.flagTooHeavilyDamped would never be set.  The engine seed is the
% right one to ask, because it is exactly the solve WHAMP would have done had
% it been on this branch.
engineStatus = 0;
for s = 1:numel(seeds)
    if ~isfinite(seeds(s)) || seeds(s) <= 0
        continue
    end
    pointDiag.attempts = pointDiag.attempts + 1;
    if s == 1
        [sol, st, info] = whampkit.solvePoint(plasma, p, z, ...
            seeds(s), opt.maxIterations, Inf);
    else
        [sol, st, info] = whampkit.solvePoint(plasma, p, z, ...
            seeds(s), opt.maxIterations, radius, xPrev);
    end
    pointDiag.clamped = pointDiag.clamped + info.clamped;
    pointDiag.tried{end+1} = struct('seed',seeds(s),'status',st, ...
        'root',local_safeRoot(sol),'displacement',info.displacement); %#ok<AGROW>

    if s == 1
        engineStatus = st;
    end
    if st ~= 1
        continue
    end
    if local_continuous(sol.X, xPrev, scale, opt)
        solution = sol;
        status = 1;
        pointDiag.seedUsed = seeds(s);
        pointDiag.rescued = s > 1;
        return
    end
    % Converged but discontinuous: remember it as a last resort and keep
    % looking for something better.
    pointDiag.rejected = pointDiag.rejected + 1;
    if isempty(solution)
        solution = sol;
        status = 1;
        pointDiag.seedUsed = seeds(s);
        pointDiag.flagged = true;
    end
end

% Nothing credible was found.  There is no point in a further "just do what
% WHAMP does" fallback: the first seed above already *is* WHAMP's solve,
% unconstrained from the same extrapolated seed, so any extra attempt would be
% a byte-for-byte repeat.  When this routine fails at a point where WHAMP
% succeeds, the two have already diverged onto different branches -- WHAMP
% jumped and this solver did not -- and the coverage counts in the comparison
% are therefore not directly comparable from that point on.
%
% Preserve WHAMP's distinction between "too heavily damped" and "did not
% converge": if the engine's own seed was refused on damping grounds, the
% point is physically damped, not merely numerically stubborn.
if engineStatus == 2
    status = 2;
end
end

function scale = local_stepScale(history, xPrev, opt)
%STEPSCALE  How far this branch normally moves between neighbouring points.
%
%   The median of the recent steps is robust to a jump that has already been
%   recorded, which a plain mean or a last-step estimate is not.  The
%   relative floor matters: on a nearly flat line the median tends to zero
%   and would then reject ordinary round-off.
if isempty(xPrev)
    scale = Inf;
    return
end
if numel(history) >= 3
    tail = history(max(1, end - 7):end);
    scale = median(abs(diff(tail)));
else
    scale = 0;
end
scale = max(scale, opt.relFloor * abs(xPrev));
end

function seeds = local_seedList(xseed, xPrev, scale, opt)
%SEEDLIST Candidate starting frequencies, best first.
%
%   Ordered so that the cheap, familiar case runs first.  The scan offsets
%   are applied to the *branch* value and scaled by the branch's own step,
%   so they sample the neighbourhood the root is expected to be found in
%   rather than the neighbourhood of a seed that may already be lost.
seeds = xseed;
if ~isempty(xPrev)
    seeds(end+1) = xPrev; %#ok<AGROW>
    if isfinite(scale)
        for k = 1:numel(opt.scanOffsets)
            seeds(end+1) = xPrev + opt.scanOffsets(k) * scale; %#ok<AGROW>
        end
    end
end
seeds = seeds(:).';
end

function ok = local_continuous(x, xPrev, scale, opt)
%CONTINUITY  Is this root within a plausible step of the branch?
ok = true;
if isempty(xPrev)
    return
end
ok = abs(x - xPrev) <= opt.jumpFactor * scale;
end

function r = local_safeRoot(sol)
if isstruct(sol) && isfield(sol, 'X')
    r = sol.X;
else
    r = NaN;
end
end

function Output = local_assemble(InputParameters, PlasmaModel, ...
        kpPhysical, kzPhysical, values, flagFound, flagDamped, ...
        flagNoConvergence, flagSuspect)
resultNames = {'f','Ex','Ey','Ez','Bx','By','Bz','Sx','Sy','Sz', ...
    'EB','VGP','VGZ','SGP','SGZ','u'};
Output.InputParameters = InputParameters;
Output.PlasmaModel = PlasmaModel;
Output.kperp = kpPhysical;
Output.kpar = kzPhysical;
for k = 1:numel(resultNames)
    name = resultNames{k};
    Output.(name) = values.(name).';
end
Output.flagSolutionFound = flagFound.';
Output.flagTooHeavilyDamped = flagDamped.';
Output.flagNoConvergence = flagNoConvergence.';
Output.flagSuspectBranch = flagSuspect.';
missing = Output.flagNoConvergence == 1;
for k = 1:numel(resultNames)
    name = resultNames{k};
    Output.(name)(missing) = NaN;
end
end

function p = local_packPoint(ip, iz, p, z, seed, status, pointDiag)
p = struct('ip',ip,'iz',iz,'p',p,'z',z,'seed',seed,'status',status, ...
    'attempts',pointDiag.attempts,'clamped',pointDiag.clamped, ...
    'rejected',pointDiag.rejected,'rescued',pointDiag.rescued, ...
    'flagged',pointDiag.flagged,'seedUsed',pointDiag.seedUsed, ...
    'tried',{pointDiag.tried});
end
