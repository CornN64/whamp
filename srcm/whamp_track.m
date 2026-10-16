function Output = whamp_track(PlasmaModel, InputParameters, Options)
%WHAMP_TRACK Branch-tracking WHAMP solver, as a single file.
%   Self-contained: the whampkit routines this solver needs are inlined below
%   as kit_* local functions, so the +whampkit package need not be on the
%   path.  Behaviour is identical to the modular whamp_track_test.m.
%   OUTPUT = WHAMP_TRACK(PlasmaModel, InputParameters) walks the same
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
%   OUTPUT = WHAMP_TRACK(...) with no arguments uses the default plasma.
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
%   The physics is the whampkit code inlined below, so every
%   version evaluates exactly the same dispersion function.  WHAMP itself is
%   left untouched and its bit-for-bit regression baseline still applies.

if nargin < 1 || isempty(PlasmaModel)
    PlasmaModel = kit_complete_plasma(kit_default_plasma());
else
    PlasmaModel = kit_complete_plasma(PlasmaModel);
end
if nargin < 2 || isempty(InputParameters)
    InputParameters = kit_complete_parameters(kit_default_parameters());
else
    InputParameters = kit_complete_parameters(InputParameters);
end
if nargin < 3 || isempty(Options)
    Options = struct();
end

opt = local_options(Options, InputParameters);

[mass, density, temperature, anisotropy, delta, drift] = ...
    kit_species_arrays(PlasmaModel.Species);
plasma = kit_prepare_plasma(mass, density, temperature, anisotropy, ...
    delta, drift, 0.0279928 * PlasmaModel.B);

kpSpec = kit_axis_spec(InputParameters.kperp, 'kperp');
kzSpec = kit_axis_spec(InputParameters.kpar, 'kpar');
kp = kit_fill_axis(kpSpec);
kz = kit_fill_axis(kzSpec);
kpPhysical = kit_physical_axis_value(kp, InputParameters.useLog);
kzPhysical = kit_physical_axis_value(kz, InputParameters.useLog);

shape = [numel(kp), numel(kz)];
values = kit_emptyValues(shape);
flagFound = zeros(shape, 'int32');
flagDamped = zeros(shape, 'int32');
flagNoConvergence = zeros(shape, 'int32');
flagSuspect = zeros(shape, 'int32');

cycleZFirst = InputParameters.varyKzFirst ~= 0;
sweep = kit_sweep_state(kpSpec, kzSpec, InputParameters.useLog);
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
        error('whamp_track:Runaway', 'The sweep did not terminate.');
    end

    ip = kit_grid_index(sweep.PLG, kpSpec, numel(kp));
    iz = kit_grid_index(sweep.ZLG, kzSpec, numel(kz));

    [solution, status, pointDiag] = local_solve_point( ...
        plasma, sweep.p, sweep.z, x, opt, history);

    if status == 1
        values = kit_record_point(values, ip, iz, solution);
        flagFound(ip, iz) = 1;
        if pointDiag.flagged
            flagSuspect(ip, iz) = 1;
            diag.flagged = diag.flagged + 1;
        end
        if pointDiag.rescued
            diag.rescued = diag.rescued + 1;
        end
        history(end+1) = solution.X; %#ok<AGROW>
        sweep = kit_remember_point(sweep, solution);
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

    sweep = kit_advance_sweep(sweep, kpSpec, kzSpec, cycleZFirst, status ~= 1);
    if sweep.done
        break
    end
    x = kit_next_frequency(sweep);
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
        [sol, st, info] = kit_solvePoint(plasma, p, z, ...
            seeds(s), opt.maxIterations, Inf);
    else
        [sol, st, info] = kit_solvePoint(plasma, p, z, ...
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

% ====== Inlined from the +whampkit package (behaviour unchanged) ======

%% ---------------- from +whampkit/advance_sweep.m ----------------
function sweep = kit_advance_sweep(sweep, kpSpec, kzSpec, cycleZFirst, stopCycling)
%ADVANCE_SWEEP Move to the next point of the two-axis sweep.
%   The outer axis is walked by repeated addition of its step in the sweep
%   coordinate, exactly as the engine does.  A point that did not converge,
%   or that was too heavily damped, abandons the rest of the *inner* axis;
%   the remaining outer points are still visited and the sweep ends once the
%   outer axis is exhausted.  Reaching a new outer point resets KV so that
%   the line is re-anchored.
sweepEnd = 1e99;
if stopCycling
    if cycleZFirst
        sweep.ZLG = sweepEnd;
    else
        sweep.PLG = sweepEnd;
    end
end
if ~cycleZFirst
    sweep.PLG = sweep.PLG + kpSpec(3);
    if sweep.PLG >= kpSpec(1) && sweep.PLG <= kpSpec(2)
        sweep.p = sweep.PLG + sweep.PZL * (10^sweep.PLG - sweep.PLG);
    else
        sweep.ZLG = sweep.ZLG + kzSpec(3);
        if sweep.ZLG < kzSpec(1) || sweep.ZLG > kzSpec(2)
            sweep.done = true;
            return
        end
        sweep.KV = 1;
        sweep.PLG = sweep.PLO;
        sweep.p = sweep.PVO;
        sweep.z = sweep.ZLG + sweep.PZL * (10^sweep.ZLG - sweep.ZLG);
    end
else
    sweep.ZLG = sweep.ZLG + kzSpec(3);
    if sweep.ZLG >= kzSpec(1) && sweep.ZLG <= kzSpec(2)
        sweep.z = sweep.ZLG + sweep.PZL * (10^sweep.ZLG - sweep.ZLG);
    else
        sweep.PLG = sweep.PLG + kpSpec(3);
        if sweep.PLG < kpSpec(1) || sweep.PLG > kpSpec(2)
            sweep.done = true;
            return
        end
        sweep.KV = 1;
        sweep.ZLG = sweep.ZLO;
        sweep.z = sweep.ZVO;
        sweep.p = sweep.PLG + sweep.PZL * (10^sweep.PLG - sweep.PLG);
    end
end
end

%% ---------------- from +whampkit/axis_spec.m ----------------
function spec = kit_axis_spec(axis, name)
%AXIS_SPEC Reduce a kperp/kpar request to the engine's [start end step].
%   The caller writes [start step end]; the engine keeps the first and third
%   numbers as the bounds and the second as the step, so the triple returned
%   here is [PM1 PM2 PM3].  The sweep walks PM3 at a time from PM1, or from
%   PM2 when PM3 is negative.  A bare scalar means one point, and the engine
%   default step of 10 is kept for it.
if ~isnumeric(axis) || ~isreal(axis) || isempty(axis) || ...
        ~(isscalar(axis) || numel(axis) == 3) || any(~isfinite(axis(:)))
    error('whamp:InvalidWaveNumber', ...
        'InputParameters.%s must be a finite scalar or [start step end].', name);
end
if isscalar(axis)
    spec = [axis axis 10];
else
    spec = [axis(1), axis(3), axis(2)];
end
if spec(1) ~= spec(2) && spec(3) == 0
    error('whamp:InvalidWaveNumberStep', ...
        'The %s step must be nonzero unless start equals end.', name);
end
end

%% ---------------- from +whampkit/complete_parameters.m ----------------
function parameters = kit_complete_parameters(parameters)
%COMPLETE_PARAMETERS Fill in omitted parameter fields and validate them.
defaults = kit_default_parameters();
if ~isstruct(parameters) || ~isscalar(parameters)
    error('whamp:InvalidInputParameters', ...
        'InputParameters must be a scalar structure.');
end
fieldNames = fieldnames(defaults);
for k = 1:numel(fieldNames)
    if ~isfield(parameters, fieldNames{k})
        parameters.(fieldNames{k}) = defaults.(fieldNames{k});
    end
end

for name = {'fstart','varyKzFirst','useLog'}
    if ~kit_is_finite_real_scalar(parameters.(name{1}))
        error('whamp:InvalidInputParameters', ...
            'InputParameters.%s must be a finite real scalar.', name{1});
    end
end
if parameters.fstart <= 0
    error('whamp:InvalidStartFrequency', 'fstart must be positive.');
end
if ~ismember(parameters.useLog, [0 1])
    error('whamp:InvalidLogOption', 'useLog must be 0 or 1.');
end
end

%% ---------------- from +whampkit/complete_plasma.m ----------------
function model = kit_complete_plasma(model)
%COMPLETE_PLASMA Fill in omitted plasma fields and validate what remains.
%   The magnetic field is checked last so that a malformed species reports
%   its own error first.
defaults = kit_default_plasma();
if ~isstruct(model) || ~isscalar(model)
    error('whamp:InvalidPlasmaModel', 'PlasmaModel must be a scalar structure.');
end
if ~isfield(model, 'B')
    model.B = defaults.B;
end
model.Species = kit_normalise_species(model, defaults);

speciesDefaults = struct('m',0,'n',1,'t',1,'a',1,'b',0,'d',1,'vd',0);
fieldNames = fieldnames(speciesDefaults);
for j = 1:numel(model.Species)
    if ~isstruct(model.Species{j}) || ~isscalar(model.Species{j})
        error('whamp:InvalidSpecies', 'Each species must be a scalar structure.');
    end
    for k = 1:numel(fieldNames)
        name = fieldNames{k};
        if ~isfield(model.Species{j}, name)
            model.Species{j}.(name) = speciesDefaults.(name);
        end
        if ~kit_is_finite_real_scalar(model.Species{j}.(name))
            error('whamp:InvalidSpecies', ...
                'Species field %s must be a finite real scalar.', name);
        end
    end
end

if ~isnumeric(model.B) || ~isscalar(model.B) || ~isfinite(model.B) || model.B <= 0
    error('whamp:InvalidMagneticField', 'PlasmaModel.B must be a positive finite scalar.');
end
end

%% ---------------- from +whampkit/dcTerm.m ----------------
function dc = kit_dcTerm(E, col)
%DCTERM The DC term of the derivative of the dispersion determinant.
%   Column 1 of E holds the dielectric tensor and each later column its
%   derivative with respect to one variable (2: frequency, 3: k-parallel,
%   4: k-perpendicular).  The three call sites differ only in that column, so
%   the expression is written once here.
dc = (E(1,col)*E(6,1)+E(1,1)*E(6,col)-2*E(3,1)*E(3,col))*E(4,1) ...
    + (E(1,1)*E(6,1)-E(3,1)^2)*E(4,col) ...
    + E(5,1)*(E(1,col)*E(5,1)+2*E(1,1)*E(5,col) ...
    + E(2,1)*E(3,col)+2*E(2,col)*E(3,1)) ...
    + E(2,1)*(E(6,col)*E(2,1)+2*E(6,1)*E(2,col) ...
    + E(5,1)*E(3,col)+2*E(5,col)*E(3,1));
end

% Susceptibility of one species, summed over the eight poles of the
% Gaussian (B- and C-) representation of the plasma dispersion function.
%
% Returns result.XSI(6,order) accumulated over the poles and result.IERR.
% The six rows are the tensor components 11, 12, 13, 22, 23, 33 and the
% columns are successive derivatives (1 = value, 2 = d/dx, 3 = d/dz,
% 4 = d/dp^2).
%
% In CHI() the relativistic factor is read as AA(J,IB), so the caller's IB
% loop picks which AA column applies; here that choice has already been made
% and passed in as `a`.  The remaining C arguments B and DELTA are likewise
% folded into the caller's weighting and are not needed again.

%% ---------------- from +whampkit/default_parameters.m ----------------
function parameters = kit_default_parameters()
%DEFAULT_PARAMETERS One point at kperp = 0, kpar = 0.0022, from f = 0.1.
parameters = struct('fstart',0.1,'kperp',[0 10 0], ...
    'kpar',[0.0022 10 0.0022],'varyKzFirst',1,'useLog',0, ...
    'maxIterations',50);
end

%% ---------------- from +whampkit/default_plasma.m ----------------
function model = kit_default_plasma()
%DEFAULT_PLASMA The built-in model: 100 nT with an O+ beam on electrons.
%   Density is in cm^-3, temperature in eV, mass in proton masses (zero
%   denotes electrons), and vd the field-aligned drift in thermal speeds.
model = struct('B', 100);
model.Species = { ...
    struct('m',16,'n',1,'t',10,'a',5,'b',0,'d',1,'vd',1), ...
    struct('m',0,'n',1,'t',100,'a',1,'b',0,'d',1,'vd',0)};
end

%% ---------------- from +whampkit/dispersion.m ----------------
function d = kit_dispersion(state, order)
%DISPERSION Dielectric tensor and dispersion determinant at one frequency.
%   Transcribes DIFU() in difu.c.  Column 1 of E is the 6-element symmetric
%   dielectric tensor (11,12,13,22,23,33) built by summing the species
%   susceptibilities; columns 2..ORDER hold its derivatives with respect to
%   the frequency, k-parallel and k-perpendicular.  ORDER selects how much of
%   that work is done, since the Newton iteration only needs D and DX.
%
%   Returns D (the determinant of the wave-normal form), the derivatives
%   DX/DZ/DP, and -- at ORDER 4 -- the polarisation and energy quantities.
%   IERR is set when a species is too heavily damped for the expansion.
plasma = state.plasma;
E = complex(zeros(6,4));
E([1 4 6],1) = 1;
IERR = false;
for j = 1:plasma.n
    xp = state.XP(j);
    E([1 4 6],1) = E([1 4 6],1) - xp;
    a = plasma.aa(j,1);
    b = plasma.aa(j,2);
    delta = plasma.delta(j);
    if a == b
        b = 0;
        delta = 1;
    end
    factor = xp / (a * (a - b));
    q = a - delta * b;
    for ib = 1:2
        xsi = kit_susceptibility(state, j, order, a);
        if xsi.IERR
            IERR = true;
            break
        end
        E(:,1:order) = E(:,1:order) + factor * q * xsi.XSI(:,1:order);
        if ib == 2
            break
        end
        % Second (loss-cone) term; absent for an isotropic species.
        q = (delta - 1) * a;
        if q == 0
            break
        end
    end
    if IERR
        break
    end
end
if IERR
    d = struct('IERR',true);
    return
end

% Refractive-index components of the first species, in units of the thermal
% speed, and the three quadratic forms Av/Bv/Cv of the determinant.
u1 = state.PP(1) * plasma.cv / state.XX(1);
u3 = state.ZZ(1) * plasma.cv / state.XX(1);
u12 = u1^2;
u32 = u3^2;
u2 = u12 + u32;
u13 = 2 * u1 * u3;
Av = u12 * E(1,1) + u13 * E(3,1) + u32 * E(6,1);
Bv = u2 * (E(1,1) * E(6,1) - E(3,1)^2) ...
    + (u3 * E(5,1) - u1 * E(2,1))^2;
Cv = (E(1,1) * E(6,1) - E(3,1)^2) * E(4,1) ...
    + E(6,1) * E(2,1)^2 ...
    + (E(1,1) * E(5,1) + 2 * E(2,1) * E(3,1)) * E(5,1);
D = (u2 - E(4,1)) * Av - Bv + Cv;
if order == 1
    d = struct('D',D,'IERR',false);
    return
end

% Frequency derivative.  The -2*[...] term converts the frequency column of
% E into the derivative of the cold-plasma part that the determinant uses.
E(:,2) = E(:,2) - 2 * [E(1,1)-1; E(2,1); E(3,1); ...
    E(4,1)-1; E(5,1); E(6,1)-1];
da = (E(1,2)-2*E(1,1))*u12 + u13*(E(3,2)-2*E(3,1)) ...
    + (E(6,2)-2*E(6,1))*u32;
db = 2*(u3*E(5,1)-u1*E(2,1)) * ...
    (u3*(E(5,2)-E(5,1))-u1*(E(2,2)-E(2,1))) ...
    + u2*((E(1,2)-E(1,1))*E(6,1) - 2*(E(3,2)-E(3,1))*E(3,1) ...
    + (E(6,2)-E(6,1))*E(1,1));
dc = kit_dcTerm(E, 2);
DX = ((u2-E(4,1))*da-(2*u2+E(4,2))*Av-db+dc)/state.XX(1);
if order == 2
    d = struct('D',D,'DX',DX,'IERR',false);
    return
end
DZ = 0;
if state.ZZ(1) ~= 0
    da = u12*E(1,3)+u13*(E(3,3)+E(3,1))+u32*(E(6,3)+2*E(6,1));
    db = 2*(u3*E(5,1)-u1*E(2,1)) * ...
        (u3*(E(5,3)+E(5,1))-u1*E(2,3)) ...
        + 2*u32*(E(1,1)*E(6,1)-E(3,1)^2) ...
        + u2*(E(1,3)*E(6,1)+E(1,1)*E(6,3)-2*E(3,1)*E(3,3));
    dc = kit_dcTerm(E, 3);
    DZ = ((u2-E(4,1))*da+(2*u32-E(4,3))*Av-db+dc)/state.ZZ(1);
end
if order == 3
    d = struct('D',D,'DX',DX,'DZ',DZ,'IERR',false);
    return
end
DP = 0;
if state.PP(1) ~= 0
    da = u12*(E(1,4)+2*E(1,1))+u13*(E(3,4)+E(3,1))+u32*E(6,4);
    db = 2*(u3*E(5,1)-u1*E(2,1)) * ...
        (u3*E(5,4)-u1*(E(2,4)+E(2,1))) ...
        + 2*u12*(E(1,1)*E(6,1)-E(3,1)^2) ...
        + u2*(E(1,4)*E(6,1)+E(1,1)*E(6,4)-2*E(3,1)*E(3,4));
    dc = kit_dcTerm(E, 4);
    DP = ((u2-E(4,1))*da+(2*u12-E(4,4))*Av-db+dc)/state.PP(1);
end
[EFL,BFL,u] = kit_wave_energy(E,u1,u3,u2,u12,u32);
d = struct('D',D,'DX',DX,'DZ',DZ,'DP',DP,'EFL',EFL, ...
    'BFL',BFL,'u',u,'IERR',false);
end

%% ---------------- from +whampkit/dispersion_state.m ----------------
function state = kit_dispersion_state(plasma, p, z, x)
%DISPERSION_STATE Frequency-dependent arguments of the dispersion function.
%   XP is the species density normalised by the plasma frequency, XX the
%   frequency in units of each species gyrofrequency, and PP/ZZ the
%   wavenumber in units of the corresponding thermal speed.
omega = (x * plasma.xa)^2;
state = struct('X',x, ...
    'XP', plasma.density ./ plasma.dek ./ plasma.ren ./ omega, ...
    'XX', x .* plasma.ren, ...
    'PP', p .* plasma.st, ...
    'ZZ', z .* plasma.st, ...
    'plasma', plasma);
end

%% ---------------- from +whampkit/emptyValues.m ----------------
function values = kit_emptyValues(shape)
%EMPTYVALUES Allocate the result grids that whamp.m fills.
names = {'f','Ex','Ey','Ez','Bx','By','Bz','Sx','Sy','Sz', ...
    'EB','VGP','VGZ','SGP','SGZ','u'};
values = struct();
for k = 1:numel(names)
    values.(names{k}) = complex(zeros(shape));
end
end

%% ---------------- from +whampkit/fill_axis.m ----------------
function axis = kit_fill_axis(spec)
%FILL_AXIS Wave numbers of the stored solutions, in sweep order.
%   Matches the engine's axis_size()/kit_fill_axis() pair: the axis begins where
%   the sweep begins and advances by the step, so axis(i) is the wave number
%   of the i-th stored solution.
if spec(1) == spec(2)
    axis = spec(1);
    return
end
count = 1 + floor(abs(spec(2) - spec(1)) / abs(spec(3)));
if count > 100000
    error('whamp:InvalidWaveNumberRange', 'Wave-number axis has an invalid size.');
end
axis = kit_sweep_origin(spec) + (0:count-1) * spec(3);
end

%% ---------------- from +whampkit/fnint.m ----------------
function n = kit_fnint(x)
%FNINT Fortran NINT: nearest integer, ties away from zero.
if x >= 0
    n = floor(x + 0.5);
else
    n = ceil(x - 0.5);
end
end

%% ---------------- from +whampkit/gauss_legendre_16.m ----------------
function [abscissa, weight] = kit_gauss_legendre_16()
%GAUSS_LEGENDRE_16 Abscissas and weights for the 16-point Gauss-Legendre rule.
%   The two halves of the interval are interleaved in the order RINT() uses.
abscissa = [-.989400934991649,-.944575023073232,-.865631202387831, ...
    -.755404408355003,-.617876244402643,-.458016777657227, ...
    -.281603550779258,-.095012509837637,.989400934991649, ...
    .944575023073232,.865631202387831,.755404408355003, ...
    .617876244402643,.458016777657227,.281603550779258,.095012509837637];
weight = [.027152459411754,.062253523938647,.095158511682492, ...
    .124628971255533,.149595988816576,.169156519395002, ...
    .182603415044923,.189450610455068,.027152459411754, ...
    .062253523938647,.095158511682492,.124628971255533, ...
    .149595988816576,.169156519395002,.182603415044923,.189450610455068];
end

%% ---------------- from +whampkit/grid_index.m ----------------
function index = kit_grid_index(value, spec, count)
% Nearest grid index of a walked wave number, measured from where the sweep
% started (the engine measures from PM1, which runs off the end of its array
% for a negative step).  Clamped so a rounding slip cannot index out of
% range.
if spec(3) == 0
    index = 1;
    return
end
index = 1 + kit_fnint((value - kit_sweep_origin(spec)) / spec(3));
index = min(max(index, 1), count);
end

%% ---------------- from +whampkit/is_finite_real_scalar.m ----------------
function tf = kit_is_finite_real_scalar(value)
%IS_FINITE_REAL_SCALAR True for a real, finite numeric scalar.
%   Written with short-circuit && so that isfinite is never handed a cell or
%   character array, which would raise instead of reporting false.
tf = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end

%% ---------------- from +whampkit/newtonRoot.m ----------------
function [state, rootFound, status, info] = kit_newtonRoot(plasma, state, maxIterations, radius, center)
%NEWTONROOT Newton search for a dispersion root, with an optional trust region.
%   [STATE, ROOTFOUND, STATUS] = kit_newtonRoot(PLASMA, STATE, MAXIT)
%   reproduces whamp.m's NEWTON_ROOT iteration exactly (radius = Inf).
%
%   [..., RADIUS] restricts the search to a disc of radius RADIUS centred on
%   the *starting* frequency STATE.X.  The cap is applied to the total
%   displacement from that seed rather than to each individual correction:
%   a per-step cap would let the repeated outer restarts creep arbitrarily
%   far from the seed, which is precisely the behaviour that lets the
%   iteration escape into a neighbouring root's basin.
%
%   STATUS is 2 when the susceptibility flags the point as too heavily
%   damped, otherwise it is unchanged from the reference routine.
%
%   INFO is a struct with the diagnostics used by the continuity checks:
%     .steps        Newton iterations used
%     .clamped      number of corrections limited by the trust region
%     .displacement |X_root - X_seed|
%     .residual     |D| at the returned point
%
%   Why a trust region helps: D(w) is meromorphic, with a pole for every
%   species gyroharmonic (bl/y with y = x - dl*z, and cot(pi*y) inside
%   r_asymptotic).  Newton started near such a pole takes an enormous step
%   and lands in whichever basin it happens to fall in.  Capping the walk
%   keeps the solve on the branch the seed came from.

if nargin < 4 || isempty(radius)
    radius = Inf;
end
% The trust disc is normally centred on the seed.  A caller that already
% knows where the branch is can centre it somewhere else -- which matters
% because the whole point of a trust region here is to hold the iteration
% near the *branch*, and the seed it is given may itself be off the branch.
if nargin < 5 || isempty(center)
    center = state.X;
end

freqRelTol = 1e-6;
freqAbsTol = 1e-6;
minimaTol = 1e-5;

seed = state.X;
rootFound = false;
status = 0;
endIteration = false;
cycleOuter = false;
clamped = 0;
steps = 0;

for iteration = 1:maxIterations
    oldMagnitude = abs(state.d.D);
    irk = 1;
    correction = state.d.D / state.d.DX;
    while true
        steps = steps + 1;
        [correction, wasClamped, exhausted] = kit_limitStep(correction, state.X, center, radius);
        clamped = clamped + wasClamped;
        % The trust budget is spent.  Reporting the resulting zero correction
        % as convergence would be wrong: it is small because it was truncated,
        % not because D is.  Without this test the iteration "converges" while
        % pinned against the wall of the disc, and the caller records a point
        % on the wall -- which is how a fallback seed two steps from the
        % branch could be accepted as a root.  Give up instead and let the
        % caller try another seed.
        if exhausted
            rootFound = false;
            endIteration = true;
            break
        end

        state = kit_evaluateAt(plasma, state, state.X - correction);
        if state.d.IERR
            status = 2;
            info = kit_makeInfo(steps, clamped, state, seed);
            return
        end
        if abs(state.d.D) < oldMagnitude
            % A clamped correction is small because it was truncated, so its
            % size says nothing about how close we are to a zero.  Keep
            % iterating (the trust region still holds us on the branch) and
            % let the residual decide.
            convergedStep = ~wasClamped && (abs(correction) <= freqRelTol * abs(state.X) || ...
                abs(correction) < freqAbsTol);
            if convergedStep
                rootFound = true;
                % At least two steps must have been made before accepting.
                if iteration >= 2
                    endIteration = true;
                else
                    cycleOuter = true;
                end
            else
                cycleOuter = true;
            end
            break
        end
        % Uphill: step back to where we came from and halve the correction.
        % Only X is restored -- D still describes the trial point, which is
        % what the local-minima test below compares against.
        state.X = state.X + correction;
        correction = correction / 2;
        if irk > 3 && (abs(state.d.D) - oldMagnitude) / oldMagnitude < minimaTol
            rootFound = false;
            endIteration = true;
            break
        end
        irk = irk + 1;
        if irk > maxIterations
            endIteration = true;
            break
        end
    end
    if endIteration
        break
    end
    if ~cycleOuter
        break
    end
end

info = kit_makeInfo(steps, clamped, state, seed);
end

function [correction, wasClamped, exhausted] = kit_limitStep(correction, x, seed, radius)
%LIMITSTEP Keep the trial frequency inside the trust disc around SEED.
%   The budget is what is left of RADIUS after the distance already
%   travelled, so the cap bounds the *total* excursion from the seed.
%   EXHAUSTED reports that no budget remains, which tells the caller the
%   iteration is pinned against the wall rather than converging.
wasClamped = 0;
exhausted = false;
if ~isfinite(radius)
    return
end
budget = radius - abs(x - seed);
if budget <= 0
    correction = 0;
    wasClamped = 1;
    exhausted = true;
    return
end
magnitude = abs(correction);
if magnitude > budget
    correction = correction / magnitude * budget;
    wasClamped = 1;
end
end

function state = kit_evaluateAt(plasma, state, x)
%EVALUATEAT Move to a new frequency and recompute the dispersion function.
state.X = x;
omega = (x * plasma.xa)^2;
state.XP = plasma.density ./ plasma.dek ./ plasma.ren ./ omega;
state.XX = x .* plasma.ren;
state.d = kit_dispersion(state, 2);
end

function info = kit_makeInfo(steps, clamped, state, seed)
%MAKEINFO  Pack the Newton diagnostics.
%   On the heavily-damped exit state.d carries only IERR, so the residual is
%   unknown rather than zero.
if isfield(state.d, 'D')
    residual = abs(state.d.D);
else
    residual = Inf;
end
info = struct('steps',steps,'clamped',clamped, ...
    'displacement',abs(state.X - seed),'residual',residual);
end

%% ---------------- from +whampkit/next_frequency.m ----------------
function x = kit_next_frequency(sweep)
%NEXT_FREQUENCY Linear extrapolation of the frequency to the next point.
if sweep.KV ~= 0
    x = sweep.XVO - ((sweep.p - sweep.PVO) * sweep.DOP ...
        + (sweep.z - sweep.ZVO) * sweep.DOZ) / sweep.DOX;
else
    x = sweep.XO - ((sweep.p - sweep.PO) * sweep.DP ...
        + (sweep.z - sweep.ZO) * sweep.DZ) / sweep.DX;
end
end

%% ---------------- from +whampkit/normalise_species.m ----------------
function species = kit_normalise_species(model, defaults)
%NORMALISE_SPECIES Present PlasmaModel.Species as a cell array of structs.
if ~isfield(model, 'Species') || isempty(model.Species)
    species = defaults.Species;
elseif isstruct(model.Species)
    species = num2cell(model.Species);
elseif iscell(model.Species)
    species = model.Species;
else
    error('whamp:InvalidSpecies', ...
        'PlasmaModel.Species must be a structure or cell array.');
end
end

%% ---------------- from +whampkit/physical_axis_value.m ----------------
function value = kit_physical_axis_value(axisValue, PZL)
%PHYSICAL_AXIS_VALUE Convert sweep coordinates to wave numbers.
if PZL == 1
    value = 10 .^ axisValue;
else
    value = axisValue;
end
end

%% ---------------- from +whampkit/pole_constants.m ----------------
function [B, C] = kit_pole_constants()
%POLE_CONSTANTS Residues and pole positions of the eight-pole approximant.
%   The Gaussian representation of the R-function used by CHI() is built from
%   these two fixed tables; B holds the residues and C the pole positions.
B = [-1.734012457471826e-2-4.630639291680322e-2i, ...
    -1.734012457471826e-2+4.630639291680322e-2i, ...
    -7.399169923225014e-1+8.395179978099844e-1i, ...
    -7.399169923225014e-1-8.395179978099844e-1i, ...
    5.840628642184073+9.536009057643667e-1i, ...
    5.840628642184073-9.536009057643667e-1i, ...
    -5.583371525286853-1.120854319126599e1i, ...
    -5.583371525286853+1.120854319126599e1i];
C = [ 2.237687789201900-1.625940856173727i, ...
    -2.237687789201900-1.625940856173727i, ...
    1.465234126106004-1.789620129162444i, ...
    -1.465234126106004-1.789620129162444i, ...
    0.8392539817232638-1.891995045765206i, ...
    -0.8392539817232638-1.891995045765206i, ...
    0.2739362226285564-1.941786875844713i, ...
    -0.2739362226285564-1.941786875844713i];
end

%% ---------------- from +whampkit/prepare_plasma.m ----------------
function plasma = kit_prepare_plasma(mass, density, temperature, aa, delta, drift, gyro)
%PREPARE_PLASMA Normalise the model into the engine's working quantities.
%   Transcribes plasma_setup() in whamp_engine.c.  Everything is expressed
%   relative to the first species: REN is the mass ratio to it, T the
%   temperature ratio scaled by REN, and ST = sqrt(T).  By the time this runs
%   the densities are in m^-3 and the temperatures in keV.
protonMassRatio = 1836.1;
debyeKm = 12405;
electronRestEnergy = 511;

% Electrons are entered with mass 0, which would divide by zero; the engine
% gives them the unit mass ratio here and lets the charge sign do the work.
ren = protonMassRatio .* mass;
ren(ren == 0) = 1;

last = find(density ~= 0, 1, 'last');
if isempty(last)
    error('whamp:EmptyPlasma', 'At least one species must have nonzero density.');
end
if temperature(1) <= 0
    error('whamp:InvalidTemperature', ...
        'The first species temperature must be positive.');
end

% A species whose two loss-cone parameters coincide has no loss cone at all.
for j = 1:last
    if aa(j,1) == aa(j,2)
        aa(j,2) = 0;
        delta(j) = 1;
    end
end

rn = ren(1);
t = temperature ./ temperature(1) .* (ren ./ rn);
if any(t(1:last) <= 0)
    error('whamp:InvalidTemperature', 'Active species temperatures must be positive.');
end
plasma.n = last;
plasma.mass = mass(1:last);
plasma.density = density(1:last);
plasma.ren = ren(1:last) ./ rn;
plasma.st = sqrt(t(1:last));
plasma.aa = aa(1:last,:);
plasma.delta = delta(1:last);
plasma.drift = drift(1:last);
plasma.dek = debyeKm * rn;
plasma.pfq = sum(density(1:last) ./ ren(1:last)) / debyeKm;
plasma.xa = gyro / rn;
% CV is the electron thermal speed in units of c, from the first species
% temperature TR = T(1)/RN with the special-relativistic factor applied.
plasma.cv = 1 / sqrt((temperature(1) / rn) * ...
    (2 * electronRestEnergy + temperature(1) / rn) / ...
    (electronRestEnergy + temperature(1) / rn)^2);
plasma.xc = gyro;
end

%% ---------------- from +whampkit/record_point.m ----------------
function values = kit_record_point(values, ip, iz, solution)
%RECORD_POINT Write one converged solution into the output grids.
%   The Poynting flux carries the engine's 10/(4*pi/2) unit conversion.
poynting = 10 / 4 / pi / 2;
values.f(ip, iz) = solution.X;
values.Ex(ip, iz) = solution.EFL(1);
values.Ey(ip, iz) = solution.EFL(2);
values.Ez(ip, iz) = solution.EFL(3);
values.Bx(ip, iz) = solution.BFL(1);
values.By(ip, iz) = solution.BFL(2);
values.Bz(ip, iz) = solution.BFL(3);
values.Sx(ip, iz) = real(solution.EFL(2) * conj(solution.BFL(3)) ...
    - solution.EFL(3) * conj(solution.BFL(2))) * poynting;
values.Sy(ip, iz) = real(solution.EFL(3) * conj(solution.BFL(1)) ...
    - solution.EFL(1) * conj(solution.BFL(3))) * poynting;
values.Sz(ip, iz) = real(solution.EFL(1) * conj(solution.BFL(2)) ...
    - solution.EFL(2) * conj(solution.BFL(1))) * poynting;
values.EB(ip, iz) = sqrt(sum(abs(solution.EFL).^2) ...
    / sum(abs(solution.BFL).^2));
values.VGP(ip, iz) = solution.VG(1);
values.VGZ(ip, iz) = solution.VG(2);
values.SGP(ip, iz) = solution.SG(1);
values.SGZ(ip, iz) = solution.SG(2);
values.u(ip, iz) = solution.u;
end

%% ---------------- from +whampkit/remember_point.m ----------------
function sweep = kit_remember_point(sweep, solution)
%REMEMBER_POINT Store a converged point as the previous-point predictor.
sweep.PO = sweep.p;
sweep.ZO = sweep.z;
sweep.XO = solution.X;
sweep.DP = solution.DP;
sweep.DZ = solution.DZ;
sweep.DX = solution.DX;
if sweep.KV ~= 0
    % First converged point of this line: every later point on it is
    % predicted from this anchor.
    sweep.XVO = solution.X;
    sweep.PVO = sweep.p;
    sweep.ZVO = sweep.z;
    sweep.PLO = sweep.PLG;
    sweep.ZLO = sweep.ZLG;
    sweep.DOX = sweep.DX;
    sweep.DOP = sweep.DP;
    sweep.DOZ = sweep.DZ;
    sweep.KV = 0;
end
end

%% ---------------- from +whampkit/r_asymptotic.m ----------------
function rc = kit_r_asymptotic(y, al)
%R_ASYMPTOTIC Asymptotic expansion of R_FUNCTION for large AL, moderate |Y|.
%   Sums the P and Q series until the terms start growing again.
y2 = y^2;
coterm = cot(pi*y);
pn = -y/al;
pyn = pn;
a = 1/(al*sqrt(2*pi*al));
qn = pi*y2*coterm*a;
qyn = qn*(2-y*pi*coterm)-y*pi^2*y2*a;
p = pn+qn;
py = pyn+qyn;
pp = -pn-1.5*qn;
ppy = -pyn-1.5*qyn;
ay = abs(y)+2;
t = 0;
c = 7.2e35;
for n = 1:100
    m = n-1;
    pyn = (pyn*(m^2-y2)-2*y2*pn)/((2*m+1)*al);
    pn = pn*(m^2-y2)/((2*m+1)*al);
    qyn = (qyn*((m+0.5)^2-y2)-2*y2*qn)/(2*n*al);
    qn = qn*((m+0.5)^2-y2)/(2*n*al);
    if m >= ay
        c = n*(abs(pn)+abs(qn));
        if c <= 1e-7*abs(pp) || c >= t
            break
        end
    end
    p = p+pn+qn;
    py = py+pyn+qyn;
    pp = pp-(n+1)*pn-(n+1.5)*qn;
    ppy = ppy-(n+1)*pyn-(n+1.5)*qyn;
    t = c;
end
rc = [p+pn+qn, py+pyn+qyn, pp+p, ppy+py];
end

%% ---------------- from +whampkit/r_function.m ----------------
function rc = kit_r_function(y, al)
%R_FUNCTION Relativistic plasma dispersion function and its derivatives.
%   RC = [R, dR/dy, R', d/dy R'] evaluated at complex Y for the real
%   perpendicular argument AL, following the method selection of RYLA():
%   a Taylor expansion for small AL or large |Y|, an asymptotic expansion
%   for large AL with moderate |Y|, and numerical integration otherwise.
ay = abs(y);
if al < 4 || ay^2 > 75*al || ay > 40+al/3
    rc = kit_r_taylor(y, al);
elseif 3*(al-10) > ay && ay^2 < 15*al
    rc = kit_r_asymptotic(y, al);
else
    rc = kit_r_integral(y, al);
end
end

%% ---------------- from +whampkit/r_integral.m ----------------
function rc = kit_r_integral(y, al)
%R_INTEGRAL Numerical integration of R_FUNCTION.
%   Two 16-point Gauss-Legendre quadratures -- over the finite interval UL
%   and over the exponential tail XO -- are combined and then corrected
%   analytically.  R is odd in Y up to the sign SIG, so a negative argument
%   is reflected and the sign restored at the end.
[abscissa, weight] = kit_gauss_legendre_16();
if real(y) < 0
    y = -y;
    sig = -1;
else
    sig = 1;
end
ya = imag(y);
yr = real(y);
ul = pi-real(2.8*y/(36+y));
coterm = cot(pi*y);
dv = pi*(1+coterm^2);          % d(cot(pi*y))/dy contribution
c = yr/al;
xo = log(c+sqrt(1+c^2));       % asinh(yr/al)
rc = complex(zeros(1,4));
for i = 1:16
    x = ul/2*(1+abscissa(i));
    z = sin(x);
    c = cos(x);
    g = yr/al*x/z;
    t = sqrt(1+g^2);
    b = log(g+t);
    g = (1/x-c/z)*g/t;
    t = al*(t*c-1);
    z = real(exp(x*ya));
    c = 0.5*(z+1/z);
    s = 0.5i*(z-1/z);
    f = coterm+g;
    h = 1-g*coterm;
    exf = exp(t-y*b);
    o = b*c+x*s;
    pv = x*c-b*s;
    xy = x*yr;
    r = (f*c+h*s)*exf;
    ry = (f*o-h*pv+dv*(c-g*s))*exf;
    rp = ((f*t-h*xy)*c+(h*t+f*xy)*s)*exf;
    rpy = (f*(t*o-xy*pv)-h*(t*pv+xy*o) ...
        +((t+xy*g)*c-(g*t-xy)*s)*dv)*exf;
    x = xo/2*(1+abscissa(i));
    z = exp(x);
    c = (z+1/z)/2-1;
    pv = exp(al*c-y*x);
    rc(1) = rc(1)+weight(i)*(ul*r+xo*pv);
    rc(2) = rc(2)-weight(i)*(ul*ry+xo*x*pv);
    rc(3) = rc(3)+weight(i)*(ul*rp+xo*al*c*pv);
    rc(4) = rc(4)-weight(i)*(ul*rpy+xo*al*x*c*pv);
end
% Analytic endpoint corrections, with the SIG reflection undone.
o = y/al;
pv = y^2/2;
rc(1) = o*(y*rc(1)/2-1)*sig;
rc(2) = 2*rc(1)+o*(pv*rc(2)+1)*sig;
rc(3) = y*o*rc(3)/2*sig;
rc(4) = 2*rc(3)+o*pv*rc(4)*sig;
end

%% ---------------- from +whampkit/r_taylor.m ----------------
function rc = kit_r_taylor(y, al)
%R_TAYLOR Taylor expansion of R_FUNCTION for small AL or large |Y|.
%   Summation stops once the term is eight orders below the leading one.
y2 = y^2;
pn = y/(y2-1);
pyn = -y*(y2+1)/(y2-1)^2;
rc = [pn pyn pn pyn];
for i = 2:100
    coterm = (2*i-1)/(y2-i^2)*al;
    pyn = coterm*(pyn-2*y2/(y2-i^2)*pn);
    pn = coterm*pn;
    rc = rc + [pn pyn i*pn i*pyn];
    if abs(pn)*1e8 < abs(rc(1))
        break
    end
end
end

%% ---------------- from +whampkit/solvePoint.m ----------------
function [solution, status, info] = kit_solvePoint(plasma, p, z, x, maxIterations, radius, center)
%SOLVEPOINT Solve the dispersion relation at one (p, z) pair.
%   [SOLUTION, STATUS] = kit_solvePoint(PLASMA, P, Z, X, MAXIT) is the
%   kit's counterpart of whamp.m's SOLVE_POINT and returns the same fields.
%   STATUS is 1 for a converged root, 2 for a point the susceptibility
%   reports as too heavily damped, and 0 when the iteration gave up.
%
%   [..., RADIUS] passes a Newton trust-region radius through to
%   kit_newtonRoot.  RADIUS = Inf (the default) reproduces the
%   reference iteration exactly.
%
%   [..., INFO] returns the Newton diagnostics from kit_newtonRoot
%   (.steps .clamped .displacement .residual), used by the continuity check.

if nargin < 5 || isempty(maxIterations)
    maxIterations = 100;
end
if nargin < 6 || isempty(radius)
    radius = Inf;
end
if nargin < 7
    center = [];
end

solution = struct();
state = kit_dispersion_state(plasma, p, z, x);
state.d = kit_dispersion(state, 2);
if state.d.IERR
    % The C engine only raises the heavily-damped flag here and carries on
    % with the D and DX that the previous point left in those globals.
    % Stopping is the deterministic reading of the same intent.
    status = 2;
    info = struct('steps',0,'clamped',0,'displacement',0,'residual',inf);
    return
end

[state, rootFound, status, info] = ...
    kit_newtonRoot(plasma, state, maxIterations, radius, center);
if ~rootFound
    status = 0;
    return
end

% A converged root: differentiate once more for the group velocity.
d = kit_dispersion(state, 4);
x = state.X;
vg = [-real(d.DP / d.DX), -real(d.DZ / d.DX)];
sg = [0 0];
if vg(1) ~= 0
    sg(1) = imag(x) / vg(1);
end
if vg(2) ~= 0
    sg(2) = imag(x) / vg(2);
end
solution = struct('X',x,'EFL',d.EFL,'BFL',d.BFL,'VG',vg, ...
    'SG',sg,'u',d.u,'DP',d.DP,'DZ',d.DZ,'DX',d.DX);
status = 1;
end

%% ---------------- from +whampkit/species_arrays.m ----------------
function [mass, density, temperature, anisotropy, delta, drift] = kit_species_arrays(species)
%SPECIES_ARRAYS Expand the species structs into the engine's component arrays.
%   Converts to engine units (density m^-3, temperature keV) and pads every
%   array to the engine's ten component slots.  As in the C driver the unused
%   slots are zero, except that delta is one.
nSpecies = numel(species);
mass = zeros(1, 10);
density = zeros(1, 10);
temperature = zeros(1, 10);
anisotropy = zeros(10, 2);
delta = ones(1, 10);
drift = zeros(1, 10);
for j = 1:nSpecies
    mass(j) = species{j}.m;
    density(j) = species{j}.n * 1e6;
    temperature(j) = species{j}.t / 1e3;
    anisotropy(j, :) = [species{j}.a, species{j}.b];
    delta(j) = species{j}.d;
    drift(j) = species{j}.vd;
end
end

%% ---------------- from +whampkit/susceptibility.m ----------------
function result = kit_susceptibility(state, j, order, a)
plasma = state.plasma;
x = state.XX(j);
z = state.ZZ(j);
p = state.PP(j);
vd = plasma.drift(j);
ai = 1i * a;
if plasma.mass(j) == 0
    % Electrons carry the opposite sign of charge.
    ai = -ai;
end
al = 0.5 * p^2;          % perpendicular argument of the R-function
ala = a * al;
xr = abs(real(x) - z * vd);
abz = abs(z);
xi = imag(x);
result.IERR = false;
result.XSI = complex(zeros(6,4));

% Validity check for the pole representation, copied from CHI().  It is
% written as nested negations of the C "if (ok) break" chain so that a NaN
% frequency fails every test and falls through to the error exit.
if ~(xi >= -abz)
    rx = xr - fix(xr);          % distance to the nearest harmonic
    if rx > 0.5
        rx = 1 - rx;
    end
    if ~(xr > 1 + 5 * abz && xr > 1 + 5 * p) && ~(xi >= -0.6 * rx)
        result.IERR = true;
        return
    end
end

[B, C] = kit_pole_constants();
Ysi = result.XSI;
Ysi(1,1) = a;                   % start from the cold-plasma terms
Ysi(6,1) = a + 2*a*vd^2;
for l = 1:8
    bl = B(l);
    cl = C(l);
    if z < 0
        cl = -cl;               % poles are mirrored for propagation backwards
    end
    dl = cl + vd;               % pole position including the drift
    y = x - dl*z;
    bly = bl / y;
    rc = kit_r_function(y, ala);
    xy = 1 + a*z*cl/y;
    ps = xy * rc(1);
    psp = xy * rc(3);
    Ysi(1,1) = Ysi(1,1) + a*bl*y*ps;
    Ysi(2,1) = Ysi(2,1) + ai*bl*psp;
    Ysi(3,1) = Ysi(3,1) + a*p*bl*dl*ps;
    Ysi(4,1) = Ysi(4,1) + bly*psp;
    Ysi(5,1) = Ysi(5,1) - ai*p*bly*dl*psp;
    Ysi(6,1) = Ysi(6,1) + 2*a*bly*dl^2*(x-vd*z+al*ps);
    if order <= 1
        continue
    end
    % ---- d/dx ----
    psy = xy * rc(2);
    ppy = xy * rc(4);
    xy = x / y;
    dps = xy * (psy - a*z*cl/y*rc(1));
    dpp = xy * (ppy - a*z*cl/y*rc(3));
    Ysi(1,2) = Ysi(1,2) + a*bl*(y*dps+x*ps);
    Ysi(2,2) = Ysi(2,2) + ai*bl*dpp;
    Ysi(3,2) = Ysi(3,2) + a*p*bl*dl*dps;
    Ysi(4,2) = Ysi(4,2) + bly*(dpp-xy*psp);
    Ysi(5,2) = Ysi(5,2) - ai*p*bly*dl*(dpp-xy*psp);
    Ysi(6,2) = Ysi(6,2) + 2*a*bly*dl^2*(x+al*dps-xy*(x-vd*z+al*ps));
    if order <= 2
        continue
    end
    % ---- d/dz ----
    zy = z / y;
    dps = zy * (a*cl*xy*rc(1)-dl*psy);
    Ysi(5,2) = Ysi(5,2) - ai*p*bly*dl*(dpp-xy*psp);
    dpp = zy * (a*cl*xy*rc(3)-dl*ppy);
    zy = dl * zy;
    Ysi(1,3) = Ysi(1,3) + a*bl*y*(dps-zy*ps);
    Ysi(2,3) = Ysi(2,3) + ai*bl*dpp;
    Ysi(3,3) = Ysi(3,3) + a*p*bl*dl*dps;
    Ysi(4,3) = Ysi(4,3) + bly*(dpp+zy*psp);
    Ysi(5,3) = Ysi(5,3) - ai*p*bly*dl*(dpp+zy*psp);
    Ysi(6,3) = Ysi(6,3) + 2*a*bly*dl^2*(al*dps-vd*z+zy*(x-vd*z+al*ps));
    if order <= 3
        continue
    end
    % ---- d/dp^2 ----
    rc = kit_r_function(y-1, ala);
    dps = 2*(psp-ps);
    yratio = y/(y-1);
    dpp = 2*ala*(yratio^2*(1+a*z*cl/y)*rc(3)-psp)-y*dps;
    Ysi(1,4) = Ysi(1,4) + a*bl*y*dps;
    Ysi(2,4) = Ysi(2,4) + ai*bl*dpp;
    Ysi(3,4) = Ysi(3,4) + a*p*bl*dl*(dps+ps);
    Ysi(4,4) = Ysi(4,4) + bly*(2*psp+dpp);
    Ysi(5,4) = Ysi(5,4) - ai*p*bly*dl*(psp+dpp);
    Ysi(6,4) = Ysi(6,4) + 4*a*bly*dl^2*al*psp;
end
Ysi(4,:) = Ysi(1,:) - 2*a^2*al*Ysi(4,:);
result.XSI = Ysi;
end

%% ---------------- from +whampkit/sweep_origin.m ----------------
function origin = kit_sweep_origin(spec)
%SWEEP_ORIGIN First value of PLG/ZLG: PM1, or PM2 for a negative step.
if spec(3) < 0
    origin = spec(2);
else
    origin = spec(1);
end
end

%% ---------------- from +whampkit/sweep_state.m ----------------
function sweep = kit_sweep_state(kpSpec, kzSpec, PZL)
%SWEEP_STATE Initial state for the C engine's two-axis sweep walk.
%   PLG/ZLG are the current positions in the sweep coordinate (logarithmic
%   when PZL is set) and P/Z the physical wave numbers derived from them.
%   KV marks that no anchor has been established for the current line yet, so
%   the starting frequency has to come from the sweep origin.  The anchor
%   (VO) and previous-point (O) predictors start at zero, matching the file
%   scope they occupy in the C engine; a zero DOX on the very first point is
%   what turns a non-convergent start into the NaN rows the reference shows.
sweep.PLG = kit_sweep_origin(kpSpec);
sweep.ZLG = kit_sweep_origin(kzSpec);
sweep.PZL = PZL;
sweep.p = kit_physical_axis_value(sweep.PLG, PZL);
sweep.z = kit_physical_axis_value(sweep.ZLG, PZL);
sweep.KV = 1;
sweep.done = false;
sweep.XVO = 0; sweep.PVO = 0; sweep.ZVO = 0;
sweep.PLO = 0; sweep.ZLO = 0;
sweep.DOX = 0; sweep.DOP = 0; sweep.DOZ = 0;
sweep.XO = 0; sweep.PO = 0; sweep.ZO = 0;
sweep.DX = 0; sweep.DP = 0; sweep.DZ = 0;
end

%% ---------------- from +whampkit/wave_energy.m ----------------
function [EFL,BFL,u] = kit_wave_energy(E,u1,u3,u2,u12,u32)
%WAVE_ENERGY Polarisation vectors and wave energy density at a root.
%   Transcribes ENERGY() in energy.c.  E is the 6x4 array from kit_dispersion();
%   u1/u3/u2/u12/u32 are the products of the refractive-index components in
%   units of the electron thermal speed.
%
%   EFL is the normalised electric-field polarisation, BFL the corresponding
%   magnetic-field polarisation (c = 299792.458 km/s in the engine's units),
%   and u the wave energy density.
%
%   The polarisation is obtained from the matrix of cofactors of the wave
%   operator, as in Melrose, Plasma Astrophysics Vol I p.43.
% d(omega) of the wave operator, symmetrised by the antisymmetric rules.
DDER = complex(zeros(3,3));
DDER(1,1) = E(1,1)+E(1,2)+u32;
DDER(1,2) = E(2,1)+E(2,2);
DDER(1,3) = E(3,1)+E(3,2)-u1*u3;
DDER(2,2) = E(4,1)+E(4,2)+u2;
DDER(2,3) = E(5,1)+E(5,2);
DDER(3,3) = E(6,1)+E(6,2)+u12;
DDER(2,1) = -DDER(1,2);
DDER(3,1) = DDER(1,3);
DDER(3,2) = -DDER(2,3);

% The wave operator itself: the dielectric tensor with the vacuum part
% replaced by the refractive-index terms.
EE = E(:,1);
EE(1) = EE(1)-u32;
EE(3) = EE(3)+u1*u3;
EE(4) = EE(4)-u2;
EE(6) = EE(6)-u12;

CM = complex(zeros(3,3));
CM(1,1) = EE(4)*EE(6)+EE(5)^2;
CM(1,2) = EE(2)*EE(6)+EE(5)*EE(3);
CM(1,3) = EE(2)*EE(5)-EE(4)*EE(3);
CM(2,2) = EE(1)*EE(6)-EE(3)^2;
CM(2,3) = EE(5)*EE(1)+EE(3)*EE(2);
CM(3,3) = EE(1)*EE(4)+EE(2)^2;
CM(2,1) = -CM(1,2);
CM(3,1) = CM(1,3);
CM(3,2) = -CM(2,3);

% Summing each column of the cofactor matrix gives a vector proportional to
% the polarisation, whichever row of the (singular) operator is degenerate.
efl = complex(zeros(1,3));
for i = 1:3
    for j = 1:3
        efl(i) = efl(i) + CM(j,i);
    end
end
% Normalise against whichever of Ex or Ez is better determined, so that
% Im Ey < 0 keeps meaning an L-mode wave (correction of 1994-04-13).
if abs(efl(1)) < 1000*abs(efl(3))
    econ = conj(efl(3));
else
    econ = conj(efl(1));
end
A = efl(1)*econ;
B = efl(2)*econ;
C = efl(3)*econ;
q = sqrt(real(A*conj(A)+B*conj(B)+C*conj(C)));
EFL = [A B C] / q;
cLight = 1/299.792458;
BFL = [-cLight*u3*EFL(2), cLight*(u3*EFL(1)-u1*EFL(3)), cLight*u1*EFL(2)];

erg = 0;
for i = 1:3
    for j = 1:3
        erg = erg + conj(EFL(i))*DDER(i,j)*EFL(j);
    end
end
u = real(erg*conj(erg)/(erg+conj(erg)));
end
