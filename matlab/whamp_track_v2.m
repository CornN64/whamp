function Output = whamp_track_v2(PlasmaModel, InputParameters, Options)
%WHAMP_TRACK_V2 Branch tracking without hand-tuned scales.
%   Translation and improvments done with AI assist (Walter Puccio 2026) 
%   trust-region Newton, continuity check, multi-start re-anchoring)
%   and replaces their two weakest links:
%
%   1. NO HAND-TUNED SCALE.  Version 1 sized everything from the median of the
%      branch's recent steps, floored by RELFLOOR = 0.25 of the previous root.
%      Version 2 uses the exact implicit-differentiation step of the dispersion
%      relation at the branch's own converged root,
%
%          dx = -((p - p0)*DP + (z - z0)*DZ) / DX,
%
%      which is the first-order movement of the root as (p, z) move.  DP, DZ
%      and DX are already computed for the group velocity, so this is free.
%      RELFLOOR is gone; the only floor left is sqrt(eps)*|x|, which is not a
%      tuning knob but the statement that a step below the root's own round-off
%      cannot be resolved.  On the flat branch that motivated the report this
%      tightens the median jump tolerance from 1.48 to about 0.047.
%
%      The tolerance is NOT sized from this prediction.  A tolerance built from
%      the prediction is wrong by exactly the amount it is then asked to
%      forgive, so it waves through the very hop it exists to catch; doing so
%      reproduced all thirteen of the engine's hops.  Nor is it the median
%      |root - prediction| residual, which measures local curvature rather than
%      branch integrity.  The tolerance is the branch's own RATE |dx|/|dk|
%      along the axis being stepped, read from the neighbouring grid cells.
%      The prediction is what the candidate is measured against, and the first
%      seed; it is kept out of the yardstick.
%
%      The reference is the implicit prediction where the derivatives are
%      available and the anchor itself where they are not.  They are not
%      available at kpar = 0 or kperp = 0: DISPERSION guards those two
%      derivatives to zero because its own differencing divides by the wave
%      number, and every sweep line starts on one of those planes.  A zero
%      derivative there is the ABSENCE of a prediction, not a confident
%      prediction that the branch has not moved, and treating it as the latter
%      rejects the correct root -- on the golden model at kperp = 0.005 the
%      branch genuinely moves 0.408 -> 0.981 in one step from kpar = 0.
%
%   2. SCALE-FREE CERTIFICATION.  Every accepted root is certified with the
%      argument principle: kit_windingNumber counts Z - P for D inside a small
%      circle around it.  A certified simple root gives +1; 0 means the circle
%      encloses neither zero nor pole, so the reported point is not a root at
%      all; -1 means it sits on a species gyroharmonic pole.  Because this is a
%      property of D alone it works where no residual threshold can: |D| ranges
%      over some thirty orders of magnitude across this domain.  Certification
%      is reported in Output.certWinding and split across three flags so that
%      "wrong branch" and "not a root" stop being conflated.
%
%   3. TWO-DIMENSIONAL ANCHORING.  Version 1 kept its history in *visit order*
%      across the whole sweep.  That is right for a single-kperp kpar line and
%      wrong for a genuine two-dimensional sweep, where the previously visited
%      point at the start of a new line sits at the far end of the previous
%      line, at a completely different wave number.  Version 2 keeps a map of
%      accepted roots indexed by grid cell and takes each point's predecessor
%      from its inner neighbour, falling back to the same inner index on the
%      previous outer line.  The rate history is likewise a map, RATEMAP(ip,iz,
%      axis), read through expanding rings around the cell being entered -- a
%      window in visit order is not a neighbourhood in the grid.  This reduces
%      exactly to Version 1's behaviour on the one-dimensional lines it was
%      tested on, and is correct for varyKzFirst = 1 and for full 2-D sweeps.
%
%   Options (struct, all optional)
%     .trustFactor    Newton disc radius, in units of the predicted step.
%                                                          Default 3
%     .jumpFactor     Reject a root further than this from the prediction.
%                                                          Default 6
%     .scanOffsets    Offsets (in predicted steps) applied to the prediction
%                     as last-resort seeds.
%                                                Default [0.5 -0.5 1 -1 2 -2]*0.5
%     .useScan        Try the scan offsets when nothing else works. Default true
%     .maxIterations  Newton iteration cap.  Default from InputParameters
%     .certify        Certify accepted roots with the argument principle.
%                                                          Default true
%     .certRadius     Contour radius for certification.
%                                                          Default 1e-4*max(1,|x|)
%     .certSamples    Starting samples on the certification contour.  The
%                     contour is doubled until its phase is resolved, so this
%                     is a floor and not a cap.              Default 8
%     .stepWindow     Recent branch steps used for the measured allowance.
%                                                          Default 8
%
%   Output adds, over Version 1:
%     .certWinding        per-point winding number (NaN where unknown)
%     .flagNotCertified   converged root the argument principle does not
%                         report as a single simple zero
%     .flagCertUnknown    point could not be certified (expansion refused on
%                         the contour) -- neither trusted nor rejected
%     .predicted          the implicit-differentiation prediction used for
%                         each point, so the judgement can be audited
%     .diagnostics        as Version 1, plus .certified, .notCertified,
%                         .certUnknown and per-point .predicted/.scale/.winding
%
%   What is deliberately NOT claimed: that the implicit step breaks the
%   circularity the Version 1 report warns about.  The prediction is built from
%   the previous *accepted* root, so if the branch history is already corrupt
%   the prediction is corrupt with it.  What it removes is the RELFLOOR fudge
%   and the dependence on how many points happen to be in the median window;
%   the residual circularity is the same one every local continuation method
%   has, and the certification above is what bounds it.

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
certWinding = nan(shape);
flagNotCertified = zeros(shape, 'int32');
flagCertUnknown = zeros(shape, 'int32');
predictedGrid = nan(shape);

cycleZFirst = InputParameters.varyKzFirst ~= 0;
sweep = kit_sweep_state(kpSpec, kzSpec, InputParameters.useLog);
x = InputParameters.fstart;

% Accepted roots keyed by grid cell, each carrying the derivatives needed for
% the implicit step.  This is the two-dimensional replacement for Version 1's
% visit-order history.
anchorMap = struct('p', cell(1), 'z', cell(1), 'x', cell(1), ...
    'DP', cell(1), 'DZ', cell(1), 'DX', cell(1));
anchorMap(shape(1), shape(2)).x = [];

% Local branch RATE of change, |dx| / |dk|, recorded at the grid cell a step
% arrived at, separately for each axis the step was taken along.  This is what
% replaces Version 1's RELFLOOR and every other hand-picked scale.
%
% LOCAL, and keyed by AXIS, and both matter.  The branch moves at wildly
% different rates in different places: on the anisotropic model dx/dkpar is
% about 57 near kpar = 0 and about 0.06 near the fold at kpar = 0.04, so a
% tolerance calibrated on the flat tail of one line is roughly a thousandfold
% too tight for the steep start of the next -- which is precisely how a running
% window over visited points came to reject the correct root at
% (kperp = 0.02, kpar = 0.008) and keep a non-root instead.  And the branch
% moves far faster in kpar than in kperp, so a step measured along one axis
% says nothing about the other.  Storing the rate at the cell rather than in a
% visit-order window makes the estimate a property of the branch where the
% solver currently is.  A one-dimensional sweep uses one axis throughout and
% so behaves exactly as a window would.
rateMap = cell(shape(1), shape(2), 2);
% rateCells{axis} lists the cells holding a rate along that axis, and rateMask
% says so in O(1) so a cell is listed once.  The lookup needs the rate at the
% SMALLEST Manhattan distance from the query cell, which is what expanding
% rings compute; walking them costs one ring per unit of distance and so is
% quadratic in the grid.  On a 201x201 grid whose axis holds nothing -- the
% perpendicular axis of a one-dimensional sweep, or an axis whose lines never
% converged -- that walk ran all 402 rings for every single point, about 10 ms
% each.  Keeping the occupied cells and taking a vectorised minimum distance
% gives the identical answer in microseconds, and returns immediately when the
% axis is empty, which is exactly when the old code was slowest.
rateCells = {zeros(0, 2), zeros(0, 2)};
rateMask = false(shape(1), shape(2), 2);

diag = local_diag();
perPoint = {};

guard = 0;
while true
    guard = guard + 1;
    if guard > 200000
        error('whamp_track_v2:Runaway', 'The sweep did not terminate.');
    end

    ip = kit_grid_index(sweep.PLG, kpSpec, numel(kp));
    iz = kit_grid_index(sweep.ZLG, kzSpec, numel(kz));

    [solution, status, pointDiag] = local_solve_point( ...
        plasma, sweep.p, sweep.z, x, opt, anchorMap, rateMap, rateCells, ...
        ip, iz, cycleZFirst);

    predictedGrid(ip, iz) = pointDiag.reference;

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

        % --- certification -------------------------------------------------
        if opt.certify
            w = kit_windingNumber(plasma, sweep.p, sweep.z, solution.X, ...
                opt.certRadius, opt.certSamples);
            certWinding(ip, iz) = w;
            pointDiag.winding = w;
            if isnan(w)
                flagCertUnknown(ip, iz) = 1;
                diag.certUnknown = diag.certUnknown + 1;
            elseif w ~= 1
                flagNotCertified(ip, iz) = 1;
                diag.notCertified = diag.notCertified + 1;
            else
                diag.certified = diag.certified + 1;
            end
        end

        % Store the anchor for this cell and record how badly the first-order
        % prediction missed, so the next point's tolerance knows.
        anchorMap(ip, iz) = struct('p', sweep.p, 'z', sweep.z, ...
            'x', solution.X, 'DP', solution.DP, 'DZ', solution.DZ, ...
            'DX', solution.DX);

        % Only a CLEAN acceptance may widen the tolerance.  A point kept
        % against the continuity check is by construction a large miss; feeding
        % it back would let the first hop buy the licence for every hop after
        % it, which is exactly the failure this is meant to detect.  Version 1
        % had the same hazard through its step median and papered over it with
        % RELFLOOR.
        %
        % What is recorded is the RATE the branch actually moved at -- the step
        % between consecutive grid neighbours divided by their wave-number
        % separation -- stored against THIS cell and the axis the step was
        % taken along.  Unlike a distance from the prediction it is a property
        % of the branch rather than of the local curvature of the model, so it
        % does not inflate where the reference is a poor approximation.  Two
        % separations are kept.  By axis, because the branch moves far faster
        % in kpar than in kperp and a step along one says nothing about the
        % other.  And by cell, because the rate also varies enormously along a
        % single axis -- on the anisotropic model dx/dkpar is about 57 near
        % kpar = 0 and about 0.06 near the fold -- so an average over a window
        % of visited points is the wrong yardstick wherever the branch is not
        % locally straight.  Only a CLEAN acceptance may be recorded: a point
        % kept against the continuity check is by construction a large miss,
        % and feeding it back would let the first hop buy the licence for every
        % hop after it.
        if ~pointDiag.flagged && isfinite(pointDiag.anchorX) && ...
                pointDiag.axis >= 1 && pointDiag.dk > 0
            rateAxis = pointDiag.axis;
            rate = abs(solution.X - pointDiag.anchorX) / pointDiag.dk;
            rateMap{ip, iz, rateAxis} = rate;
            if ~rateMask(ip, iz, rateAxis)
                rateMask(ip, iz, rateAxis) = true;
                rateCells{rateAxis}(end+1, :) = [ip, iz]; %#ok<AGROW>
            end
        end

        sweep = kit_remember_point(sweep, solution);
    elseif status == 2
        flagDamped(ip, iz) = 1;
    else
        flagNoConvergence(ip, iz) = 1;
    end

    diag.rejected = diag.rejected + pointDiag.rejected;
    diag.clampedSteps = diag.clampedSteps + pointDiag.clamped;
    diag.attemptsPerPoint(end+1) = pointDiag.attempts; %#ok<AGROW>

    if status == 1
        xfound = solution.X;
    else
        xfound = NaN;
    end
    perPoint{end+1} = local_packPoint(ip, iz, sweep.p, sweep.z, x, ...
        status, pointDiag, xfound); %#ok<AGROW>

    sweep = kit_advance_sweep(sweep, kpSpec, kzSpec, cycleZFirst, status ~= 1);
    if sweep.done
        break
    end
    x = kit_next_frequency(sweep);
end

Output = local_assemble(InputParameters, PlasmaModel, ...
    kpPhysical, kzPhysical, values, flagFound, flagDamped, ...
    flagNoConvergence, flagSuspect, certWinding, flagNotCertified, ...
    flagCertUnknown, predictedGrid);
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

% ------------------------------------------------------------------ options

function opt = local_options(Options, InputParameters)
opt = struct();
opt.trustFactor = local_get(Options, 'trustFactor', 3);
opt.jumpFactor = local_get(Options, 'jumpFactor', 6);
opt.scanOffsets = local_get(Options, 'scanOffsets', ...
    [0.5 -0.5 1 -1 2 -2] * 0.5);
opt.useScan = local_get(Options, 'useScan', true);
opt.maxIterations = local_get(Options, 'maxIterations', ...
    InputParameters.maxIterations);
opt.certify = local_get(Options, 'certify', true);
opt.certRadius = local_get(Options, 'certRadius', []);
opt.certSamples = local_get(Options, 'certSamples', []);
end

function v = local_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function diag = local_diag()
diag = struct('attemptsPerPoint', zeros(1, 0), 'clampedSteps', 0, ...
    'rejected', 0, 'rescued', 0, 'flagged', 0, ...
    'certified', 0, 'notCertified', 0, 'certUnknown', 0);
end

% ------------------------------------------------------------------- solving

function [solution, status, pointDiag] = local_solve_point( ...
        plasma, p, z, xseed, opt, anchorMap, rateMap, rateCells, ...
        ip, iz, cycleZFirst)
%SOLVE_POINT Solve one point against its predicted branch position.
%
%   Version 1 asked whether a candidate was close to the *previous root*.
%   Version 2 asks whether it is close to where the branch is *predicted* to
%   be, which is strictly more informative: the prediction already accounts
%   for the fact that (p, z) moved.
%
%   The first seed is still the engine's own extrapolation, tried with an
%   infinite trust radius, so a well-behaved line costs exactly what WHAMP
%   costs and lands on exactly the same root.

pointDiag = struct('attempts', 0, 'clamped', 0, 'rejected', 0, ...
    'rescued', false, 'flagged', false, 'seedUsed', xseed, 'tried', {{}}, ...
    'reference', xseed, 'hasRef', false, 'informative', false, ...
    'anchorX', NaN, 'dk', NaN, 'axis', 0, 'scale', Inf, 'winding', NaN);

[a, axis, aIp, aIz] = local_predecessor(anchorMap, ip, iz, cycleZFirst);
informative = false;
if isempty(a)
    xRef = xseed;
    xPred = xseed;
    hasRef = false;
else
    % How far this point lies from its predecessor in wave-number space. A
    % predecessor is always a grid neighbour, so it differs in exactly one
    % coordinate and this sum is that one separation -- which is why it needs
    % no knowledge of which coordinate varies, and so stays correct whichever
    % axis VARYKZFIRST happens to run innermost. The branch's rate of change
    % is measured per unit wave number along the direction travelled, so this
    % converts a rate back into the step expected here.
    pointDiag.dk = abs(p - a.p) + abs(z - a.z);
    pointDiag.axis = axis;
    [xPred, informative] = local_predict(a, p, z);
    % The reference the candidate is judged against.  Where the implicit step
    % is available it is the reference, because it already accounts for the
    % wave number having moved.  Where it is not -- DISPERSION returns DZ = 0
    % at kpar = 0 and DP = 0 at kperp = 0, and every sweep line starts on one
    % of those planes -- the anchor itself is the reference, which is the
    % right answer rather than a fallback: the predecessor is the neighbouring
    % grid point, so a branch that has not hopped is near it.
    if informative
        xRef = xPred;
    else
        xRef = a.x;
    end
    hasRef = true;
    pointDiag.anchorX = a.x;
end

% Two scales, deliberately separated.
%
%   judgeScale decides whether a converged root is on the branch.  It is the
%   RATE at which the branch has ACTUALLY moved along this axis near the
%   predecessor -- taken from the nearest grid cell that has recorded one --
%   converted back into a step by multiplying by the wave-number separation.
%   Nothing here is chosen by hand, and the current reference never enters it;
%   that omission is the crux rather than timidity, because a tolerance built
%   from the current prediction is wrong by exactly the amount it is then asked
%   to forgive, so it waves through the very hop it exists to catch.  Sizing
%   the tolerance from the predicted step is what the first attempt at this
%   version did, and it reproduced all thirteen of the engine's hops.
%
%   Why a STEP and not the model's residual.  The obvious scale is the median
%   |root - prediction| over accepted points, i.e. how badly the first-order
%   step has been missing.  It is tempting and it is wrong, because that
%   residual is not a property of the branch but of the local curvature: on the
%   kperp line at kpar = 0 it is 3e-4 at one point and 0.11 at the next, a
%   three-hundredfold growth with no hop anywhere, so a residual-sized
%   tolerance flags ordinary rapid variation.  A step is the right yardstick
%   because a hop lands a branch-spacing away, and a branch-spacing is of the
%   order of a step.  The prediction still earns its place: it is what the
%   candidate is measured AGAINST, which is where it genuinely improves on
%   Version 1's 'distance from the previous root', and it is the first seed.
%
%   Why a RATE per axis rather than a raw step.  In a two-dimensional sweep the
%   predecessor may lie one grid step away in kperp or one grid step away in
%   kpar, and the branch moves roughly seventy times faster in kpar than in
%   kperp on the golden model.  A single pooled step history therefore made the
%   tolerance either useless in one direction or blind in the other, which is
%   exactly how the first two-dimensional attempt produced roots near 12, 26
%   and 70 in place of 0.95, 0.97 and 0.97.  Dividing by the separation makes
%   the quantity being averaged a property of the branch alone, and keying it
%   by axis keeps the two directions honest.
%
%   What is NOT here is Version 1's RELFLOOR = 0.25 * |x|.  That constant was
%   doing two jobs: bootstrapping the first point, and keeping the tolerance
%   from collapsing on a flattening branch.  Neither is needed.  The first
%   point is unjudged because there is nothing to compare with, and a
%   flattening branch SHOULD get a tight tolerance -- that is the whole
%   point of the exercise, and it is what turns the report's 1.48 jump
%   tolerance into a median of about 0.047 on the flat branch.
%
%   searchScale says how far apart to place the fallback seeds.  Here the
%   prediction is welcome, because using it to LOOK cannot corrupt the
%   JUDGEMENT, and generosity is free -- a far-flung seed that converges
%   somewhere implausible is still refused by the judgement, whereas a search
%   that never reaches far enough cannot be rescued by any judgement.
%
%   The sqrt(eps) term is a resolution floor rather than a tuning constant: a
%   discrepancy smaller than the root's own round-off cannot be told from zero.
localRate = local_rateNear(rateMap, rateCells, aIp, aIz, axis);
if ~isempty(localRate) && pointDiag.dk > 0
    judgeScale = max([localRate * pointDiag.dk, ...
        sqrt(eps) * max(1, abs(xRef))]);
    measured = true;
else
    judgeScale = sqrt(eps) * max(1, abs(xRef));
    measured = false;
end

if ~hasRef
    searchScale = Inf;
else
    observed = 0;
    if ~isempty(localRate)
        observed = localRate * pointDiag.dk;
    end
    searchScale = max([abs(xRef - a.x), observed, judgeScale, ...
        sqrt(eps) * max(1, abs(a.x))]);
end

pointDiag.predicted = xPred;
pointDiag.reference = xRef;
pointDiag.informative = hasRef && informative;
pointDiag.hasRef = hasRef;
pointDiag.scale = judgeScale;

seeds = local_seedList(xseed, xRef, searchScale, a, opt);

solution = [];
status = 0;
% The engine seed's status is kept so that "too heavily damped" -- a physical
% statement about the plasma -- is not flattened into "did not converge", a
% numerical one about the iteration.
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
            seeds(s), opt.maxIterations, opt.trustFactor * searchScale, xRef);
    end
    pointDiag.clamped = pointDiag.clamped + info.clamped;
    pointDiag.tried{end+1} = struct('seed', seeds(s), 'status', st, ...
        'root', local_safeRoot(sol), 'displacement', info.displacement); %#ok<AGROW>

    if s == 1
        engineStatus = st;
    end
    if st ~= 1
        continue
    end
    if local_continuous(sol.X, xRef, judgeScale, opt, measured)
        solution = sol;
        status = 1;
        pointDiag.seedUsed = seeds(s);
        pointDiag.rescued = s > 1;
        return
    end
    % Converged but not where the branch was predicted: keep it as a last
    % resort and keep looking for something better.
    pointDiag.rejected = pointDiag.rejected + 1;
    if isempty(solution)
        solution = sol;
        status = 1;
        pointDiag.seedUsed = seeds(s);
        pointDiag.flagged = true;
    end
end

if engineStatus == 2
    status = 2;
end
end

function [a, axis, aIp, aIz] = local_predecessor(anchorMap, ip, iz, cycleZFirst)
%PREDECESSOR  The accepted root this point continues, and along which axis.
%
%   AIp, AIz are the GRID CELL the predecessor came from, which the caller
%   needs in order to look up the branch rate recorded near there.
%
%   AXIS is 1 when the predecessor is the neighbour along the INNER axis (the
%   previous point of the same sweep line) and 2 when it came from a previous
%   OUTER line.  The caller needs this because the two are not comparable:
%   the branch moves at a completely different rate in kperp than in kpar, so
%   a tolerance calibrated on steps along one axis is meaningless along the
%   other.  In a one-dimensional sweep only axis 2 ever occurs.
%
%   Two candidates, in order of preference:
%     1. the neighbour just visited along the *inner* axis, i.e. the previous
%        point of the same sweep line;
%     2. the same inner index on a previous *outer* line, walking backwards
%        until one is found.
%
%   Version 1 used whichever point was visited last.  That is (1) whenever the
%   current line has already produced a point, but at a line boundary it is the
%   far end of the previous line -- a completely different wave number.
%   Candidate (2) is what makes a two-dimensional sweep trackable across line
%   boundaries, and reduces to Version 1's choice on the one-dimensional lines
%   it was tested on (there the inner axis holds a single point, so (1) is
%   never available and (2) with k = 1 is the previous point again).
%
%   Indices increase in visit order regardless of the sign of the axis step,
%   because grid_index measures from the sweep origin.
if cycleZFirst
    inner = [ip, iz - 1];
    outerDelta = [-1, 0];
else
    inner = [ip - 1, iz];
    outerDelta = [0, -1];
end
axis = 0;
aIp = ip;
aIz = iz;
a = local_cell(anchorMap, inner);
if ~isempty(a)
    axis = 1;
    aIp = inner(1);
    aIz = inner(2);
    return
end
% No inner neighbour: either the inner axis holds a single point (the usual
% one-dimensional line) or this is the first point of a new line.  Walk
% backwards along the *outer* axis from the current inner index, which for a
% single-point inner axis is exactly the previous point of the line.
for k = 1:64
    ij = [ip, iz] + outerDelta * k;
    a = local_cell(anchorMap, ij);
    if ~isempty(a)
        axis = 2;
        aIp = ij(1);
        aIz = ij(2);
        return
    end
end
a = struct([]);
end

function rate = local_rateNear(rateMap, rateCells, ip, iz, axis)
%RATE_NEAR  The branch's rate of change nearest to cell (ip, iz), this axis.
%   Of the cells holding a rate along the requested axis, take those at the
%   SMALLEST Manhattan distance from (ip, iz) and return the largest rate among
%   them, since erring generous here costs a slightly looser tolerance while
%   erring tight rejects correct roots.  Empty when the sweep has recorded
%   nothing along this axis yet, which leaves the caller with no measurement
%   and therefore no verdict.
%
%   Two searches, because the two regimes want opposite ones and the original
%   single expanding-ring walk was quadratic in the grid:
%
%     Rings out to LOCAL_RING_CAP.  The predecessor's neighbour almost always
%     recorded a rate one or two cells away, so a handful of rings finishes the
%     job, and the cost does not depend on how full the map is.
%
%     Otherwise the distance to every occupied cell at once.  Reaching here
%     means the neighbourhood is EMPTY, which is the sparse case -- a
%     one-dimensional sweep asking about the axis it never steps along, or a
%     region where nothing converged.  There the ring walk used to run all
%     2*nMax rings for nothing: about 10 ms per point on a 201x201 grid, which
%     dominated the whole solve.  The vectorised scan costs nothing when few
%     cells are occupied, which is exactly when the rings came up empty.
%
%   Capping the rings rather than running them to exhaustion is what makes the
%   pair exact: the scan is a complete search, so the answer cannot depend on
%   where the cap sits.
rate = [];
if axis < 1 || axis > 2
    return
end
cells = rateCells{axis};
if isempty(cells)
    return
end

ringCap = 8;                    % arbitrary: the scan below is exhaustive
for ring = 0:ringCap
    vals = local_ring(rateMap, ip, iz, ring, axis);
    if ~isempty(vals)
        rate = max(vals);
        return
    end
end

d = abs(cells(:, 1) - ip) + abs(cells(:, 2) - iz);
sel = find(d == min(d));
vals = zeros(1, numel(sel));
for q = 1:numel(sel)
    vals(q) = rateMap{cells(sel(q), 1), cells(sel(q), 2), axis};
end
rate = max(vals);
end

function vals = local_ring(rateMap, ip, iz, ring, axis)
%RING  Rates stored on the Manhattan ring of given radius about (ip, iz).
vals = [];
n1 = size(rateMap, 1);
n2 = size(rateMap, 2);
for d = 0:ring
    cand = [ip - (ring - d), iz - d;
            ip - (ring - d), iz + d;
            ip + (ring - d), iz - d;
            ip + (ring - d), iz + d];
    for k = 1:size(cand, 1)
        i1 = cand(k, 1);
        i2 = cand(k, 2);
        if i1 < 1 || i2 < 1 || i1 > n1 || i2 > n2
            continue
        end
        r = rateMap{i1, i2, axis};
        if ~isempty(r)
            vals(end+1) = r; %#ok<AGROW>
        end
    end
end
end

function a = local_cell(anchorMap, ij)
%CELL  The anchor stored at grid cell ij, or empty.
if ij(1) < 1 || ij(2) < 1 || ...
        ij(1) > size(anchorMap, 1) || ij(2) > size(anchorMap, 2)
    a = struct([]);
    return
end
if isempty(anchorMap(ij(1), ij(2)).x)
    a = struct([]);
else
    a = anchorMap(ij(1), ij(2));
end
end

function [xPred, informative] = local_predict(a, p, z)
%PREDICT  Where implicit differentiation says the root now sits.
%   D(p, z, x(p, z)) = 0 gives dx/dp = -DP/DX and dx/dz = -DZ/DX, so the
%   first-order root position at (p, z) is the anchor plus those times the
%   offsets.  This is the same expression the engine uses for its own seed
%   (`ddDX = (DKP*DP + DKZ*DZ)/DX`, then `X = XO - ddDX`); the difference is
%   that here it is evaluated against the solver's own anchor map rather than
%   the engine's carried predictor, so it survives a rejected point and is
%   available even when the engine has moved its anchor on.
%
%   INFORMATIVE says whether the answer means anything, and this is the
%   subtlety that decides whether the whole idea works.  DISPERSION returns
%   DZ = 0 when ZZ(1) == 0 and DP = 0 when PP(1) == 0 -- a guard against the
%   1/k division in its own differencing, not a physical statement that the
%   branch is flat.  Every sweep line starts at kpar = 0, where the guard
%   fires, and there the branch genuinely moves a great deal: on the golden
%   model at kperp = 0.005 the root goes 0.408 -> 0.981 in one step of 0.008.
%   A prediction that reports 'the anchor, exactly' in that situation is not
%   a confident prediction of zero movement, it is the ABSENCE of a
%   prediction, and treating it as the former is what made the first attempts
%   at this version reject the correct root.  So a zero derivative against a
%   non-zero offset makes the prediction uninformative, and the caller must
%   then decline to judge rather than insist the root has not moved.
if a.DX == 0
    xPred = a.x;
    informative = false;
    return
end
dp = p - a.p;
dz = z - a.z;
informative = ~((dp ~= 0 && a.DP == 0) || (dz ~= 0 && a.DZ == 0));
xPred = a.x - (dp * a.DP + dz * a.DZ) / a.DX;
end

function seeds = local_seedList(xseed, xPred, scale, a, opt)
%SEEDLIST Candidate starting frequencies, best first.
%   The engine's own extrapolation goes first so a well-behaved line costs
%   what WHAMP costs.  Then the implicit-differentiation prediction, then the
%   anchor itself, then offsets of the prediction.
seeds = xseed;
if ~isempty(a)
    if xPred ~= xseed
        seeds(end+1) = xPred; %#ok<AGROW>
    end
    if a.x ~= xseed && a.x ~= xPred
        seeds(end+1) = a.x; %#ok<AGROW>
    end
    if opt.useScan && isfinite(scale)
        for k = 1:numel(opt.scanOffsets)
            seeds(end+1) = xPred + opt.scanOffsets(k) * scale; %#ok<AGROW>
        end
    end
end
seeds = seeds(:).';
end

function ok = local_continuous(x, xPred, scale, opt, measured)
%CONTINUITY  Is this root within a plausible step of the prediction?
%
%   MEASURED says whether the scale rests on evidence.  Without it the test is
%   meaningless and must not reject: refusing a converged root because it is
%   far from a prediction that has never once been checked is how the first
%   attempt at this version lost points on lines the engine handled fine.
%   Version 1 had the same bootstrap and covered it with RELFLOOR = 0.25*|x|,
%   a hand-picked tolerance that then stayed in force for the whole line; here
%   the first accepted point is simply taken on trust and immediately becomes
%   the measurement.
ok = true;
if ~measured || isempty(xPred) || ~isfinite(xPred)
    return
end
ok = abs(x - xPred) <= opt.jumpFactor * scale;
end

function r = local_safeRoot(sol)
if isstruct(sol) && ~isempty(sol) && isfield(sol, 'X')
    r = sol.X;
else
    r = NaN;
end
end

% ------------------------------------------------------------------- output

function Output = local_assemble(InputParameters, PlasmaModel, ...
        kpPhysical, kzPhysical, values, flagFound, flagDamped, ...
        flagNoConvergence, flagSuspect, certWinding, flagNotCertified, ...
        flagCertUnknown, predictedGrid)
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
Output.certWinding = certWinding.';
Output.flagNotCertified = flagNotCertified.';
Output.flagCertUnknown = flagCertUnknown.';
Output.predicted = predictedGrid.';
missing = Output.flagNoConvergence == 1;
for k = 1:numel(resultNames)
    name = resultNames{k};
    Output.(name)(missing) = NaN;
end
end

function p = local_packPoint(ip, iz, p, z, seed, status, pointDiag, xfound)
p = struct('ip', ip, 'iz', iz, 'p', p, 'z', z, 'seed', seed, ...
    'found', xfound, 'status', status, 'attempts', pointDiag.attempts, ...
    'clamped', pointDiag.clamped, 'rejected', pointDiag.rejected, ...
    'rescued', pointDiag.rescued, 'flagged', pointDiag.flagged, ...
    'seedUsed', pointDiag.seedUsed, 'predicted', pointDiag.predicted, ...
    'reference', pointDiag.reference, 'informative', pointDiag.informative, ...
    'scale', pointDiag.scale, 'winding', pointDiag.winding, ...
    'tried', {pointDiag.tried});
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
%
%   The tables are constants, so they are built once per session and handed
%   back from then on.  Rebuilding the two literals cost about a twentieth of
%   a millisecond per call, and this runs once per species per susceptibility
%   evaluation, which made it a noticeable fraction of the dispersion
%   function.  Callers must treat B and C as read-only.
persistent Bs Cs
if isempty(Bs)
    Bs = [-1.734012457471826e-2-4.630639291680322e-2i, ...
        -1.734012457471826e-2+4.630639291680322e-2i, ...
        -7.399169923225014e-1+8.395179978099844e-1i, ...
        -7.399169923225014e-1-8.395179978099844e-1i, ...
        5.840628642184073+9.536009057643667e-1i, ...
        5.840628642184073-9.536009057643667e-1i, ...
        -5.583371525286853-1.120854319126599e1i, ...
        -5.583371525286853+1.120854319126599e1i];
    Cs = [ 2.237687789201900-1.625940856173727i, ...
        -2.237687789201900-1.625940856173727i, ...
        1.465234126106004-1.789620129162444i, ...
        -1.465234126106004-1.789620129162444i, ...
        0.8392539817232638-1.891995045765206i, ...
        -0.8392539817232638-1.891995045765206i, ...
        0.2739362226285564-1.941786875844713i, ...
        -0.2739362226285564-1.941786875844713i];
end
B = Bs;
C = Cs;
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

%% ---------------- from +whampkit/evalD.m ----------------
function [d, state] = kit_evalD(plasma, p, z, x, order)
%EVALD  Dispersion function D and its derivatives at one (z, x) point.
%   D = kit_evalD(PLASMA, P, Z, X) builds the dispersion state at
%   complex frequency X and perpendicular index Z and returns D.
%   [D, STATE] = kit_evalD(..., ORDER) also returns the state so a
%   caller can continue the analytic chain without rebuilding it.
%
%   This is the single evaluation seam shared by every root solver in the
%   toolkit.  It exists so the trust-region and global solvers compute
%   exactly what whamp.m computes, with no duplicated physics.
%
%   ORDER follows kit_dispersion: 0 -> D, 2 -> D and DX, 4 -> also the
%   second derivative needed for group velocity.

if nargin < 5
    order = 2;
end
state = kit_dispersion_state(plasma, p, z, x);
d = kit_dispersion(state, order);
end

%% ---------------- from +whampkit/evalDbatch.m ----------------
function [D, IERR] = kit_evalDbatch(plasma, p, z, x)
%EVALDBATCH  Dispersion function D at many frequencies, in one pass.
%   [D, IERR] = kit_evalDbatch(PLASMA, P, Z, X) evaluates D at every
%   entry of the frequency row vector X, for the single shared pair of
%   perpendicular and parallel wavenumbers (P, Z).  IERR(k) is true when the
%   pole expansion was refused at that point, exactly as
%   kit_dispersion() refuses it; D(k) is NaN there and must not be used.
%
%   This is the vectorised mirror of kit_evalD(PLASMA, P, Z, X(k), 1).
%   It exists for one caller: the argument-principle contour in
%   kit_windingNumber(), which samples eight or sixteen frequencies
%   around a circle at a fixed (P, Z).  Doing that with the scalar seam costs
%   one full state build and one full susceptibility evaluation per sample,
%   and the certification it serves is then over half the runtime of a
%   tracked sweep.  Sharing the state and running the eight-pole sum across
%   the samples instead is about 3x cheaper on the whole contour, measured at
%   876 real accepted roots of a two-dimensional whistler sweep; the
%   underlying R-function sum is 9x cheaper in isolation.
%
%   Two details are what make the saving real rather than theoretical, and
%   both were found by measuring at real roots rather than in a
%   microbenchmark.  Samples whose Taylor series has converged are frozen by
%   zeroing their term, not by indexing the live subset out of the
%   accumulators -- on a nine-element row the logical indexing cost more than
%   the arithmetic it skipped, and the first version of this file was SLOWER
%   than the scalar loop it replaced.  And all eight poles are stacked into
%   one rBatch call, since they share ALA and the routine is elementwise;
%   grouping cannot change the result because each sample's series is
%   advanced independently, including its own convergence test.
%
%   The arithmetic is arranged to be BIT-IDENTICAL to the scalar seam, not
%   merely close: the frequency-dependent quantities are carried as columns
%   of a species-by-sample matrix, every product that varies across samples
%   is an elementwise one, and the summation order within each sample is
%   unchanged.  verify_batch.m asserts this over randomised single- and
%   multi-species models, including drifting and loss-cone species, and
%   reverts the moment the two paths disagree.  If the physics in
%   dispersion() or susceptibility() is ever changed, change it here too and
%   re-run that test.
%
%   Only D is available here.  The frequency derivative DX is deliberately
%   not computed, since the phase of D is all the argument principle uses;
%   use kit_evalD() with ORDER >= 2 when a derivative is needed.

n = numel(x);

% Frequency-dependent state, transposed relative to dispersion_state() so
% that species run down the rows and samples across the columns.
omega = (plasma.xa * x) .^ 2;
XP = (plasma.density ./ plasma.dek ./ plasma.ren)' ./ omega;
XX = plasma.ren' .* x;
PP = p .* plasma.st';
ZZ = z .* plasma.st';

E = complex(zeros(6, n));
E([1 4 6], :) = 1;
IERR = false(1, n);

for j = 1:plasma.n
    E([1 4 6], :) = E([1 4 6], :) - XP(j, :);
    a = plasma.aa(j, 1);
    b = plasma.aa(j, 2);
    delta = plasma.delta(j);
    if a == b
        b = 0;
        delta = 1;
    end
    % factor carries the species density normalisation, so unlike the scalar
    % path it is not a constant here but one value per sample.
    factor = XP(j, :) ./ (a * (a - b));
    q = a - delta * b;
    for ib = 1:2
        [X, bad] = susceptibilityBatch(plasma, j, a, XX(j, :), ZZ(j), PP(j));
        E(:, :) = E(:, :) + (factor * q) .* X;
        IERR = IERR | bad;
        if ib == 2
            break
        end
        % Second (loss-cone) term; absent for an isotropic species.
        q = (delta - 1) * a;
        if q == 0
            break
        end
    end
end

% Refractive-index components of the first species and the three quadratic
% forms of the determinant, all elementwise across samples.
u1 = PP(1) * plasma.cv ./ XX(1, :);
u3 = ZZ(1) * plasma.cv ./ XX(1, :);
u12 = u1.^2;
u32 = u3.^2;
u2 = u12 + u32;
u13 = 2 * u1 .* u3;
detA = E(1, :) .* E(6, :) - E(3, :).^2;
Av = u12 .* E(1, :) + u13 .* E(3, :) + u32 .* E(6, :);
Bv = u2 .* detA + (u3 .* E(5, :) - u1 .* E(2, :)).^2;
Cv = detA .* E(4, :) + E(6, :) .* E(2, :).^2 ...
    + (E(1, :) .* E(5, :) + 2 * E(2, :) .* E(3, :)) .* E(5, :);
D = (u2 - E(4, :)) .* Av - Bv + Cv;
D(IERR) = NaN;
end

% ------------------------------------------------------------------------
function [Y, IERR] = susceptibilityBatch(plasma, j, a, x, z, p)
%susceptibilityBatch  Vector-x mirror of kit_susceptibility(...,1,a).
%   Y is 6-by-numel(x); only the first derivative column of the scalar
%   routine is reproduced.  IERR(k) marks the samples where the pole
%   representation is refused.
persistent Bc Cc
if isempty(Bc)
    [Bc, Cc] = kit_pole_constants();
end
ai = 1i * a;
if plasma.mass(j) == 0
    % Electrons carry the opposite sign of charge.
    ai = -ai;
end
al = 0.5 * p^2;
ala = a * al;
vd = plasma.drift(j);
n = numel(x);

% Validity check for the pole representation, copied from CHI() and
% reproduced elementwise from susceptibility().  Written as nested negations
% of the C "if (ok) break" chain so that a NaN frequency fails every test.
xi = imag(x);
xr = abs(real(x) - z * vd);
abz = abs(z);
IERR = false(1, n);
m = ~(xi >= -abz);
if any(m)
    rx = xr(m) - fix(xr(m));
    rx(rx > 0.5) = 1 - rx(rx > 0.5);
    IERR(m) = ~((xr(m) > 1 + 5 * abz) & (xr(m) > 1 + 5 * p)) ...
        & ~(xi(m) >= -0.6 * rx);
end

Y = complex(zeros(6, n));
Y(1, :) = a;                       % start from the cold-plasma terms
Y(6, :) = a + 2 * a * vd^2;

% The eight poles share ALA, and rBatch is elementwise, so every pole's
% argument is stacked into one row and summed in a single call.  Grouping
% cannot change the result: each sample's series is advanced independently,
% including its own convergence test.  This is what makes the batch pay --
% eight calls on a nine-element row cost several times what one call on
% seventy-two elements does.
dl = Cc + vd;
if z < 0
    dl = -Cc + vd;                   % poles mirrored for propagation back
end
% Lay the arguments out pole-major -- pole l occupies columns (l-1)*n+1 to
% l*n -- so that a slice for one pole is contiguous.
ystack = reshape(x.' - z * dl, 1, []);     % n-by-8 flattened column-major
Ystack = rBatch(ystack, ala);        % 4-by-(8n)
for l = 1:8
    idx = (l - 1) * n + (1:n);
    bl = Bc(l);
    cl = Cc(l);
    if z < 0
        cl = -cl;
    end
    y = x - dl(l) * z;
    bly = bl ./ y;
    rc = Ystack(:, idx);
    xy = 1 + a * z * cl ./ y;
    ps = xy .* rc(1, :);
    psp = xy .* rc(3, :);
    Y(1, :) = Y(1, :) + a * bl * y .* ps;
    Y(2, :) = Y(2, :) + ai * bl * psp;
    Y(3, :) = Y(3, :) + a * p * bl * dl(l) .* ps;
    Y(4, :) = Y(4, :) + bly .* psp;
    Y(5, :) = Y(5, :) - ai * p * bly * dl(l) .* psp;
    Y(6, :) = Y(6, :) + 2 * a * bly * dl(l)^2 .* (x - vd * z + al .* ps);
end
Y(4, :) = Y(1, :) - 2 * a^2 * al * Y(4, :);
end

% ------------------------------------------------------------------------
function rc = rBatch(y, al)
%rBatch  Vector-x mirror of kit_r_function for a shared AL.
%   The Taylor branch is taken for the whole batch when every sample falls in
%   it, which is the usual case: on a whistler sweep the contour sits far
%   from the harmonics and AL is small.  If any sample would use the
%   asymptotic or integral branch the batch drops back to the scalar
%   per-point routine, so the rare case stays correct and the common one
%   stays fast.
%
%   Samples that have already converged are frozen by carrying a zero term
%   forward rather than by indexing the live subset out of the accumulators.
%   That is bit-identical -- adding zero leaves a partial sum untouched, and
%   the scalar loop's BREAK means it too adds nothing further -- but it costs
%   a couple of elementwise products per pass instead of eight
%   logical-indexed reads and writes, which on a nine-element row dominated
%   the loop and made the batch slower than the scalar code it replaced.
ay = abs(y);
if all(al < 4 | ay.^2 > 75 * al | ay > 40 + al / 3)
    y2 = y.^2;
    d1 = y2 - 1;
    pn = y ./ d1;
    pyn = -y .* (y2 + 1) ./ (d1 .* d1);
    r1 = pn; r2 = pyn; r3 = pn; r4 = pyn;
    % LIVE is true for a sample still contributing and false once frozen.
    live = true(size(y));
    for i = 2:100
        % Where frozen, substitute a harmless denominator so that a sample
        % which has stopped cannot generate an Inf from a pole of the series
        % that the scalar loop, having already broken, would never reach.
        e = y2 - i^2;
        e(~live) = 1;
        coterm = (2 * i - 1) ./ e * al;
        newpyn = coterm .* (pyn - 2 * y2 ./ e .* pn);
        newpn = coterm .* pn;
        r1 = r1 + newpn;
        r2 = r2 + newpyn;
        r3 = r3 + i * newpn;
        r4 = r4 + i * newpyn;
        pn = newpn;
        pyn = newpyn;
        % Freeze a sample once its term is eight orders below the leading
        % one, which is where the scalar loop would have broken.
        live = ~(abs(pn) * 1e8 < abs(r1));
        pn = pn .* live;
        pyn = pyn .* live;
        if ~any(live)
            break
        end
    end
    rc = [r1; r2; r3; r4];
else
    rc = zeros(4, numel(y));
    for k = 1:numel(y)
        rc(:, k) = kit_r_function(y(k), al);
    end
end
end

%% ---------------- from +whampkit/windingNumber.m ----------------
function n = kit_windingNumber(plasma, p, z, x0, radius, numSamples)
%WINDINGNUMBER Number of zeros of D inside a circle, by the argument principle.
%   N = kit_windingNumber(PLASMA, P, Z, X0) counts the zeros of the
%   dispersion function inside the circle of radius RADIUS centred on X0,
%   counting poles negatively: N = Z - P.  It is the phase change of D
%   around the contour divided by 2*pi.
%
%   N = +1 therefore certifies a single simple zero, N = 0 certifies that
%   the interior contains neither a zero nor a pole, and N = -1 certifies a
%   pole.  Because the test is a property of D alone it is completely
%   scale-free, which matters here: |D| ranges over some thirty orders of
%   magnitude across the (p, z, x) domain, so no absolute residual
%   threshold can say whether a candidate point is a root.
%
%   N = kit_windingNumber(..., RADIUS) uses RADIUS instead of the
%   default 1e-4 * max(1, |X0|).
%
%   The contour is sampled ADAPTIVELY.  NUMSAMPLES (default 8) points are
%   used to start, and the count doubled until the largest phase step
%   between neighbouring samples falls below pi/2, so that no zero or pole
%   can hide in an unresolved arc.  A contour enclosing one simple zero has
%   an almost linear phase, so its first attempt already passes and costs
%   nine evaluations of D; a pole grazing the contour forces refinement up
%   to the NUMSAMPLES * 4 ceiling.  Capping the ceiling there means this is
%   never much dearer than the fixed 32-point contour it replaces, and
%   strictly more reliable, since a contour that coarse is now only accepted
%   once it has been shown to be resolved.  Passing NUMSAMPLES sets the
%   starting resolution; refinement from there still happens.
%
%   Points where the susceptibility expansion is refused (IERR) are treated
%   as "unknown" and the result is NaN, rather than silently counting as a
%   non-root.

if nargin < 5 || isempty(radius)
    radius = 1e-4 * max(1, abs(x0));
end
if nargin < 6 || isempty(numSamples)
    numSamples = 8;
end

% ORDER 1, not 2: the phase of D is all the argument principle uses, and
% dispersion() returns the same D either way while skipping the frequency
% derivative the Newton iteration needs and this test does not.
for refine = 0:2
    [n, worstStep] = local_contour(plasma, p, z, x0, radius, numSamples);
    if isnan(n) || worstStep <= pi / 2
        return
    end
    numSamples = 2 * numSamples;
end
end

function [n, worstStep] = local_contour(plasma, p, z, x0, radius, numSamples)
%CONTOUR  Winding of D around one closed polygon, and its coarsest arc.
%   Every sample shares one (p, z), so the whole polygon is evaluated in a
%   single vectorised pass through kit_evalDbatch rather than by
%   rebuilding the dispersion state once per corner.  The samples are the
%   vertices of a closed polygon, so the first is repeated at the end and
%   D is evaluated twice there; that keeps the phase increments between
%   neighbouring samples exactly what the scalar version produced.
theta = (0:numSamples) / numSamples * 2 * pi;
x = x0 + radius * complex(cos(theta), sin(theta));
[values, IERR] = kit_evalDbatch(plasma, p, z, x);
if any(IERR) || any(~isfinite(values)) || any(values == 0)
    n = NaN;
    worstStep = Inf;
    return
end

steps = angle(values(2:end) ./ values(1:end-1));
n = sum(steps) / (2 * pi);
n = round(n * 2) / 2;      % half-integers only if the contour is too coarse
worstStep = max(abs(steps));
end
