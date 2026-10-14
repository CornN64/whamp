function [state, rootFound, status, info] = newtonRoot(plasma, state, maxIterations, radius, center)
%NEWTONROOT Newton search for a dispersion root, with an optional trust region.
%   [STATE, ROOTFOUND, STATUS] = whampkit.newtonRoot(PLASMA, STATE, MAXIT)
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
        [correction, wasClamped, exhausted] = limitStep(correction, state.X, center, radius);
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

        state = evaluateAt(plasma, state, state.X - correction);
        if state.d.IERR
            status = 2;
            info = makeInfo(steps, clamped, state, seed);
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

info = makeInfo(steps, clamped, state, seed);
end

function [correction, wasClamped, exhausted] = limitStep(correction, x, seed, radius)
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

function state = evaluateAt(plasma, state, x)
%EVALUATEAT Move to a new frequency and recompute the dispersion function.
state.X = x;
omega = (x * plasma.xa)^2;
state.XP = plasma.density ./ plasma.dek ./ plasma.ren ./ omega;
state.XX = x .* plasma.ren;
state.d = whampkit.dispersion(state, 2);
end

function info = makeInfo(steps, clamped, state, seed)
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
