function [solution, status, info] = solvePoint(plasma, p, z, x, maxIterations, radius, center)
%SOLVEPOINT Solve the dispersion relation at one (p, z) pair.
%   [SOLUTION, STATUS] = whampkit.solvePoint(PLASMA, P, Z, X, MAXIT) is the
%   kit's counterpart of whamp.m's SOLVE_POINT and returns the same fields.
%   STATUS is 1 for a converged root, 2 for a point the susceptibility
%   reports as too heavily damped, and 0 when the iteration gave up.
%
%   [..., RADIUS] passes a Newton trust-region radius through to
%   whampkit.newtonRoot.  RADIUS = Inf (the default) reproduces the
%   reference iteration exactly.
%
%   [..., INFO] returns the Newton diagnostics from whampkit.newtonRoot
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
state = whampkit.dispersion_state(plasma, p, z, x);
state.d = whampkit.dispersion(state, 2);
if state.d.IERR
    % The C engine only raises the heavily-damped flag here and carries on
    % with the D and DX that the previous point left in those globals.
    % Stopping is the deterministic reading of the same intent.
    status = 2;
    info = struct('steps',0,'clamped',0,'displacement',0,'residual',inf);
    return
end

[state, rootFound, status, info] = ...
    whampkit.newtonRoot(plasma, state, maxIterations, radius, center);
if ~rootFound
    status = 0;
    return
end

% A converged root: differentiate once more for the group velocity.
d = whampkit.dispersion(state, 4);
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
