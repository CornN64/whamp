function [d, state] = evalD(plasma, p, z, x, order)
%EVALD  Dispersion function D and its derivatives at one (z, x) point.
%   D = whampkit.evalD(PLASMA, P, Z, X) builds the dispersion state at
%   complex frequency X and perpendicular index Z and returns D.
%   [D, STATE] = whampkit.evalD(..., ORDER) also returns the state so a
%   caller can continue the analytic chain without rebuilding it.
%
%   This is the single evaluation seam shared by every root solver in the
%   toolkit.  It exists so the trust-region and global solvers compute
%   exactly what whamp.m computes, with no duplicated physics.
%
%   ORDER follows whampkit.dispersion: 0 -> D, 2 -> D and DX, 4 -> also the
%   second derivative needed for group velocity.

if nargin < 5
    order = 2;
end
state = whampkit.dispersion_state(plasma, p, z, x);
d = whampkit.dispersion(state, order);
end
