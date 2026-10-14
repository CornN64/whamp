function n = windingNumber(plasma, p, z, x0, radius, numSamples)
%WINDINGNUMBER Number of zeros of D inside a circle, by the argument principle.
%   N = whampkit.windingNumber(PLASMA, P, Z, X0) counts the zeros of the
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
%   whampkit.windingNumber(..., RADIUS) uses RADIUS instead of the default
%   1e-4 * max(1, |X0|).  whampkit.windingNumber(..., RADIUS, NUMSAMPLES)
%   changes the contour resolution from the default 32.
%
%   Points where the susceptibility expansion is refused (IERR) are treated
%   as "unknown" and the result is NaN, rather than silently counting as a
%   non-root.

if nargin < 5 || isempty(radius)
    radius = 1e-4 * max(1, abs(x0));
end
if nargin < 6 || isempty(numSamples)
    numSamples = 32;
end

theta = (0:numSamples) / numSamples * 2 * pi;
values = nan(1, numSamples + 1);
for k = 1:numSamples + 1
    x = x0 + radius * complex(cos(theta(k)), sin(theta(k)));
    d = whampkit.evalD(plasma, p, z, x, 2);
    if d.IERR ~= 0 || ~isfinite(d.D) || d.D == 0
        n = NaN;
        return
    end
    values(k) = d.D;
end

n = sum(angle(values(2:end) ./ values(1:end-1))) / (2 * pi);
n = round(n * 2) / 2;      % half-integers only if the contour is too coarse
end
