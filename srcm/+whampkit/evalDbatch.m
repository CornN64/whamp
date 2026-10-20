function [D, IERR] = evalDbatch(plasma, p, z, x)
%EVALDBATCH  Dispersion function D at many frequencies, in one pass.
%   [D, IERR] = whampkit.evalDbatch(PLASMA, P, Z, X) evaluates D at every
%   entry of the frequency row vector X, for the single shared pair of
%   perpendicular and parallel wavenumbers (P, Z).  IERR(k) is true when the
%   pole expansion was refused at that point, exactly as
%   whampkit.dispersion() refuses it; D(k) is NaN there and must not be used.
%
%   This is the vectorised mirror of whampkit.evalD(PLASMA, P, Z, X(k), 1).
%   It exists for one caller: the argument-principle contour in
%   whampkit.windingNumber(), which samples eight or sixteen frequencies
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
%   use whampkit.evalD() with ORDER >= 2 when a derivative is needed.

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
%susceptibilityBatch  Vector-x mirror of whampkit.susceptibility(...,1,a).
%   Y is 6-by-numel(x); only the first derivative column of the scalar
%   routine is reproduced.  IERR(k) marks the samples where the pole
%   representation is refused.
persistent Bc Cc
if isempty(Bc)
    [Bc, Cc] = whampkit.pole_constants();
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
%rBatch  Vector-x mirror of whampkit.r_function for a shared AL.
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
        rc(:, k) = whampkit.r_function(y(k), al);
    end
end
end
