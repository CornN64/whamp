function result = susceptibility(state, j, order, a)
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

[B, C] = whampkit.pole_constants();
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
    rc = whampkit.r_function(y, ala);
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
    rc = whampkit.r_function(y-1, ala);
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

