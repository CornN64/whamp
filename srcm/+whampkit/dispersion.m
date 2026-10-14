function d = dispersion(state, order)
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
        xsi = whampkit.susceptibility(state, j, order, a);
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
dc = whampkit.dcTerm(E, 2);
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
    dc = whampkit.dcTerm(E, 3);
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
    dc = whampkit.dcTerm(E, 4);
    DP = ((u2-E(4,1))*da+(2*u12-E(4,4))*Av-db+dc)/state.PP(1);
end
[EFL,BFL,u] = whampkit.wave_energy(E,u1,u3,u2,u12,u32);
d = struct('D',D,'DX',DX,'DZ',DZ,'DP',DP,'EFL',EFL, ...
    'BFL',BFL,'u',u,'IERR',false);
end

