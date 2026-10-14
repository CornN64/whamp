function [EFL,BFL,u] = wave_energy(E,u1,u3,u2,u12,u32)
%WAVE_ENERGY Polarisation vectors and wave energy density at a root.
%   Transcribes ENERGY() in energy.c.  E is the 6x4 array from whampkit.dispersion();
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

