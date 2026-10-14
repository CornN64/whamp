function dc = dcTerm(E, col)
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
