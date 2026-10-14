function rc = r_function(y, al)
%R_FUNCTION Relativistic plasma dispersion function and its derivatives.
%   RC = [R, dR/dy, R', d/dy R'] evaluated at complex Y for the real
%   perpendicular argument AL, following the method selection of RYLA():
%   a Taylor expansion for small AL or large |Y|, an asymptotic expansion
%   for large AL with moderate |Y|, and numerical integration otherwise.
ay = abs(y);
if al < 4 || ay^2 > 75*al || ay > 40+al/3
    rc = whampkit.r_taylor(y, al);
elseif 3*(al-10) > ay && ay^2 < 15*al
    rc = whampkit.r_asymptotic(y, al);
else
    rc = whampkit.r_integral(y, al);
end
end

