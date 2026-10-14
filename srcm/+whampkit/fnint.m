function n = fnint(x)
%FNINT Fortran NINT: nearest integer, ties away from zero.
if x >= 0
    n = floor(x + 0.5);
else
    n = ceil(x - 0.5);
end
end

