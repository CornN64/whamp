function tf = is_finite_real_scalar(value)
%IS_FINITE_REAL_SCALAR True for a real, finite numeric scalar.
%   Written with short-circuit && so that isfinite is never handed a cell or
%   character array, which would raise instead of reporting false.
tf = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
end

