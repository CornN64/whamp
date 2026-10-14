function spec = axis_spec(axis, name)
%AXIS_SPEC Reduce a kperp/kpar request to the engine's [start end step].
%   The caller writes [start step end]; the engine keeps the first and third
%   numbers as the bounds and the second as the step, so the triple returned
%   here is [PM1 PM2 PM3].  The sweep walks PM3 at a time from PM1, or from
%   PM2 when PM3 is negative.  A bare scalar means one point, and the engine
%   default step of 10 is kept for it.
if ~isnumeric(axis) || ~isreal(axis) || isempty(axis) || ...
        ~(isscalar(axis) || numel(axis) == 3) || any(~isfinite(axis(:)))
    error('whamp:InvalidWaveNumber', ...
        'InputParameters.%s must be a finite scalar or [start step end].', name);
end
if isscalar(axis)
    spec = [axis axis 10];
else
    spec = [axis(1), axis(3), axis(2)];
end
if spec(1) ~= spec(2) && spec(3) == 0
    error('whamp:InvalidWaveNumberStep', ...
        'The %s step must be nonzero unless start equals end.', name);
end
end

