function axis = fill_axis(spec)
%FILL_AXIS Wave numbers of the stored solutions, in sweep order.
%   Matches the engine's axis_size()/whampkit.fill_axis() pair: the axis begins where
%   the sweep begins and advances by the step, so axis(i) is the wave number
%   of the i-th stored solution.
if spec(1) == spec(2)
    axis = spec(1);
    return
end
count = 1 + floor(abs(spec(2) - spec(1)) / abs(spec(3)));
if count > 100000
    error('whamp:InvalidWaveNumberRange', 'Wave-number axis has an invalid size.');
end
axis = whampkit.sweep_origin(spec) + (0:count-1) * spec(3);
end

