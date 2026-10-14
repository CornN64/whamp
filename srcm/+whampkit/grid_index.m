function index = grid_index(value, spec, count)
% Nearest grid index of a walked wave number, measured from where the sweep
% started (the engine measures from PM1, which runs off the end of its array
% for a negative step).  Clamped so a rounding slip cannot index out of
% range.
if spec(3) == 0
    index = 1;
    return
end
index = 1 + whampkit.fnint((value - whampkit.sweep_origin(spec)) / spec(3));
index = min(max(index, 1), count);
end

