function origin = sweep_origin(spec)
%SWEEP_ORIGIN First value of PLG/ZLG: PM1, or PM2 for a negative step.
if spec(3) < 0
    origin = spec(2);
else
    origin = spec(1);
end
end

