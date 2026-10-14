function value = physical_axis_value(axisValue, PZL)
%PHYSICAL_AXIS_VALUE Convert sweep coordinates to wave numbers.
if PZL == 1
    value = 10 .^ axisValue;
else
    value = axisValue;
end
end

