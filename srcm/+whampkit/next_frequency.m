function x = next_frequency(sweep)
%NEXT_FREQUENCY Linear extrapolation of the frequency to the next point.
if sweep.KV ~= 0
    x = sweep.XVO - ((sweep.p - sweep.PVO) * sweep.DOP ...
        + (sweep.z - sweep.ZVO) * sweep.DOZ) / sweep.DOX;
else
    x = sweep.XO - ((sweep.p - sweep.PO) * sweep.DP ...
        + (sweep.z - sweep.ZO) * sweep.DZ) / sweep.DX;
end
end

