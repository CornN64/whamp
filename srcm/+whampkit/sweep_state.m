function sweep = sweep_state(kpSpec, kzSpec, PZL)
%SWEEP_STATE Initial state for the C engine's two-axis sweep walk.
%   PLG/ZLG are the current positions in the sweep coordinate (logarithmic
%   when PZL is set) and P/Z the physical wave numbers derived from them.
%   KV marks that no anchor has been established for the current line yet, so
%   the starting frequency has to come from the sweep origin.  The anchor
%   (VO) and previous-point (O) predictors start at zero, matching the file
%   scope they occupy in the C engine; a zero DOX on the very first point is
%   what turns a non-convergent start into the NaN rows the reference shows.
sweep.PLG = whampkit.sweep_origin(kpSpec);
sweep.ZLG = whampkit.sweep_origin(kzSpec);
sweep.PZL = PZL;
sweep.p = whampkit.physical_axis_value(sweep.PLG, PZL);
sweep.z = whampkit.physical_axis_value(sweep.ZLG, PZL);
sweep.KV = 1;
sweep.done = false;
sweep.XVO = 0; sweep.PVO = 0; sweep.ZVO = 0;
sweep.PLO = 0; sweep.ZLO = 0;
sweep.DOX = 0; sweep.DOP = 0; sweep.DOZ = 0;
sweep.XO = 0; sweep.PO = 0; sweep.ZO = 0;
sweep.DX = 0; sweep.DP = 0; sweep.DZ = 0;
end

