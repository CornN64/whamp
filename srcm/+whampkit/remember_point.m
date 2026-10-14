function sweep = remember_point(sweep, solution)
%REMEMBER_POINT Store a converged point as the previous-point predictor.
sweep.PO = sweep.p;
sweep.ZO = sweep.z;
sweep.XO = solution.X;
sweep.DP = solution.DP;
sweep.DZ = solution.DZ;
sweep.DX = solution.DX;
if sweep.KV ~= 0
    % First converged point of this line: every later point on it is
    % predicted from this anchor.
    sweep.XVO = solution.X;
    sweep.PVO = sweep.p;
    sweep.ZVO = sweep.z;
    sweep.PLO = sweep.PLG;
    sweep.ZLO = sweep.ZLG;
    sweep.DOX = sweep.DX;
    sweep.DOP = sweep.DP;
    sweep.DOZ = sweep.DZ;
    sweep.KV = 0;
end
end

