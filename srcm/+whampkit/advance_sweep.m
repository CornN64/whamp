function sweep = advance_sweep(sweep, kpSpec, kzSpec, cycleZFirst, stopCycling)
%ADVANCE_SWEEP Move to the next point of the two-axis sweep.
%   The outer axis is walked by repeated addition of its step in the sweep
%   coordinate, exactly as the engine does.  A point that did not converge,
%   or that was too heavily damped, abandons the rest of the *inner* axis;
%   the remaining outer points are still visited and the sweep ends once the
%   outer axis is exhausted.  Reaching a new outer point resets KV so that
%   the line is re-anchored.
sweepEnd = 1e99;
if stopCycling
    if cycleZFirst
        sweep.ZLG = sweepEnd;
    else
        sweep.PLG = sweepEnd;
    end
end
if ~cycleZFirst
    sweep.PLG = sweep.PLG + kpSpec(3);
    if sweep.PLG >= kpSpec(1) && sweep.PLG <= kpSpec(2)
        sweep.p = sweep.PLG + sweep.PZL * (10^sweep.PLG - sweep.PLG);
    else
        sweep.ZLG = sweep.ZLG + kzSpec(3);
        if sweep.ZLG < kzSpec(1) || sweep.ZLG > kzSpec(2)
            sweep.done = true;
            return
        end
        sweep.KV = 1;
        sweep.PLG = sweep.PLO;
        sweep.p = sweep.PVO;
        sweep.z = sweep.ZLG + sweep.PZL * (10^sweep.ZLG - sweep.ZLG);
    end
else
    sweep.ZLG = sweep.ZLG + kzSpec(3);
    if sweep.ZLG >= kzSpec(1) && sweep.ZLG <= kzSpec(2)
        sweep.z = sweep.ZLG + sweep.PZL * (10^sweep.ZLG - sweep.ZLG);
    else
        sweep.PLG = sweep.PLG + kpSpec(3);
        if sweep.PLG < kpSpec(1) || sweep.PLG > kpSpec(2)
            sweep.done = true;
            return
        end
        sweep.KV = 1;
        sweep.ZLG = sweep.ZLO;
        sweep.z = sweep.ZVO;
        sweep.p = sweep.PLG + sweep.PZL * (10^sweep.PLG - sweep.PLG);
    end
end
end

