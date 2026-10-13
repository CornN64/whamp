function whamp_vs_repl_golden()
%WHAMP_VS_REPL_GOLDEN Compare whamp.m against the srccopt64 REPL transcripts.
%   Parses the "kp kz fre fim" rows from tools_golden/golden64/*.txt and
%   compares them with whamp.m for the same model and wavenumber settings.
%   The transcripts print 7 significant digits, so agreement is checked at
%   that level.

root = fileparts(fileparts(mfilename('fullpath')));
golden = fullfile(root, 'tools_golden', 'golden64');

gyroToB = @(gyro) gyro / 0.0279928;
ex1 = struct('B', gyroToB(2.8));
ex1.Species = { ...
    struct('m',16,'n',1,'t',10,'a',5,'b',0,'d',1,'vd',1), ...
    struct('m',0,'n',1,'t',1,'a',1,'b',0,'d',1,'vd',0)};
ex2 = struct('B', gyroToB(0.008978));
ex2.Species = { ...
    struct('m',0,'n',1,'t',10,'a',1,'b',0,'d',1,'vd',0), ...
    struct('m',1,'n',1,'t',1,'a',1,'b',0,'d',1,'vd',0)};

% name, model, InputParameters, useLog
cases = {
    'doc_pzf', ex1, struct('fstart',0.1,'kperp',0,'kpar',0.0022, ...
        'varyKzFirst',1,'useLog',0), 0;
    'sweep_p', ex1, struct('fstart',0.1,'kperp',[0 0.01 0.05], ...
        'kpar',0.0022,'varyKzFirst',1,'useLog',0), 0;
    'sweep_neg_step', ex1, struct('fstart',0.1,'kperp',[0 -0.01 0.05], ...
        'kpar',0.0022,'varyKzFirst',1,'useLog',0), 0;
    'log_scale', ex1, struct('fstart',0.1,'kperp',[-2.3 0.1 -1.8], ...
        'kpar',[-3 0.1 -2.7],'varyKzFirst',1,'useLog',1), 1;
    'no_convergence', ex1, struct('fstart',0.0001,'kperp',0,'kpar',0.0022, ...
        'varyKzFirst',1,'useLog',0), 0;
    };
% 'ex2_sweep' is deliberately not compared: the reference lets a NaN
% frequency and the SG components leak between points through file scope, and
% whamp.m documents that it starts each point from a deterministic state
% instead.  See the comments in solve_point and whamp.

fprintf('%-16s %8s %12s  %s\n', 'case', 'rows', 'max rel', 'verdict');
for k = 1:size(cases,1)
    name = cases{k,1};
    path = fullfile(golden, [name '.txt']);
    if exist(path, 'file') ~= 2
        fprintf('%-16s MISSING golden file\n', name);
        continue
    end
    rows = read_rows(path);
    if isempty(rows)
        fprintf('%-16s no result rows, skipped\n', name);
        continue
    end
    O = whamp(cases{k,2}, cases{k,3});
    maxRel = 0;
    bad = 0;
    for r = 1:size(rows,1)
        kp = rows(r,1);
        kz = rows(r,2);
        refF = complex(rows(r,3), rows(r,4));
        gotF = lookup(O, kp, kz);
        if isnan(refF)
            if ~isnan(gotF)
                bad = bad + 1;
            end
            continue
        end
        if isnan(gotF)
            bad = bad + 1;
            continue
        end
        rel = abs(gotF - refF) / max(abs(refF), realmin);
        if rel > maxRel
            maxRel = rel;
        end
    end
    if maxRel < 2e-7 && bad == 0
        verdict = 'MATCH';
    elseif maxRel < 1e-4 && bad == 0
        verdict = 'close';
    else
        verdict = sprintf('DIFFERS (nan mismatches: %d)', bad);
    end
    fprintf('%-16s %8d %12.3e  %s\n', name, size(rows,1), maxRel, verdict);
end
end

function v = lookup(O, kp, kz)
%LOOKUP Pick the stored solution whose axes match a transcript row.
[~, ip] = min(abs(O.kperp(:).' - kp));
[~, iz] = min(abs(O.kpar(:).' - kz));
v = O.f(iz, ip);
end

function rows = read_rows(path)
%READ_ROWS Parse "  kp  kz  fre  fim" result rows out of a transcript.
lines = strsplit(fileread(path), newline);
rows = zeros(0, 4);
for k = 1:numel(lines)
    s = strtrim(lines{k});
    if isempty(s) || s(1) == '#'
        continue
    end
    parts = regexp(s, '\s+', 'split');
    if numel(parts) < 4
        continue
    end
    vals = zeros(1, 4);
    ok = true;
    for j = 1:4
        if strcmpi(parts{j}, 'NaN')
            vals(j) = NaN;
        else
            v = str2double(parts{j});
            if isnan(v) && ~strcmpi(parts{j}, 'NaN')
                ok = false;
                break
            end
            vals(j) = v;
        end
    end
    if ok
        rows(end+1,:) = vals; %#ok<AGROW>
    end
end
end
