function whamp_vs_srcc64()
%WHAMP_VS_SRCC64 Compare whamp.m against the srcc64 engine golden output.
%   Parses the "fOUT=" lines from tools_golden/golden64/engine_*.txt and
%   compares them with the frequencies whamp.m computes for the same
%   settings.  Reports the largest relative deviation per case.

root = fileparts(fileparts(mfilename('fullpath')));
golden = fullfile(root, 'tools_golden', 'golden64');

% The engine cases load Models/Ex1 and Models/Ex2, whose values differ from
% whamp.m's built-in defaults (electron temperature 1 eV, and B set so that
% the gyro frequency is exactly the value in the model file).
gyroToB = @(gyro) gyro / 0.0279928;

% Models/Ex1: n 1e6/1e6 m-3, T .01/.001 keV, d 1/1, a 5/1, b 0/0,
%             masses 16/0, vd 1/0, XC 2.8 kHz.
model = struct('B', gyroToB(2.8));
model.Species = { ...
    struct('m',16,'n',1,'t',10,'a',5,'b',0,'d',1,'vd',1), ...
    struct('m',0,'n',1,'t',1,'a',1,'b',0,'d',1,'vd',0)};

% Models/Ex2: same densities/temperatures, a 1/1, masses 0/1 (electron first),
%             vd 0/0, XC 0.008978 kHz.
ex2 = struct('B', gyroToB(0.008978));
ex2.Species = { ...
    struct('m',0,'n',1,'t',10,'a',1,'b',0,'d',1,'vd',0), ...
    struct('m',1,'n',1,'t',1,'a',1,'b',0,'d',1,'vd',0)};

cases = {
    'engine_doc', model, struct('fstart',0.1,'kperp',0,'kpar',0.0022, ...
        'varyKzFirst',1,'useLog',0,'maxIterations',50);
    'engine_sweep', model, struct('fstart',0.1,'kperp',[0 0.01 0.04], ...
        'kpar',0.0022,'varyKzFirst',1,'useLog',0,'maxIterations',50);
    'engine_log', model, struct('fstart',0.1,'kperp',[-2.3 0.1 -1.8], ...
        'kpar',[-3 0.1 -2.7],'varyKzFirst',1,'useLog',1,'maxIterations',50);
    'engine_change', model, struct('fstart',0.1,'kperp',0,'kpar',0.0022, ...
        'varyKzFirst',1,'useLog',0,'maxIterations',50);
    'engine_ex2', ex2, struct('fstart',0.1,'kperp',[0 0.01 0.02], ...
        'kpar',[0.001 0.0022 0.005],'varyKzFirst',1,'useLog',0, ...
        'maxIterations',50);
    };

fprintf('%-16s %10s %12s %12s  %s\n', 'case', 'points', 'max rel', 'max abs', 'verdict');
for k = 1:size(cases,1)
    name = cases{k,1};
    params = cases{k,3};
    path = fullfile(golden, [name '.txt']);
    if exist(path, 'file') ~= 2
        fprintf('%-16s MISSING golden file\n', name);
        continue
    end
    ref = read_fout(path);
    % engine_change re-runs after 't.02'; the last fOUT block is the one
    % matching the modified model, so compare against the first block, which
    % is the unmodified default plasma.
    if isempty(ref)
        fprintf('%-16s no fOUT block in golden file, skipped\n', name);
        continue
    end
    firstBlock = ref{1};
    ref = [firstBlock{:}];
    O = whamp(cases{k,2}, params);
    got = O.f.';            % engine stores [kperp x kpar]
    got = got(:).';
    refv = ref(:).';
    n = min(numel(refv), numel(got));
    if n ~= numel(refv)
        fprintf('%-16s point count %d vs golden %d\n', name, numel(got), numel(refv));
    end
    maxRel = 0;
    maxAbs = 0;
    for j = 1:n
        r = refv(j);
        g = got(j);
        if isnan(r) || isnan(g)
            if ~(isnan(r) && isnan(g))
                fprintf('%-16s NaN mismatch at %d\n', name, j);
            end
            continue
        end
        scale = max(abs(r), realmin);
        maxRel = max(maxRel, abs(g - r) / scale);
        maxAbs = max(maxAbs, abs(g - r));
    end
    if maxRel < 1e-9
        verdict = 'MATCH';
    elseif maxRel < 1e-4
        verdict = 'close';
    else
        verdict = 'DIFFERS';
    end
    fprintf('%-16s %10d %12.3e %12.3e  %s\n', name, n, maxRel, maxAbs, verdict);
end
end

function blocks = read_fout(path)
%READ_FOUT Extract the complex fOUT vectors from an engine test transcript.
txt = fileread(path);
blocks = {};
start = 1;
while true
    i = strpos(txt, ' fOUT=', start);
    if i == 0
        break
    end
    j = strpos(txt, newline, i);
    if j == 0
        j = numel(txt) + 1;
    end
    blocks{end+1} = parse_complex(txt(i+6 : j-1)); %#ok<AGROW>
    start = j + 1;
end
end

function v = parse_complex(s)
v = {};
i = 1;
n = numel(s);
while true
    o = strpos(s, '(', i);
    if o == 0
        break
    end
    c = strpos(s, ')', o);
    if c == 0
        break
    end
    parts = strsplit(strtrim(s(o+1:c-1)), ',');
    v{end+1} = complex(str2double(strtrim(parts{1})), ...
        str2double(strtrim(parts{2}))); %#ok<AGROW>
    i = c + 1;
end
end

function idx = strpos(hay, needle, from)
idx = strfind(hay, needle);
if isempty(idx)
    idx = 0;
else
    idx = idx(idx >= from);
    if isempty(idx)
        idx = 0;
    else
        idx = idx(1);
    end
end
end
