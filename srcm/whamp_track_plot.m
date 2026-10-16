function varargout = whamp_track_plot(Output, varargin)
%WHAMP_TRACK_PLOT  Dispersion diagram with branch confidence drawn in.
%   P = WHAMP_TRACK_PLOT(OUTPUT) plots the frequency branch from a
%   WHAMP_TRACK (or WHAMP_TRACK_V2) result and marks, on the curve
%   itself, how much each point can be trusted.  The solver's flags were
%   previously only countable, not visible: looking at a dispersion diagram it
%   is impossible to tell a confidently tracked root from one the continuity
%   check admitted only as a last resort.
%
%   Marks, in increasing order of concern:
%     solid line, filled dots   converged and certified (winding number +1)
%     open circle               converged but NOT certified by the argument
%                               principle -- the point is not a single simple
%                               zero of D (winding 0 or -1)
%     open square               could not be certified: the susceptibility
%                               expansion was refused on the contour.  Neither
%                               trust nor rejection; the solver declined to
%                               answer.
%     red X                     flagSuspectBranch: the root converged to
%                               something the continuity check rejected and
%                               was kept only so the line stays drawn.  This
%                               is the hop, and the mark to look for.
%     gap in the line           no solution.  Gaps are never bridged -- a
%                               straight segment across a hole is exactly the
%                               fiction the solver avoids.
%
%   WHAMP_TRACK_PLOT(OUTPUT,'ShowPredicted',true) overlays the
%   implicit-differentiation prediction (Output.predicted) as a faint dashed
%   line.  Where the solid curve leaves the dashed one is where the branch
%   stopped being a first-order continuation of its neighbour: a hop shows up
%   as a kink before it shows up as a jump.
%
%   WHAMP_TRACK_PLOT(OUTPUT,'Mode','plain') draws an ordinary plot with no
%   confidence marks.  'Frequency','f' selects the reported quantity, and any
%   other field of Output.works too (e.g. 'Frequency','EB').
%
%   Each column of OUTPUT is plotted as one line against Output.kpar.  For a
%   sweep whose fast axis is kperp, pass 'Transpose',true.
%
%   P = WHAMP_TRACK_PLOT(...) returns a struct of the graphics handles.
%
%   Also see: whamp_track.m, whamp_track_v2.m

% ---------------------------------------------------------------- arguments
opt = struct('Mode', 'flags', 'ShowPredicted', false, 'Frequency', 'f', ...
    'YScale', 'linear', 'MarkerSize', 6, 'LineWidth', 1.4, ...
    'Transpose', false);
k = 1;
while k <= numel(varargin)
    a = varargin{k};
    if ischar(a) && isfield(opt, a) && k < numel(varargin)
        opt.(a) = varargin{k+1};
        k = k + 2;
    elseif isstruct(a)
        fn = fieldnames(a);
        for j = 1:numel(fn)
            if isfield(opt, fn{j}), opt.(fn{j}) = a.(fn{j}); end
        end
        k = k + 1;
    else
        error('whamp_track_plot:BadArgument', ...
            'Expected ''Name'',Value pairs; got %s.', mat2str(a));
    end
end

if ~isstruct(Output) || ~isfield(Output, opt.Frequency)
    error('whamp_track_plot:BadOutput', ...
        'OUTPUT must be a whamp_track result struct with a %s field.', ...
        opt.Frequency);
end

% WHAMP_TRACK returns every result array as (kpar x kperp): the ROW index runs
% over kpar and the COLUMN index over kperp.  A dispersion curve is a line at
% FIXED kperp, so it is a COLUMN of the raw array and its abscissa is kpar --
% the raw orientation is already 'one curve per column', which is what the
% plotting loop below wants.  'Transpose',true asks for the other reading:
% curves at fixed kpar, drawn against kperp, which is the raw array transposed
% so that its columns are again the curves.
f = Output.(opt.Frequency);
if opt.Transpose
    f = f.';
    kx = real(Output.kperp(:));
else
    kx = real(Output.kpar(:));
end
[nrow, ncol] = size(f);
if numel(kx) ~= nrow
    error('whamp_track_plot:Shape', ...
        'the %s axis has %d values but the curves have length %d.', ...
        tern(opt.Transpose, 'kperp', 'kpar'), numel(kx), nrow);
end

% ------------------------------------------------------------- confidence
% Same orientation as f, so local_field applies the same transpose.
found = local_field(Output, 'flagSolutionFound', [nrow ncol], 0, opt);
suspect = local_field(Output, 'flagSuspectBranch', [nrow ncol], 0, opt);
notCert = local_field(Output, 'flagNotCertified', [nrow ncol], 0, opt);
unknown = local_field(Output, 'flagCertUnknown', [nrow ncol], 0, opt);
pred = local_field(Output, 'predicted', [nrow ncol], NaN, opt);

if strcmpi(opt.Mode, 'plain')
    suspect(:) = 0;
    notCert(:) = 0;
    unknown(:) = 0;
end

% ------------------------------------------------------------------- plot
ax = gca;
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

cGood = [0 0.35 0.75];
cSuspect = [0.85 0.1 0.1];
cNotCert = [0.95 0.55 0.1];
cUnknown = [0.5 0.5 0.5];
cPred = [0.72 0.72 0.72];

objs = struct();
objs.lines = gobjects(0);
objs.suspect = gobjects(0);
objs.notCertified = gobjects(0);
objs.certUnknown = gobjects(0);
objs.predicted = gobjects(0);

for c = 1:ncol
    fc = f(:, c);
    ok = found(:, c) ~= 0 & isfinite(fc);
    if ~any(ok)
        continue
    end
    % Insert NaN where a point is missing so no segment crosses a hole.
    y = real(fc);
    y(~ok) = NaN;
    objs.lines(end+1) = plot(ax, kx, y, '-', ...
        'Color', cGood, 'LineWidth', opt.LineWidth, ...
        'HandleVisibility', 'off'); %#ok<AGROW>

    % Filled dots only where there are few enough to read; on a dense line
    % the curve alone is clearer.
    good = ok & suspect(:, c) == 0 & notCert(:, c) == 0 & unknown(:, c) == 0;
    rows = find(good);
    if nnz(good) <= 150
        plot(ax, kx(rows), y(rows), '.', ...
            'Color', cGood, 'MarkerSize', opt.MarkerSize * 0.8, ...
            'HandleVisibility', 'off');
    end

    if opt.ShowPredicted
        py = real(pred(:, c));
        have = isfinite(py) & isfinite(y);
        if any(have)
            pyv = nan(nrow, 1);
            pyv(have) = py(have);
            objs.predicted(end+1) = plot(ax, kx, pyv, ...
                '--', 'Color', cPred, 'LineWidth', 0.9, ...
                'HandleVisibility', 'off'); %#ok<AGROW>
        end
    end

    % Marks are drawn with no connecting line: they sit on top of the curve at
    % the individual points they refer to, and a line through them would
    % re-bridge exactly the gaps the curve deliberately leaves open.
    rows = find(ok & suspect(:, c) ~= 0);
    if ~isempty(rows)
        objs.suspect(end+1) = plot(ax, kx(rows), y(rows), 'x', ...
            'Color', cSuspect, 'MarkerSize', opt.MarkerSize * 1.6, ...
            'LineWidth', 2.0, 'LineStyle', 'none'); %#ok<AGROW>
    end
    rows = find(ok & notCert(:, c) ~= 0);
    if ~isempty(rows)
        objs.notCertified(end+1) = plot(ax, kx(rows), y(rows), 'o', ...
            'Color', cNotCert, 'MarkerSize', opt.MarkerSize * 1.5, ...
            'LineWidth', 1.4, 'LineStyle', 'none'); %#ok<AGROW>
    end
    rows = find(ok & unknown(:, c) ~= 0);
    if ~isempty(rows)
        objs.certUnknown(end+1) = plot(ax, kx(rows), y(rows), 's', ...
            'Color', cUnknown, 'MarkerSize', opt.MarkerSize * 1.3, ...
            'LineWidth', 1.2, 'LineStyle', 'none'); %#ok<AGROW>
    end
end

if strcmpi(opt.YScale, 'log')
    set(ax, 'YScale', 'log');
end
% Plain-text axis labels: \parallel is a LaTeX macro and the default TeX
% interpreter rejects it, which turns into a warning on every render.
if opt.Transpose
    xlabel('k_{perp}', 'Interpreter', 'tex');
else
    xlabel('k_{par}', 'Interpreter', 'tex');
end
ylabel(sprintf('%s (real part, \\Omega_e units)', opt.Frequency), ...
    'Interpreter', 'tex');
title(sprintf('%d/%d points found, %d suspect, %d uncertified, %d unknown', ...
    nnz(found), nrow * ncol, nnz(suspect & found), ...
    nnz(notCert & found), nnz(unknown & found)));

local_legend(objs);

if nargout > 0
    varargout{1} = objs;
end
end

% ------------------------------------------------------------------ helpers

function local_legend(objs)
%LEGEND  List only the marks actually present, so the box stays small.
entries = {};
handles = gobjects(0);
if ~isempty(objs.suspect)
    handles(end+1) = objs.suspect(1); %#ok<AGROW>
    entries{end+1} = 'suspect branch (continuity rejected)'; %#ok<AGROW>
end
if ~isempty(objs.notCertified)
    handles(end+1) = objs.notCertified(1); %#ok<AGROW>
    entries{end+1} = 'not certified by argument principle'; %#ok<AGROW>
end
if ~isempty(objs.certUnknown)
    handles(end+1) = objs.certUnknown(1); %#ok<AGROW>
    entries{end+1} = 'certification refused (expansion)'; %#ok<AGROW>
end
if ~isempty(handles)
    legend(handles, entries, 'Location', 'best');
end
end

function v = local_field(Output, name, shape, fillWith, opt)
%FIELD  Output.(name), oriented like f, or fillWith where the field is absent.
%   Version 1 of the solver produced no certification flags at all, so an
%   absent field is normal and must read as 'no marks' rather than as an error.
if nargin < 4
    fillWith = 0;
end
if isfield(Output, name) && ~isempty(Output.(name))
    raw = Output.(name);
    % Exactly the transform f had, so the flags line up with the curves.
    if opt.Transpose
        raw = raw.';
    end
    if isequal(size(raw), shape)
        v = raw;
        return
    end
end
v = fillWith * ones(shape);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
