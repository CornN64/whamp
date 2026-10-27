function whamp_track_vs_c(varargin)
%WHAMP_TRACK_VS_C Compare srcm/whamp_track.m with srcc64/whamp_track.c.
%   Runs both implementations over the same plasma models and wave-number
%   grids and reports the largest relative deviation of every output array.
%   The C side is srcc64/whamp_track_dump.exe, whose cases are named here.
%
%   Usage:
%       whamp_track_vs_c()            % all cases
%       whamp_track_vs_c('line_kpar') % one case
%
%   The MATLAB interface takes density in cm^-3, temperature in eV and B in
%   nT; the C driver takes engine units directly, so B is rebuilt from XC with
%   the same 0.0279928 factor whamp_track.m uses.  MATLAB transposes every
%   grid on output (Output.f is [kpar x kperp]) and writes NaN where
%   flagNoConvergence is set; both are undone here before comparing.

root = fileparts(fileparts(mfilename('fullpath')));
exe = fullfile(root, 'srcc64', 'whamp_track_dump.exe');
if exist(exe, 'file') ~= 2
    error('whamp_track_vs_c:MissingDriver', ...
        ['Build the dump driver first: cd srcc64 && mingw32-make ' ...
         'whamp_track_dump']);
end

verbose = false;
if nargin >= 2 && ~isempty(varargin{2})
    verbose = logical(varargin{2});
end
cases = local_cases();
if nargin >= 1 && ~isempty(varargin{1})
    wanted = varargin{1};
    keep = false(1, numel(cases));
    for k = 1:numel(cases)
        keep(k) = any(strcmp(cases{k}.name, wanted));
    end
    if ~any(keep)
        error('whamp_track_vs_c:UnknownCase', 'No case named %s.', wanted);
    end
    cases = cases(keep);
end

gyroToB = @(gyro) gyro / 0.0279928;

fprintf('%-18s %6s %11s %11s  %-8s %s\n', ...
    'case', 'points', 'max rel', 'max abs', 'worst', 'verdict');
worstOverall = 0;
for k = 1:numel(cases)
    tc = cases{k};

    plasma.B = gyroToB(tc.xc);
    plasma.Species = cell(1, numel(tc.species));
    for j = 1:numel(tc.species)
        s = tc.species{j};
        plasma.Species{j} = struct('m',s.m,'n',s.n,'t',s.t, ...
            'a',s.a,'b',s.b,'d',s.d,'vd',s.vd);
    end
    ip = struct('fstart',tc.fstart, ...
        'kperp',[tc.kperp(1) tc.kperp(2) tc.kperp(3)], ...
        'kpar', [tc.kpar(1)  tc.kpar(2)  tc.kpar(3)], ...
        'varyKzFirst',tc.varyKzFirst,'useLog',tc.useLog, ...
        'maxIterations',50);

    M = whamp_track(plasma, ip);
    [status, out] = system(sprintf('"%s" %s', exe, tc.name));
    if status ~= 0
        error('whamp_track_vs_c:DriverFailed', ...
            '%s exited with %d:\n%s', tc.name, status, out);
    end
    C = local_parse(out);

    [maxRel, maxAbs, notes, perField] = local_compare(M, C, tc);
    worstOverall = max(worstOverall, maxRel);
    if maxRel < 1e-9
        verdict = 'MATCH';
    elseif maxRel < 1e-4
        verdict = 'close';
    else
        verdict = 'DIFFERS';
    end
    worstName = '';
    worstField = 0;
    for n = 1:numel(perField)
        if perField{n}.rel > worstField
            worstField = perField{n}.rel;
            worstName = perField{n}.name;
        end
    end
    fprintf('%-18s %6d %11.3e %11.3e  %-8s %s\n', tc.name, C.nk * C.nz, ...
        maxRel, maxAbs, worstName, verdict);
    for n = 1:numel(notes)
        fprintf('    %s\n', notes{n});
    end
    if verbose
        for n = 1:numel(perField)
            if isfield(perField{n}, 'ratio')
                fprintf('    %-12s bounded by inputs: gap / predicted = %.3f\n', ...
                    perField{n}.name, perField{n}.ratio);
            else
                fprintf('    %-12s rel %.3e  abs %.3e  cell rel %.3e\n', ...
                    perField{n}.name, perField{n}.rel, perField{n}.abs, ...
                    perField{n}.elemRel);
            end
        end
        local_dump(M, C);
    end
end
fprintf('worst relative deviation over all cases: %.3e\n', worstOverall);
end

% ------------------------------------------------------------------- cases

function cases = local_cases()
%CASES Mirrors the CASES table in srcc64/whamp_track_dump.c.
%   kperp/kpar are [start step stop] as the MATLAB interface wants them.
ex1 = {struct('m',16,'n',1,'t',10,'a',5,'b',0,'d',1,'vd',1), ...
       struct('m',0, 'n',1,'t',1, 'a',1,'b',0,'d',1,'vd',0)};
ex2 = {struct('m',0, 'n',1,'t',10,'a',1,'b',0,'d',1,'vd',0), ...
       struct('m',1, 'n',1,'t',1, 'a',1,'b',0,'d',1,'vd',0)};
iso = {struct('m',16,'n',1,'t',10,'a',1,'b',0,'d',1,'vd',0), ...
       struct('m',0, 'n',1,'t',1, 'a',1,'b',0,'d',1,'vd',0)};
% Whistlereanisotropy_.m: one hot anisotropic electron species at 50 nT.  The
% root on this line is strongly complex (imaginary part reaching -0.89), so it
% is the case that exercises carrying a COMPLEX seed from point to point.  A
% seed truncated to its real part still converges over most of the line and
% only wanders onto another branch near the top, which is exactly the failure
% this must catch.
whistler = {struct('m',0,'n',1,'t',100,'a',5,'b',0,'d',1,'vd',0)};
three = {struct('m',16,'n',1, 't',10,'a',5,'b',0,'d',1,'vd',1), ...
         struct('m',1, 'n',0.5,'t',1, 'a',1,'b',0,'d',1,'vd',0.5), ...
         struct('m',0, 'n',1.5,'t',1, 'a',1,'b',0,'d',1,'vd',0)};

cases = { ...
    local_case('line_kpar', 2.8, ex1, [0 10 0], [0.002 0.001 0.008], 1, 0, 0.1);
    local_case('line_kperp', 2.8, ex1, [0.005 0.001 0.010], [0.0022 10 0.0022], 1, 0, 0.1);
    local_case('line_kpar_pfirst', 2.8, ex1, [0 10 0], [0.002 0.001 0.008], 0, 0, 0.1);
    local_case('grid2d_zfirst', 2.8, ex1, [0.005 0.001 0.008], [0.002 0.001 0.006], 1, 0, 0.1);
    local_case('grid2d_pfirst', 2.8, ex1, [0.005 0.001 0.008], [0.002 0.001 0.006], 0, 0, 0.1);
    local_case('grid2d_negstep', 2.8, ex1, [0.008 -0.001 0.005], [0.006 -0.001 0.002], 1, 0, 0.1);
    local_case('grid2d_log', 2.8, ex1, [-2.301 0.1 -2.096], [-2.699 0.1 -2.522], 1, 1, 0.1);
    local_case('ex2_grid', 0.008978, ex2, [0 0.01 0.02], [0.001 0.002 0.005], 1, 0, 0.1);
    local_case('ex2_log', 0.008978, ex2, [-2.3 0.1 -1.8], [-3.0 0.1 -2.7], 1, 1, 0.1);
    local_case('isotropic', 2.8, iso, [0 10 0], [0.002 0.001 0.006], 1, 0, 0.1);
    local_case('three_species', 2.8, three, [0.005 0.001 0.008], [0.002 0.001 0.006], 1, 0, 0.1);
    local_case('start_high', 2.8, ex1, [0 10 0], [0.002 0.001 0.008], 1, 0, 0.9);
    local_case('line89_whistler', 1.39964, whistler, ...
        [0.01 10 0.01], [0.01 0.01 0.90], 1, 0, 0.006208543881328393);
    local_case('damped_boundary', 1.39964, whistler, ...
        [0.25 10 0.25], [0.01 0.01 0.90], 1, 0, 0.006208543881328393);
    };
end

function tc = local_case(name, xc, species, kperp, kpar, varyKzFirst, useLog, fstart)
tc = struct('name',name,'xc',xc,'species',{species}, ...
    'kperp',kperp,'kpar',kpar,'varyKzFirst',varyKzFirst, ...
    'useLog',useLog,'fstart',fstart);
end

% ------------------------------------------------------------------ parsing

function C = local_parse(text)
%PARSE Read the whamp_track_dump transcript into a struct of [kperp x kpar].
lines = strsplit(strtrim(text), newline);
C = struct();
C.kperp = [];
C.kpar = [];
for k = 1:numel(lines)
    parts = strsplit(strtrim(lines{k}));
    if isempty(parts)
        continue
    end
    switch parts{1}
        case 'SIZE'
            C.nk = str2double(parts{2});
            C.nz = str2double(parts{3});
            C.flagFound = zeros(C.nk, C.nz, 'int32');
            C.flagDamped = zeros(C.nk, C.nz, 'int32');
            C.flagNoConv = zeros(C.nk, C.nz, 'int32');
            C.flagSuspect = zeros(C.nk, C.nz, 'int32');
            C.certWinding = nan(C.nk, C.nz);
            C.flagNotCertified = zeros(C.nk, C.nz, 'int32');
            C.flagCertUnknown = zeros(C.nk, C.nz, 'int32');
            C.predicted = nan(C.nk, C.nz);
            C.grids = struct();
            names = {'f','Ex','Ey','Ez','Bx','By','Bz','Sx','Sy','Sz'};
            for n = 1:numel(names)
                C.grids.(names{n}) = nan(C.nk, C.nz);
            end
            for n = {'EB','VGP','VGZ','SGP','SGZ','u'}
                C.grids.(n{1}) = nan(C.nk, C.nz);
            end
        case 'KPERP'
            C.kperp(str2double(parts{2})) = str2double(parts{3});
        case 'KPAR'
            C.kpar(str2double(parts{2})) = str2double(parts{3});
        case 'CELL'
            ip = str2double(parts{2});
            iz = str2double(parts{3});
            v = str2double(parts(4:end));
            c = 1;
            gnames = {'f','Ex','Ey','Ez','Bx','By','Bz','Sx','Sy','Sz'};
            for n = 1:numel(gnames)
                C.grids.(gnames{n})(ip, iz) = complex(v(c), v(c+1));
                c = c + 2;
            end
            for n = {'EB','VGP','VGZ','SGP','SGZ','u'}
                C.grids.(n{1})(ip, iz) = v(c); c = c + 1;
            end
            C.flagFound(ip, iz)       = v(c); c = c + 1;
            C.flagDamped(ip, iz)      = v(c); c = c + 1;
            C.flagNoConv(ip, iz)      = v(c); c = c + 1;
            C.flagSuspect(ip, iz)     = v(c); c = c + 1;
            C.certWinding(ip, iz)     = v(c); c = c + 1;
            C.flagNotCertified(ip, iz)= v(c); c = c + 1;
            C.flagCertUnknown(ip, iz) = v(c); c = c + 1;
            C.predicted(ip, iz)       = complex(v(c), v(c+1));
        case 'DIAG'
            C.diag = str2double(parts(2:end));
        otherwise
            % ignore
    end
end
if ~isfield(C, 'nk')
    error('whamp_track_vs_c:BadTranscript', 'no SIZE line in driver output');
end
end

% ---------------------------------------------------------------- comparison

function [maxRel, maxAbs, notes, perField] = local_compare(M, C, tc)
%COMPARE MATLAB Output against the parsed C dump.
%   MATLAB transposes (Output.f is [kpar x kperp]) and NaNs out the cells it
%   flagged as non-convergent; the C arrays are [kperp x kpar] and raw.
notes = {};
perField = {};
maxRel = 0;
maxAbs = 0;

% Grids must agree first, or every later comparison is meaningless.
[maxRel, maxAbs, perField] = local_field(maxRel, maxAbs, perField, 'kperp', ...
    M.kperp(:).', C.kperp(:).', []);
[maxRel, maxAbs, perField] = local_field(maxRel, maxAbs, perField, 'kpar', ...
    M.kpar(:).', C.kpar(:).', []);

% Flags: compare raw, and report every mismatch by cell.
local_flag_check(M.flagSolutionFound.',   C.flagFound,       'flagSolutionFound', notes);
local_flag_check(M.flagTooHeavilyDamped.', C.flagDamped,      'flagTooHeavilyDamped', notes);
local_flag_check(M.flagNoConvergence.',    C.flagNoConv,      'flagNoConvergence', notes);
local_flag_check(M.flagSuspectBranch.',    C.flagSuspect,     'flagSuspectBranch', notes);
local_flag_check(M.flagNotCertified.',     C.flagNotCertified,'flagNotCertified', notes);
local_flag_check(M.flagCertUnknown.',      C.flagCertUnknown, 'flagCertUnknown', notes);

% Values.  A non-convergent cell is NaN in MATLAB and zero (never written) in
% C, so it is skipped; everything else must match.
%
% Each quantity is judged against its OWN scale, not cell by cell.  Several
% derived quantities (SGP, SGZ, Sy) are zero to machine precision on these
% cases, where a difference of 1e-17 is one ulp but reads as a relative error
% of one if divided by the value at that cell.  The yardstick is the largest
% magnitude of the quantity anywhere in the sweep, so the reported number says
% how far apart the two runs are as a fraction of what the array reaches.
%
% The Poynting components and the E/B ratio are built from the field
% components, so their floor is the field floor amplified by the cancellation
% in the cross-products; they are judged against the field arrays directly
% (see local_bilinear_bound) rather than against themselves.
valid = M.flagNoConvergence.' == 0;

% The six field components share one yardstick: the largest field magnitude
% anywhere in the sweep.  On the kperp = 0 line Ez and Sx vanish by symmetry,
% so each is ~1e-16 and judging one against its own scale calls one ulp a
% fraction of itself.  Against the field they belong to, they are zero.
fieldNames = {'Ex','Ey','Ez','Bx','By','Bz'};
fieldScale = 0;
for n = 1:numel(fieldNames)
    A = abs(M.(fieldNames{n}).');
    A(~valid) = NaN;
    fieldScale = max(fieldScale, max(A(:), [], 'all', 'omitnan'));
end
for n = 1:numel(fieldNames)
    name = fieldNames{n};
    [maxRel, maxAbs, perField] = local_field(maxRel, maxAbs, perField, name, ...
        (M.(name)).', C.grids.(name), valid, fieldScale);
end

names = {'f','EB','VGP','VGZ','u'};
for n = 1:numel(names)
    name = names{n};
    [maxRel, maxAbs, perField] = local_field(maxRel, maxAbs, perField, name, ...
        (M.(name)).', C.grids.(name), valid);
end

% Sx and Sy are differences of nearly cancelling products, so their own ulp
% sits far above the field components' and a per-cell ratio says nothing.
% Judge them instead by whether the observed gap exceeds what the gaps in the
% field components they are built from could produce.
[bound, notes] = local_bilinear(M.Sz.', C.grids.Sz, ...
    M.Ex.', M.By.', M.Ey.', M.Bx.', ...
    C.grids.Ex, C.grids.By, C.grids.Ey, C.grids.Bx, valid, 'Sz', notes);
perField{end+1} = struct('name', 'Sz', 'rel', 0, 'abs', 0, ...
    'elemRel', 0, 'ratio', bound);
[bound, notes] = local_bilinear(M.Sx.', C.grids.Sx, ...
    M.Ey.', M.Bz.', M.Ez.', M.By.', ...
    C.grids.Ey, C.grids.Bz, C.grids.Ez, C.grids.By, valid, 'Sx', notes);
perField{end+1} = struct('name', 'Sx', 'rel', 0, 'abs', 0, ...
    'elemRel', 0, 'ratio', bound);
[bound, notes] = local_bilinear(M.Sy.', C.grids.Sy, ...
    M.Ez.', M.Bx.', M.Ex.', M.Bz.', ...
    C.grids.Ez, C.grids.Bx, C.grids.Ex, C.grids.Bz, valid, 'Sy', notes);
perField{end+1} = struct('name', 'Sy', 'rel', 0, 'abs', 0, ...
    'elemRel', 0, 'ratio', bound);

% SGP and SGZ are the growth rate divided by the group velocity, and on these
% cases both are near zero, so the quotient is 0/0 and its own relative
% deviation is meaningless.  Bound the gap by the gaps in numerator and
% denominator instead.
[bound, notes] = local_quotient(M.SGP.', C.grids.SGP, ...
    imag(M.f).', C.grids.f, M.VGP.', C.grids.VGP, valid, 'SGP', notes);
perField{end+1} = struct('name', 'SGP', 'rel', 0, 'abs', 0, ...
    'elemRel', 0, 'ratio', bound);
[bound, notes] = local_quotient(M.SGZ.', C.grids.SGZ, ...
    imag(M.f).', C.grids.f, M.VGZ.', C.grids.VGZ, valid, 'SGZ', notes);
perField{end+1} = struct('name', 'SGZ', 'rel', 0, 'abs', 0, ...
    'elemRel', 0, 'ratio', bound);

% predicted and certWinding are NaN where the tracker had no reference / did
% not certify; NaN must match NaN.
[maxRel, maxAbs, perField] = local_field(maxRel, maxAbs, perField, 'predicted', ...
    M.predicted.', C.predicted, []);
[maxRel, maxAbs, perField] = local_field(maxRel, maxAbs, perField, 'certWinding', ...
    M.certWinding.', C.certWinding, []);
end

function [ratio, notes] = local_quotient(QA, QB, NA, NB, DA, DB, valid, name, notes)
%QUOTIENT  Is the gap in Q = N ./ D inherited from N and D?
%   Where D is near zero the quotient is 0/0 and its own relative deviation
%   is meaningless.  Bound |dQ| by (|dN| + |Q| |dD|) / |D|, the first-order
%   change from the gaps in numerator and denominator.
ratio = 0;
for k = 1:numel(QA)
    if ~isempty(valid) && ~valid(k)
        continue
    end
    a = QA(k);
    b = QB(k);
    if isnan(a) || isnan(b)
        continue
    end
    d = abs(DB(k));
    if d == 0
        continue
    end
    pred = (abs(NA(k) - NB(k)) + abs(a) * abs(DA(k) - DB(k))) / d;
    if pred > 0 && abs(b - a) / pred > ratio
        ratio = abs(b - a) / pred;
    end
end
if ratio > 2
    notes{end+1} = sprintf( ...
        '%s gap exceeds its inputs by %.2gx: not inherited from f and VG', ...
        name, ratio); %#ok<AGROW>
end
end

function [ratio, notes] = local_bilinear(SA, SB, A1, A2, A3, A4, ...
    B1, B2, B3, B4, valid, name, notes)
%BILINEAR  Is the gap in S = coef*Re(E1 conj(E2) - E3 conj(E4)) inherited?
%   S is a difference of two products that nearly cancel, so its own ulp is
%   many orders above the field components' and a per-cell ratio of the gap
%   to S is meaningless.  Instead ask whether the observed gap could have
%   come from the gaps in the inputs: bound |dS| by the first-order change
%   |dE1||E2| + |E1||dE2| + |dE3||E4| + |E3||dE4|.  A ratio at or below one
%   says the difference is inherited from the fields and not a logic
%   divergence; a ratio well above one is a real disagreement and is noted.
ratio = 0;
n = numel(SA);
for k = 1:n
    if ~isempty(valid) && ~valid(k)
        continue
    end
    a = SA(k);
    b = SB(k);
    if isnan(a) || isnan(b)
        continue
    end
    diff = abs(b - a);
    % First-order change from the input gaps, plus the rounding of the two
    % products themselves: they nearly cancel, so the ulp of the TERMS, not of
    % the difference, is the floor.  Without it the bound is far too tight.
    pred = abs(A1(k) - B1(k)) * abs(B2(k)) + abs(B1(k)) * abs(A2(k) - B2(k)) + ...
           abs(A3(k) - B3(k)) * abs(B4(k)) + abs(B3(k)) * abs(A4(k) - B4(k)) + ...
           2 * eps * (abs(A1(k)) * abs(B2(k)) + abs(A3(k)) * abs(B4(k)));
    if pred <= 0
        continue
    end
    if diff / pred > ratio
        ratio = diff / pred;
    end
end
if ratio > 2
    notes{end+1} = sprintf( ...
        '%s gap exceeds its inputs by %.2gx: not inherited from the fields', ...
        name, ratio); %#ok<AGROW>
end
end

function [maxRel, maxAbs, perField] = local_field(maxRel, maxAbs, perField, ...
    name, A, B, valid, scale)
%FIELD  Compare one array and record its deviation.
%   SCALE, when given, is the yardstick the deviation is expressed against;
%   by default the array's own largest magnitude.
if nargin < 8 || isempty(scale)
    scale = [];
end
[rel, abs_, elemRel] = local_scan(A, B, valid, scale);
if numel(rel) == 0
    return
end
perField{end+1} = struct('name', name, 'rel', rel, 'abs', abs_, ...
    'elemRel', elemRel); %#ok<AGROW>
maxRel = max(maxRel, rel);
maxAbs = max(maxAbs, abs_);
end

function local_dump(M, C)
%DUMP Per-cell side-by-side of the quantities most likely to diverge.
fprintf('    cell  MATLAB f                     C f                          drel        pred(M/C)                        wind M   wind C\n');
for iz = 1:size(C.grids.f, 2)
    for ip = 1:size(C.grids.f, 1)
        a = M.f(iz, ip);
        b = C.grids.f(ip, iz);
        if isnan(a) && isnan(b)
            continue
        end
        if isnan(a) || isnan(b)
            fprintf('    %d,%d  MATLAB %s  C %s  flags M %d/%d C %d/%d\n', ip, iz, ...
                mat2str(a), mat2str(b), M.flagFound(iz, ip), ...
                M.flagNoConvergence(iz, ip), C.flagFound(ip, iz), ...
                C.flagNoConv(ip, iz));
            continue
        end
        rel = abs(b - a) / max(abs(a), realmin);
        if rel < 1e-16
            continue
        end
        fprintf('    %d,%d  %.15g%+.15gi  %.15g%+.15gi  %.2e  %.12g / %.12g  %g %g\n', ...
            ip, iz, real(a), imag(a), real(b), imag(b), rel, ...
            real(M.predicted(iz, ip)), real(C.predicted(ip, iz)), ...
            M.certWinding(iz, ip), C.certWinding(ip, iz));
    end
end
end

function local_flag_check(A, B, name, notes)
A = double(A);
B = double(B);
bad = find(A ~= B);
for k = 1:numel(bad)
    notes{end+1} = sprintf('%s differs at linear cell %d: MATLAB %g, C %g', ...
        name, bad(k), A(bad(k)), B(bad(k))); %#ok<AGROW>
end
end

function [rel, absd, elemRel] = local_scan(A, B, valid, scale)
%SCAN Largest deviation over the valid cells, in three measures.
%   NaN is skipped (a cell MATLAB could not resolve has nothing to compare);
%   a size mismatch is a harness bug and raises.
%
%   absd    largest |b - a|.
%   rel     absd divided by SCALE, or by the largest |a| in the array when no
%           scale is given: the deviation as a fraction of the quantity's own
%           size.  Dividing cell by cell, as elemRel does, is meaningless for
%           the derived quantities that are zero to machine precision on these
%           cases (SGZ, Sy), where a difference of 1e-17 reads as a relative
%           error of 1e+1.
%   elemRel largest per-cell |b - a| / |a|, reported for diagnostics only.
rel = [];
absd = [];
elemRel = [];
if ~isequal(size(A), size(B))
    error('whamp_track_vs_c:SizeMismatch', 'size %s vs %s', ...
        mat2str(size(A)), mat2str(size(B)));
end
if nargin < 4 || isempty(scale)
    scale = [];
end
maxDiff = 0;
maxScale = 0;
maxElem = 0;
comparable = false;
for k = 1:numel(A)
    if ~isempty(valid) && ~valid(k)
        continue
    end
    a = A(k);
    b = B(k);
    if isnan(a) || isnan(b)
        continue
    end
    comparable = true;
    cellScale = max(abs(a), abs(b));
    maxScale = max(maxScale, cellScale);
    diff = abs(b - a);
    maxDiff = max(maxDiff, diff);
    if cellScale > 0
        maxElem = max(maxElem, diff / cellScale);
    end
end
if comparable
    absd = maxDiff;
    if isempty(scale)
        scale = maxScale;
    end
    rel = maxDiff / max(scale, realmin);
    elemRel = maxElem;
end
end
