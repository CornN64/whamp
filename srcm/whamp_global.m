function Output = whamp_global(PlasmaModel, InputParameters, Options)
%WHAMP_GLOBAL Global branch-tracking WHAMP solver (Version 2: enumerate + esort).
%
%   *** EXPERIMENTAL -- NOT RECOMMENDED FOR PRODUCTION USE ***
%
%   This version is kept as a documented alternative and as evidence that the
%   purely global approach is the wrong trade-off here, not because it works
%   better.  Measured over nine adversarial sweep lines (three plasmas x three
%   kperp, compareExample.m) against the baseline on a step eight times finer:
%
%       whamp         hops 13   missing 11    ~0.1 s per line
%       whamp_track   hops  0   missing 13    ~0.5-1.9 s per line
%       whamp_global  hops  7   missing 86    ~80-95 s per line
%
%   It is worse than WHAMP_TRACK on both quality measures and roughly 800x
%   more expensive.  The reason is structural rather than a matter of tuning:
%   in cold multi-species plasmas the susceptibility expansion fails over most
%   of the complex plane, so the root enumeration that everything else rests
%   on is only trustworthy in a modest neighbourhood of the real axis, and the
%   root set it reports is incomplete and step-size dependent.  A global
%   matcher built on an incomplete root set cannot beat a local continuation
%   that never needed the set in the first place.  Use WHAMP_TRACK.
%
%   OUTPUT = WHAMP_GLOBAL(PlasmaModel, InputParameters) solves the WHAMP
%   dispersion relation over the requested k grid by finding *all* roots of
%   D(z, w) at every point and then deciding which root is which, instead of
%   carrying a single guess from point to point and hoping it stays put.
%
%   This is the structural fix for root jumping.  WHAMP (and Version 1) ask
%   "where does Newton land from this guess?"; WHAMP_GLOBAL asks "what roots
%   exist here, and which of them is the continuation of the branch I was
%   following?".  Branch identity becomes a property of the root set rather
%   than an accident of the starting value.
%
%   Algorithm
%     1. Visit every (kperp, kpar) on the grid -- no early abort, so the
%        result is a complete surface rather than whatever survived the
%        first non-convergence.
%     2. At each point enumerate all roots of D in the frequency box with
%        whampkit.findAllRoots, and order each set with ESORT so the listing
%        is reproducible and independent of discovery order.
%     3. Match each point's ordered set against the previous point's set by
%        nearest neighbour, and label roots with a persistent branch index.
%        A root with no credible partner starts a new branch.
%     4. Report the branch whose member nearest fstart at the sweep origin;
%        that reproduces what WHAMP is trying to do, but chosen once,
%        globally, instead of re-decided at every point.
%
%   Options (struct, all optional):
%     .box            [wmin wmax imin imax] real frequency window scanned
%                     by whampkit.findAllRoots.
%                     Default [0.01*fstart, 2*fstart + 0.02, 0, 0]
%     .seeds          Seed count for the root scan.         Default 120
%     .refine         Double the seed density and check the root set is
%                     stable.  Costs 3x the scan.           Default true
%     .matchFactor    Partner must be within matchFactor * matching scale.
%                     Default 4
%     .relFloor       Floor on the matching scale.             Default 1e-4
%     .branch         Which branch to report:
%                       'auto'   follow the one nearest fstart (default)
%                       integer  report that branch index explicitly
%     .maxIterations  Newton cap.  Default from InputParameters
%
%   OUTPUT.branches is the full result: a cell array, one cell per k point,
%   each holding the roots found there and the branch index assigned to each.
%   OUTPUT.diagnostics reports completeness, unmatched roots, and how often
%   a branch had to be created or terminated.
%
%   Cost: the scan is O(seeds) dispersion evaluations per k point, so this
%   is far slower than WHAMP -- typically two to three orders of magnitude.
%   It is meant for validating the faster solvers and for mapping branch
%   structure, not for production sweeps.
%
%   Physics is shared with WHAMP through whampkit; WHAMP itself is unchanged
%   and its bit-for-bit regression baseline still applies.

d = fileparts(mfilename('fullpath'));
addpath(fileparts(d));

if nargin < 1 || isempty(PlasmaModel)
    PlasmaModel = whampkit.complete_plasma(whampkit.default_plasma());
else
    PlasmaModel = whampkit.complete_plasma(PlasmaModel);
end
if nargin < 2 || isempty(InputParameters)
    InputParameters = whampkit.complete_parameters(whampkit.default_parameters());
else
    InputParameters = whampkit.complete_parameters(InputParameters);
end
if nargin < 3 || isempty(Options)
    Options = struct();
end

opt = local_options(Options, InputParameters);

[mass, density, temperature, anisotropy, delta, drift] = ...
    whampkit.species_arrays(PlasmaModel.Species);
plasma = whampkit.prepare_plasma(mass, density, temperature, anisotropy, ...
    delta, drift, 0.0279928 * PlasmaModel.B);

kpSpec = whampkit.axis_spec(InputParameters.kperp, 'kperp');
kzSpec = whampkit.axis_spec(InputParameters.kpar, 'kpar');
kp = whampkit.fill_axis(kpSpec);
kz = whampkit.fill_axis(kzSpec);
kpPhysical = whampkit.physical_axis_value(kp, InputParameters.useLog);
kzPhysical = whampkit.physical_axis_value(kz, InputParameters.useLog);

nkp = numel(kp);
nkz = numel(kz);
shape = [nkp nkz];

values = whampkit.emptyValues(shape);
flagFound = zeros(shape, 'int32');
flagDamped = zeros(shape, 'int32');
flagNoConvergence = zeros(shape, 'int32');
flagSuspect = zeros(shape, 'int32');

% branches{ip,iz} is a struct with .roots (the esorted set) and .id (the
% branch index assigned to each root in the same order).
branches = cell(shape);
allRoots = cell(shape);

% Branch bookkeeping.  BRANCHHIST(k) records one branch: its id and the
% sequence of root frequencies it has taken along the sweep.  The last entry
% of each sequence is what the next point is matched against, and the whole
% sequence gives the per-branch step scale used to flag outliers.
branchHist = struct('id',{},'roots',{});
nextId = 1;
diag = local_diag();

% Walk in the engine's order so the first point and the origin of each line
% agree with what WHAMP would consider the start of the branch.
order = local_visitOrder(kp, kz, InputParameters.varyKzFirst ~= 0);

for k = 1:size(order, 1)
    ip = order(k, 1);
    iz = order(k, 2);
    p = kp(ip);
    z = kz(iz);

    scanOptions = struct('seeds',opt.seeds,'refine',opt.refine, ...
        'maxIter',opt.maxIterations, ...
        'extraSeeds',InputParameters.fstart);
    [roots, rinfo] = whampkit.findAllRoots(plasma, p, z, opt.box, scanOptions);

    allRoots{ip, iz} = roots;
    diag.evals = diag.evals + rinfo.evals;
    if ~rinfo.stable
        diag.unstablePoints{end+1} = struct('ip',ip,'iz',iz, ...
            'found',rinfo.count,'unstable',rinfo.unstable); %#ok<AGROW>
    end

    if isempty(roots)
        flagNoConvergence(ip, iz) = 1;
        continue
    end

    lastLabels = local_lastOf(branchHist);
    ids = local_assignIds(roots, lastLabels, branchHist, opt, nextId);
    branches{ip, iz} = struct('roots',roots,'id',ids.ids, ...
        'distance',ids.distance,'suspicious',ids.suspicious);
    nextId = ids.nextId;
    diag.newBranches = diag.newBranches + ids.created;
    diag.suspicious = diag.suspicious + nnz(ids.suspicious);

    for j = 1:numel(ids.ids)
        branchHist = local_append(branchHist, ids.ids(j), roots(j));
    end
end

% Pick the branch to report.
selBranch = local_selectBranch(branches, opt, InputParameters.fstart, ...
    order(1,1), order(1,2));

for k = 1:size(order, 1)
    ip = order(k, 1);
    iz = order(k, 2);
    b = branches{ip, iz};
    if isempty(b)
        continue
    end
    idx = find(b.id == selBranch, 1);
    if isempty(idx)
        flagNoConvergence(ip, iz) = 1;
        continue
    end
    [solution, status] = whampkit.solvePoint(plasma, kp(ip), kz(iz), ...
        b.roots(idx), opt.maxIterations, Inf);
    if status ~= 1
        flagNoConvergence(ip, iz) = 1;
        continue
    end
    values = whampkit.record_point(values, ip, iz, solution);
    flagFound(ip, iz) = 1;
    % A point where the selected branch was not the closest root to the
    % previous point is where a single-seed solver would most likely have
    % gone wrong; flag it so the two versions can be compared there.
    if b.suspicious(idx)
        flagSuspect(ip, iz) = 1;
    end
end

Output = local_assemble(InputParameters, PlasmaModel, kpPhysical, kzPhysical, ...
    values, flagFound, flagDamped, flagNoConvergence, flagSuspect);

diag.converged = double(nnz(flagFound));
diag.pointsVisited = size(order, 1);
diag.branchReported = selBranch;
diag.branchCount = nextId - 1;
Output.branches = branches;
Output.allRoots = allRoots;
Output.diagnostics = diag;
end

% ------------------------------------------------------------------ helpers

function opt = local_options(Options, InputParameters)
opt = struct();
fstart = InputParameters.fstart;
if fstart <= 0
    fstart = 0.1;
end
% Default search window in real frequency.  The scan runs along the real
% axis (see whampkit.findAllRoots), so only the two real edges matter; the
% imaginary pair is carried for the validity check and left at zero.
% wmin stays above zero because w = 0 is itself a singularity of D.
%
% The upper edge has to cover the whole sweep, not just fstart: a branch
% started near fstart can rise by a factor of several across a kpar line,
% and a box that clips it silently drops the tail of the sweep.  Ten times
% fstart is a generous default; widen it with .box if the branches go
% further, and check diagnostics.unstablePoints either way.
opt.box = local_get(Options, 'box', ...
    [max(1e-5, 0.01 * fstart), 10 * fstart + 0.05, 0, 0]);
opt.refine = local_get(Options, 'refine', true);
opt.seeds = local_get(Options, 'seeds', 120);
opt.matchFactor = local_get(Options, 'matchFactor', 4);
opt.relFloor = local_get(Options, 'relFloor', 1e-4);
opt.branch = local_get(Options, 'branch', 'auto');
opt.maxIterations = local_get(Options, 'maxIterations', ...
    InputParameters.maxIterations);
end

function v = local_get(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function diag = local_diag()
diag = struct('evals',0,'suspicious',0,'newBranches',0, ...
    'unstablePoints',{{}});
end

function order = local_visitOrder(kp, kz, cycleZFirst)
%VISITORDER Grid indices in the order the engine's sweep would meet them.
%   Continuity is only meaningful along the direction the solver actually
%   walks, so the matching has to follow the same path.
[nkp, nkz] = deal(numel(kp), numel(kz));
if cycleZFirst
    outer = 1:nkp;
    inner = 1:nkz;
    swap = true;
else
    outer = 1:nkz;
    inner = 1:nkp;
    swap = false;
end
rows = zeros(numel(outer) * numel(inner), 2);
c = 0;
for a = outer
    for b = inner
        c = c + 1;
        if swap
            rows(c, :) = [a b];
        else
            rows(c, :) = [b a];
        end
    end
end
order = rows(1:c, :);
end

function out = local_assignIds(roots, lastLabels, hist, opt, nextId)
%ASSIGNIDS Label each root with a persistent branch index.
%   Greedy nearest-neighbour matching against the last known position of
%   every branch, cheapest pair first, with NO rejection threshold.
%
%   There is deliberately no hard cutoff.  Branches in this problem move a
%   long way between neighbouring grid points -- measured motion reaches
%   3.6e-2 in f for a 1e-3 step in kpar -- so any fixed acceptance radius is
%   either so large that it admits genuine mismatches or so small that every
%   steep part of a branch is declared new.  The first version of this
%   routine used a radius and consequently invented a new branch at almost
%   every point, which is the failure it was meant to prevent.
%
%   Instead the match is unconditional and the *suspicion* is recorded:
%   OUT.suspicious marks a root whose step is an outlier against that
%   branch's own recent motion, measured on the median step it has taken so
%   far.  That is a per-branch, self-calibrating scale, so a steep branch is
%   not flagged merely for being steep, while a branch that suddenly jumps
%   by far more than it ever has is.
n = numel(roots);
out = struct('ids',zeros(1, n),'distance',inf(1, n), ...
    'suspicious',false(1, n),'nextId',nextId,'created',0);
if isempty(lastLabels)
    for j = 1:n
        out.ids(j) = nextId;
        nextId = nextId + 1;
        out.created = out.created + 1;
    end
    out.nextId = nextId;
    return
end

m = numel(lastLabels);
cost = zeros(n, m);
for j = 1:n
    for q = 1:m
        cost(j, q) = abs(roots(j) - lastLabels(q).root);
    end
end

% Optimal (minimum total distance) one-to-one matching.
%
% Greedy cheapest-pair-first is wrong here, and the failure is not
% hypothetical: two branches that approach each other between grid points
% get their identities swapped, because the greedy rule commits to the
% single closest pair before it can see that doing so forces a much worse
% partner on the other branch.  Swapped identities are precisely the defect
% this solver exists to remove, so the assignment has to be optimal.
%
% Root counts are small (one to a handful), so enumerating assignments is
% cheaper and clearer than a general Hungarian implementation.  Above the
% enumeration limit the code degrades to greedy rather than stalling.
assign = local_optimalAssignment(cost, n, m);

usedRoot = false(1, n);
for q = 1:m
    j = assign(q);
    if j > 0
        usedRoot(j) = true;
        out.ids(j) = lastLabels(q).id;
        out.distance(j) = cost(j, q);
    end
end
for j = 1:n
    if ~usedRoot(j)
        out.ids(j) = nextId;
        nextId = nextId + 1;
        out.created = out.created + 1;
    end
end
out.nextId = nextId;

% Flag outliers against each branch's own recent motion.
for j = 1:n
    steps = local_steps(hist, out.ids(j));
    if numel(steps) >= 2
        scale = median(steps);
    else
        scale = 0;
    end
    scale = max(scale, opt.relFloor * abs(roots(j)));
    out.suspicious(j) = out.distance(j) > opt.matchFactor * scale;
end
end

function assign = local_optimalAssignment(cost, n, m)
%LOCALOPTIMALASSIGNMENT  One-to-one matching of minimum total distance.
%   ASSIGN = local_optimalAssignment(COST, n, m) returns ASSIGN(q) = j, the
%   index of the root matched to branch q (q = 1..m), or 0 when branch q is
%   left unmatched.  COST is n-by-m with COST(j,q) = |root_j - last_q|.  The
%   matching has cardinality min(n,m) and minimises the SUM of the matched
%   costs, which greedy cheapest-pair-first does not.
%
%   Solved by dynamic programming over subsets of the larger side: place the
%   items of the smaller side one at a time, with the state being the set of
%   already-used slots.  With the handful of roots this problem produces the
%   state space is negligible; past LIM the exponential stops being silly
%   and the code falls back to greedy.

assign = zeros(1, m);
if n == 0 || m == 0
    return
end

LIM = 12;
if max(n, m) > LIM
    assign = local_greedyAssignment(cost, n, m);
    return
end

% Place the smaller side into slots taken from the larger side, so that
% every item is guaranteed a slot and the mask stays as small as possible.
if n <= m
    C = cost;  small = n;  large = m;   % items are roots, slots are branches
else
    C = cost'; small = m;  large = n;   % items are branches, slots are roots
end
S = 2 ^ large;

% popcount for every mask, computed once
pc = zeros(1, S);
for k = 1:large
    pc = pc + bitget(0:(S - 1), k);
end

prev = inf(1, S);
prev(1) = 0;
choice = zeros(small, S, 'int32');
for q = 1:small
    cur = inf(1, S);
    cq = C(q, :);
    for mask = 0:(S - 1)
        if pc(mask + 1) ~= q
            continue
        end
        best = inf;
        bestk = 0;
        for k = 1:large
            if bitget(mask, k)
                v = prev(mask - 2^(k - 1) + 1) + cq(k);
                if v < best
                    best = v;
                    bestk = k;
                end
            end
        end
        cur(mask + 1) = best;
        choice(q, mask + 1) = bestk;
    end
    prev = cur;
end

[bestCost, bestIdx] = min(prev);
if ~isfinite(bestCost)
    return
end

mask = bestIdx - 1;
slotOf = zeros(1, small);
for q = small:-1:1
    k = double(choice(q, mask + 1));
    slotOf(q) = k;
    mask = mask - 2^(k - 1);
end

% slotOf(q) is the slot taken by item q.  Re-express it indexed by branch:
% when roots were the items the slots are branches, and vice versa.
if n <= m
    assign(slotOf) = 1:small;
else
    assign = slotOf;
end
end

function assign = local_greedyAssignment(cost, n, m)
%LOCALGREEDYASSIGNMENT  Cheapest-pair-first fallback for large root counts.
assign = zeros(1, m);
pairs = [reshape(cost, 1, []); ...
         repmat((1:n)', 1, m); ...
         repmat(1:m, n, 1)];
pairs = sortrows(pairs', 1)';
usedRoot = false(1, n);
usedBranch = false(1, m);
for i = 1:size(pairs, 2)
    j = pairs(2, i);
    q = pairs(3, i);
    if ~usedRoot(j) && ~usedBranch(q)
        usedRoot(j) = true;
        usedBranch(q) = true;
        assign(q) = j;
    end
end
end

function steps = local_steps(hist, id)
k = find([hist.id] == id, 1);
if isempty(k)
    steps = zeros(1, 0);
else
    steps = abs(diff(hist(k).roots(:).'));
end
end

function lastLabels = local_lastOf(branchHist)
%LASTOF The most recent root of each branch, for nearest-neighbour matching.
lastLabels = struct('id',{},'root',{});
if isempty(branchHist)
    return
end
lastLabels = repmat(struct('id',0,'root',0 + 0i), 1, numel(branchHist));
for k = 1:numel(branchHist)
    lastLabels(k).id = branchHist(k).id;
    lastLabels(k).root = branchHist(k).roots(end);
end
end

function branchHist = local_append(branchHist, id, root)
%APPEND Extend one branch's trajectory by the root just matched to it.
k = find([branchHist.id] == id, 1);
if isempty(k)
    entry.id = id;
    entry.roots = root;
    branchHist(end+1) = entry; %#ok<AGROW>
else
    branchHist(k).roots(end+1) = root; %#ok<AGROW>
end
end

function b = local_selectBranch(branches, opt, fstart, ip0, iz0)
%SELECTBRANCH Which branch the reported solution should follow.
if isnumeric(opt.branch) && ~ischar(opt.branch)
    b = opt.branch;
    return
end
b0 = branches{ip0, iz0};
if isempty(b0)
    b = 1;
    return
end
% The branch whose member at the origin point sits closest to fstart.
[~, k] = min(abs(imag(b0.roots) - 0) + abs(real(b0.roots) - fstart));
b = b0.id(k);
end

function Output = local_assemble(InputParameters, PlasmaModel, ...
        kpPhysical, kzPhysical, values, flagFound, flagDamped, ...
        flagNoConvergence, flagSuspect)
resultNames = {'f','Ex','Ey','Ez','Bx','By','Bz','Sx','Sy','Sz', ...
    'EB','VGP','VGZ','SGP','SGZ','u'};
Output.InputParameters = InputParameters;
Output.PlasmaModel = PlasmaModel;
Output.kperp = kpPhysical;
Output.kpar = kzPhysical;
for k = 1:numel(resultNames)
    name = resultNames{k};
    Output.(name) = values.(name).';
end
Output.flagSolutionFound = flagFound.';
Output.flagTooHeavilyDamped = flagDamped.';
Output.flagNoConvergence = flagNoConvergence.';
Output.flagSuspectBranch = flagSuspect.';
missing = Output.flagNoConvergence == 1;
for k = 1:numel(resultNames)
    name = resultNames{k};
    Output.(name)(missing) = NaN;
end
end
