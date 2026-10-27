function varargout = mexwrap(varargin)
%MEXWRAP Run WHAMP through the compiled C engine with whamp.m's struct input.
%   OUTPUT = MEXWRAP(PLASMAMODEL, INPUTPARAMETERS) takes exactly the same
%   PlasmaModel and InputParameters structures as WHAMP and returns the same
%   Output structure, but the dispersion solving is done by the compiled C
%   engine rather than by the MATLAB implementation.
%
%   The implementation is a MEX built from mexwrap.c, which links the engine
%   in and is therefore self-contained: a MEX built in one directory uses only
%   that directory's sources and libwhamp.a, and never reaches into a sibling
%   tree.  This file only locates the MEX and puts its directory ahead of this
%   one on the MATLAB path, so pointing MATLAB at the repository is enough:
%
%       addpath(fullfile(whampRoot, 'srcm'))
%       Output = mexwrap(PlasmaModel, InputParameters)
%
%   Do not copy a built mexwrap.mex* into this directory: a MEX beside this
%   file takes precedence over it, so the variant checks below would never
%   run.  Build it in the engine directory and let this file find it.
%
%   ENGINE VARIANTS.  mexwrap.c and the engine live in two sibling
%   directories, and they differ numerically, so which one runs matters:
%
%       srcc64   all-double, tidied  (default; matches whamp_track.m)
%       srcc     bit-for-bit against the Fortran reference
%
%   Select one with the WHAMP_ENGINE environment variable; srcc64 is used
%   when it is unset.  The choice is never silently overridden: if the selected
%   variant has not been built, or if the binary found is not that variant,
%   MEXWRAP errors instead of quietly running a different engine.
%
%       setenv('WHAMP_ENGINE', 'srcc')      % or set WHAMP_ENGINE before
%                                           % starting MATLAB
%
%   ENGINE = MEXWRAP('engine') returns the name of the variant in use without
%   running anything.  Each binary stamps its own name in at compile time and
%   every result carries it in Output.engine, so a number can always be traced
%   to the tree that produced it.
%
%   Build the MEX first, from the engine directory, with MATLAB's mex on PATH:
%
%       cd <repository>/srcc64
%       make mexwrap
%
%   Switching variants within a session needs CLEAR MEXWRAP first, because
%   MATLAB keeps a MEX loaded until it is cleared.
%
%   See ../REPORT.md, section 2 for the C/Fortran comparison, and section 6
%   for what distinguishes the two trees.
%
%   See also WHAMP, WHAMP/RUN.

if ~mexwrap_is_compiled()
    % Nothing loaded yet: locate the requested variant and put it ahead of
    % this file.
    dirName = engine_directory();
    addpath(dirName, '-begin');

    % mexwrap now resolves to the MEX, because a MEX in a directory ahead of
    % this one takes precedence.  Verify that it did: if the lookup still
    % lands on a mexwrap.m the MEX cannot be loaded, and calling through
    % would recurse into this function forever.
    if ~mexwrap_is_compiled()
        error('mexwrap:NotLoaded', ...
            ['The compiled mexwrap in "%s" could not be loaded.\n' ...
             'Build it with "make mexwrap" and make sure MATLAB can load the ' ...
             'result.'], dirName);
    end
end

% The directory a binary sits in is not proof of what was linked into it, and
% whichever MEX MATLAB found first is the one that gets loaded regardless of
% what this file asked for.  Ask the loaded binary which variant it is and
% refuse to continue if that disagrees, so a copied, renamed, or un-stamped
% MEX can never pass itself off as another variant.
variant = mexwrap('engine');
if ~ischar(variant) || isempty(variant)
    % A MEX built without the Makefile's -D flag, or an older build that
    % predates the query form.  It would still run, but its provenance is
    % unknowable, which is exactly the ambiguity this check exists to remove.
    error('mexwrap:VariantUnknown', ...
        ['The mexwrap loaded from "%s" does not report its engine variant.\n' ...
         'Rebuild it with "make mexwrap" in that directory.'], ...
        fileparts(which('mexwrap')));
end
if ~strcmp(variant, requested_variant())
    error('mexwrap:VariantMismatch', ...
        ['You asked for the "%s" engine variant, but the mexwrap that was ' ...
         'loaded reports itself as "%s".\n  Loaded from: %s\n' ...
         'If you changed WHAMP_ENGINE during this session, run CLEAR MEXWRAP ' ...
         'first.  Otherwise remove or rebuild the binary that was loaded, and ' ...
         'check for a stray mexwrap.%s earlier on the MATLAB path.'], ...
        requested_variant(), variant, which('mexwrap'), mexext);
end

[varargout{1:nargout}] = mexwrap(varargin{:});
end


function compiled = mexwrap_is_compiled()
%MEXWRAP_IS_COMPILED True when the name mexwrap resolves to a MEX file.
implementation = which('mexwrap');
suffix = mexext;
compiled = ~isempty(implementation) && ...
    length(implementation) >= length(suffix) && ...
    strcmpi(implementation(end - length(suffix) + 1:end), suffix);
end

function variant = requested_variant()
%REQUESTED_VARIANT The variant the caller asked for.
%   WHAMP_ENGINE is the switch.  srcc64 is the default because it is the
%   all-double engine that whamp_track.m and the published comparison
%   were built against.
variant = getenv('WHAMP_ENGINE');
if isempty(variant)
    variant = 'srcc64';
elseif ~any(strcmp(variant, known_variants()))
    error('mexwrap:UnknownVariant', ...
        ['WHAMP_ENGINE is "%s", which is not an engine variant.\n' ...
         'Known variants: %s.'], variant, strjoin(known_variants(), ', '));
end
end

function variants = known_variants()
variants = {'srcc64', 'srcc'};
end

function dirName = engine_directory()
%ENGINE_DIRECTORY Directory holding the requested variant's mexwrap MEX.
%   The repository root is taken from the path entry that holds this very
%   file, so the result does not depend on the current directory or on how
%   the repository was added to the path.  Only the requested variant is
%   considered -- there is no fallback list, because silently running a
%   numerically different engine is worse than failing.
root = repository_root();
variant = requested_variant();
dirName = fullfile(root, variant);

if ~isfile(fullfile(dirName, ['mexwrap.' mexext]))
    built = {};
    names = known_variants();
    for k = 1:numel(names)
        if isfile(fullfile(root, names{k}, ['mexwrap.' mexext]))
            built{end+1} = names{k}; %#ok<AGROW>
        end
    end
    message = sprintf(['Could not find mexwrap.%s for the "%s" engine variant.\n' ...
        '  Expected: %s\n' ...
        '  Build it with:  cd %s  &&  make mexwrap'], ...
        mexext, variant, fullfile(dirName, ['mexwrap.' mexext]), dirName);
    if ~isempty(built)
        message = sprintf(['%s\n  Built variants: %s.  Select one with\n' ...
            '  setenv(WHAMP_ENGINE, %s%s%s)'], message, strjoin(built, ', '), ...
            char(39), built{1}, char(39));
    end
    error('mexwrap:NotBuilt', '%s', message);
end
end

function root = repository_root()
%REPOSITORY_ROOT Parent of the directory holding this mexwrap.m.
%   mfilename('fulltext') comes back empty when the function is reached
%   through run() or eval(), so fall back to searching the MATLAB path.
thisFile = which([mfilename '.m']);
if isempty(thisFile)
    entries = strsplit(path, pathsep);
    for k = 1:numel(entries)
        candidate = fullfile(entries{k}, [mfilename '.m']);
        if isfile(candidate)
            thisFile = candidate;
            break
        end
    end
end
root = fileparts(fileparts(thisFile));
end
