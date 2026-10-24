function varargout = mexwrap(varargin)
%MEXWRAP Run WHAMP through the compiled C engine with whamp.m's struct input.
%   OUTPUT = MEXWRAP(PLASMAMODEL, INPUTPARAMETERS) takes exactly the same
%   PlasmaModel and InputParameters structures as WHAMP and returns the same
%   Output structure, but the dispersion solving is done by the compiled C
%   engine rather than by the MATLAB implementation.
%
%   The implementation is the C MEX srccopt64/mexwrap.c, which links the
%   engine in and is therefore self-contained.  This file only locates that
%   MEX and puts its directory ahead of this one on the MATLAB path, so
%   pointing MATLAB at the repository is enough:
%
%       addpath(fullfile(whampRoot, 'matlab'))
%       Output = mexwrap(PlasmaModel, InputParameters)
%
%   Build the MEX first, from the engine directory, with MATLAB's mex on PATH:
%
%       cd <repository>/srccopt64
%       make mexwrap          # or: mex -O -R2017b mexwrap.c libwhamp.a
%
%   srccopt64 is the all-double engine variant.  The same wrapper can be built
%   in ../srccopt for the variant that reproduces the Fortran reference
%   bit-for-bit; see ../REPORT.md, section 2.  If both are built, the one that appears
%   first on the MATLAB path is used.
%
%   See also WHAMP, WHAMP/RUN.

if mexwrap_is_compiled()
    % Already resolved: call straight through.
    [varargout{1:nargout}] = mexwrap(varargin{:});
    return
end

dirName = engine_directory();
addpath(dirName, '-begin');

% Calling mexwrap now resolves to the MEX, because a MEX file in a directory
% ahead of this one takes precedence.  Verify that it did: if the lookup still
% lands on a mexwrap.m the MEX cannot be loaded, and calling through would
% recurse into this function forever.
if ~mexwrap_is_compiled()
    error('mexwrap:NotLoaded', ...
        ['The compiled mexwrap in "%s" could not be loaded.\n' ...
         'Build it with "make mexwrap" (or "mex -O -R2017b mexwrap.c ' ...
         'libwhamp.a") and make sure MATLAB can load the result.'], dirName);
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

function dirName = engine_directory()
%ENGINE_DIRECTORY First directory holding a mexwrap MEX.
%   The repository root is taken from the path entry that holds this very
%   file, so the result does not depend on the current directory or on how the
%   repository was added to the path.
extension = ['.' mexext];

root = repository_root();
candidates = {fullfile(root, 'srccopt64'), fullfile(root, 'srccopt'), ...
    fullfile(root, 'srcc64'), fullfile(root, 'srcc'), pwd};

for k = 1:numel(candidates)
    if isfile(fullfile(candidates{k}, ['mexwrap' extension]))
        dirName = candidates{k};
        return
    end
end

error('mexwrap:NotBuilt', ...
    ['Could not find mexwrap%s.\n' ...
     'Build it in the engine directory with "make mexwrap" (for example in ' ...
     '%s), with MATLAB''s mex on PATH, and add that directory to the MATLAB ' ...
     'path.'], extension, fullfile(root, 'srccopt64'));
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
