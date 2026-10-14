function species = normalise_species(model, defaults)
%NORMALISE_SPECIES Present PlasmaModel.Species as a cell array of structs.
if ~isfield(model, 'Species') || isempty(model.Species)
    species = defaults.Species;
elseif isstruct(model.Species)
    species = num2cell(model.Species);
elseif iscell(model.Species)
    species = model.Species;
else
    error('whamp:InvalidSpecies', ...
        'PlasmaModel.Species must be a structure or cell array.');
end
end

