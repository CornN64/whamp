function model = complete_plasma(model)
%COMPLETE_PLASMA Fill in omitted plasma fields and validate what remains.
%   The magnetic field is checked last so that a malformed species reports
%   its own error first.
defaults = whampkit.default_plasma();
if ~isstruct(model) || ~isscalar(model)
    error('whamp:InvalidPlasmaModel', 'PlasmaModel must be a scalar structure.');
end
if ~isfield(model, 'B')
    model.B = defaults.B;
end
model.Species = whampkit.normalise_species(model, defaults);

speciesDefaults = struct('m',0,'n',1,'t',1,'a',1,'b',0,'d',1,'vd',0);
fieldNames = fieldnames(speciesDefaults);
for j = 1:numel(model.Species)
    if ~isstruct(model.Species{j}) || ~isscalar(model.Species{j})
        error('whamp:InvalidSpecies', 'Each species must be a scalar structure.');
    end
    for k = 1:numel(fieldNames)
        name = fieldNames{k};
        if ~isfield(model.Species{j}, name)
            model.Species{j}.(name) = speciesDefaults.(name);
        end
        if ~whampkit.is_finite_real_scalar(model.Species{j}.(name))
            error('whamp:InvalidSpecies', ...
                'Species field %s must be a finite real scalar.', name);
        end
    end
end

if ~isnumeric(model.B) || ~isscalar(model.B) || ~isfinite(model.B) || model.B <= 0
    error('whamp:InvalidMagneticField', 'PlasmaModel.B must be a positive finite scalar.');
end
end

