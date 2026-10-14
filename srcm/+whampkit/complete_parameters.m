function parameters = complete_parameters(parameters)
%COMPLETE_PARAMETERS Fill in omitted parameter fields and validate them.
defaults = whampkit.default_parameters();
if ~isstruct(parameters) || ~isscalar(parameters)
    error('whamp:InvalidInputParameters', ...
        'InputParameters must be a scalar structure.');
end
fieldNames = fieldnames(defaults);
for k = 1:numel(fieldNames)
    if ~isfield(parameters, fieldNames{k})
        parameters.(fieldNames{k}) = defaults.(fieldNames{k});
    end
end

for name = {'fstart','varyKzFirst','useLog'}
    if ~whampkit.is_finite_real_scalar(parameters.(name{1}))
        error('whamp:InvalidInputParameters', ...
            'InputParameters.%s must be a finite real scalar.', name{1});
    end
end
if parameters.fstart <= 0
    error('whamp:InvalidStartFrequency', 'fstart must be positive.');
end
if ~ismember(parameters.useLog, [0 1])
    error('whamp:InvalidLogOption', 'useLog must be 0 or 1.');
end
end

