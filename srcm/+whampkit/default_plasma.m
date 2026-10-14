function model = default_plasma()
%DEFAULT_PLASMA The built-in model: 100 nT with an O+ beam on electrons.
%   Density is in cm^-3, temperature in eV, mass in proton masses (zero
%   denotes electrons), and vd the field-aligned drift in thermal speeds.
model = struct('B', 100);
model.Species = { ...
    struct('m',16,'n',1,'t',10,'a',5,'b',0,'d',1,'vd',1), ...
    struct('m',0,'n',1,'t',100,'a',1,'b',0,'d',1,'vd',0)};
end

