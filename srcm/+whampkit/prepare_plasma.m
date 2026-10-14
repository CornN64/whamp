function plasma = prepare_plasma(mass, density, temperature, aa, delta, drift, gyro)
%PREPARE_PLASMA Normalise the model into the engine's working quantities.
%   Transcribes plasma_setup() in whamp_engine.c.  Everything is expressed
%   relative to the first species: REN is the mass ratio to it, T the
%   temperature ratio scaled by REN, and ST = sqrt(T).  By the time this runs
%   the densities are in m^-3 and the temperatures in keV.
protonMassRatio = 1836.1;
debyeKm = 12405;
electronRestEnergy = 511;

% Electrons are entered with mass 0, which would divide by zero; the engine
% gives them the unit mass ratio here and lets the charge sign do the work.
ren = protonMassRatio .* mass;
ren(ren == 0) = 1;

last = find(density ~= 0, 1, 'last');
if isempty(last)
    error('whamp:EmptyPlasma', 'At least one species must have nonzero density.');
end
if temperature(1) <= 0
    error('whamp:InvalidTemperature', ...
        'The first species temperature must be positive.');
end

% A species whose two loss-cone parameters coincide has no loss cone at all.
for j = 1:last
    if aa(j,1) == aa(j,2)
        aa(j,2) = 0;
        delta(j) = 1;
    end
end

rn = ren(1);
t = temperature ./ temperature(1) .* (ren ./ rn);
if any(t(1:last) <= 0)
    error('whamp:InvalidTemperature', 'Active species temperatures must be positive.');
end
plasma.n = last;
plasma.mass = mass(1:last);
plasma.density = density(1:last);
plasma.ren = ren(1:last) ./ rn;
plasma.st = sqrt(t(1:last));
plasma.aa = aa(1:last,:);
plasma.delta = delta(1:last);
plasma.drift = drift(1:last);
plasma.dek = debyeKm * rn;
plasma.pfq = sum(density(1:last) ./ ren(1:last)) / debyeKm;
plasma.xa = gyro / rn;
% CV is the electron thermal speed in units of c, from the first species
% temperature TR = T(1)/RN with the special-relativistic factor applied.
plasma.cv = 1 / sqrt((temperature(1) / rn) * ...
    (2 * electronRestEnergy + temperature(1) / rn) / ...
    (electronRestEnergy + temperature(1) / rn)^2);
plasma.xc = gyro;
end

