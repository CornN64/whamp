function [mass, density, temperature, anisotropy, delta, drift] = species_arrays(species)
%SPECIES_ARRAYS Expand the species structs into the engine's component arrays.
%   Converts to engine units (density m^-3, temperature keV) and pads every
%   array to the engine's ten component slots.  As in the C driver the unused
%   slots are zero, except that delta is one.
nSpecies = numel(species);
mass = zeros(1, 10);
density = zeros(1, 10);
temperature = zeros(1, 10);
anisotropy = zeros(10, 2);
delta = ones(1, 10);
drift = zeros(1, 10);
for j = 1:nSpecies
    mass(j) = species{j}.m;
    density(j) = species{j}.n * 1e6;
    temperature(j) = species{j}.t / 1e3;
    anisotropy(j, :) = [species{j}.a, species{j}.b];
    delta(j) = species{j}.d;
    drift(j) = species{j}.vd;
end
end

