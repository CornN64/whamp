function values = emptyValues(shape)
%EMPTYVALUES Allocate the result grids that whamp.m fills.
names = {'f','Ex','Ey','Ez','Bx','By','Bz','Sx','Sy','Sz', ...
    'EB','VGP','VGZ','SGP','SGZ','u'};
values = struct();
for k = 1:numel(names)
    values.(names{k}) = complex(zeros(shape));
end
end
