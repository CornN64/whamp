function parameters = default_parameters()
%DEFAULT_PARAMETERS One point at kperp = 0, kpar = 0.0022, from f = 0.1.
parameters = struct('fstart',0.1,'kperp',[0 10 0], ...
    'kpar',[0.0022 10 0.0022],'varyKzFirst',1,'useLog',0, ...
    'maxIterations',50);
end

