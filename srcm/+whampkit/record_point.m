function values = record_point(values, ip, iz, solution)
%RECORD_POINT Write one converged solution into the output grids.
%   The Poynting flux carries the engine's 10/(4*pi/2) unit conversion.
poynting = 10 / 4 / pi / 2;
values.f(ip, iz) = solution.X;
values.Ex(ip, iz) = solution.EFL(1);
values.Ey(ip, iz) = solution.EFL(2);
values.Ez(ip, iz) = solution.EFL(3);
values.Bx(ip, iz) = solution.BFL(1);
values.By(ip, iz) = solution.BFL(2);
values.Bz(ip, iz) = solution.BFL(3);
values.Sx(ip, iz) = real(solution.EFL(2) * conj(solution.BFL(3)) ...
    - solution.EFL(3) * conj(solution.BFL(2))) * poynting;
values.Sy(ip, iz) = real(solution.EFL(3) * conj(solution.BFL(1)) ...
    - solution.EFL(1) * conj(solution.BFL(3))) * poynting;
values.Sz(ip, iz) = real(solution.EFL(1) * conj(solution.BFL(2)) ...
    - solution.EFL(2) * conj(solution.BFL(1))) * poynting;
values.EB(ip, iz) = sqrt(sum(abs(solution.EFL).^2) ...
    / sum(abs(solution.BFL).^2));
values.VGP(ip, iz) = solution.VG(1);
values.VGZ(ip, iz) = solution.VG(2);
values.SGP(ip, iz) = solution.SG(1);
values.SGZ(ip, iz) = solution.SG(2);
values.u(ip, iz) = solution.u;
end

