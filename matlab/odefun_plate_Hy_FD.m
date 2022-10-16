function dHdt = odefun_plate_Hy_FD(t, H, a, M, t_pulse)

mu0 = pi*4e-7;
sigma = 35e6;

Hsurf = current(t, t_pulse)/a;

dHdt = -1/(mu0*sigma)*M*[0; H; Hsurf];

end
