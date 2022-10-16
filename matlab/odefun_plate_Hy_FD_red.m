function dHdt_red = odefun_plate_Hy_FD_red(t, H_red, a, M_red, t_pulse)

mu0 = pi*4e-7;
sigma = 35e6;

Hsurf = current(t, t_pulse)/a;

dHdt_red = -1/(mu0*sigma)*M_red*[0; H_red; Hsurf];

end
