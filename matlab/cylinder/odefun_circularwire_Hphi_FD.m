function dHdt = odefun_circularwire_Hphi_FD(t, H, a, M, t_pulse)

mu0 = pi*4e-7;
sigma = 35e6;

Hsurf = current(t, t_pulse)/(2*pi*a);

dHdt = 1/(mu0*sigma)*M*[0; H; Hsurf];

end