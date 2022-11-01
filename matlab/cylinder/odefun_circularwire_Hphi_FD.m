function dHdt = odefun_circularwire_Hphi_FD(t, H, a, M, t_pulse, mu, sigma)
    Hsurf = current(t, t_pulse)/(2*pi*a);
    dHdt = 1/(mu*sigma)*M*[0; H; Hsurf];
end
