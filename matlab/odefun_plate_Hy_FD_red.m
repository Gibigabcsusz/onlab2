function dHdt_red = odefun_plate_Hy_FD_red(t, H_red, a, M_red, t_pulse, mu, sigma)
    Hsurf = current(t, t_pulse)/a;
    dHdt_red = -1/(mu*sigma)*M_red*[0; H_red; Hsurf];
end
