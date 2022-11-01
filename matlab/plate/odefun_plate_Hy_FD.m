function dHdt = odefun_plate_Hy_FD(t, H, a, M, t_pulse, mu, sigma)
    Hsurf = current(t, t_pulse)/a;
    dHdt = -1/(mu*sigma)*M*[0; H; Hsurf];
end
