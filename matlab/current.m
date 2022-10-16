function curr = current(t,t_pulse)

curr = (1-heaviside(t-t_pulse)).*sin(pi*t/t_pulse);