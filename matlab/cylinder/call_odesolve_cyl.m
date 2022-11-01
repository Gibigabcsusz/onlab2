clearvars
close all
clc

% A circular conductor of infinite length carries a current i(t). Herein a 
% half-sine pulse is considered.
% The diffusion equation for H_phi is solved by the method of lines. The
% spatial discretization is done by finite differences, whereas the time
% evolution is calculated by the fwd Euler scheme.

n = 50;             % discretization of r
a = 0.005;          % [m] radius of the wire
t_pulse = 5e-4;     % [s] pulse duration, worth to try:
                    % a) 5e-4 fast variation, strong transient skin effect
                    % b) 1e-2 slow variation, almost uniform current
                    %         distribution within the wire
                    % 
Tmax = 2*t_pulse;   % time window of simulation
dt   = 2e-7;        % [s] time step
animate = 0;        % 1: animate soltion, otherwise: don't


r  = linspace(0, a, n);
dr = r(2)-r(1);

mu0   = pi*4e-7;    % [Vs/Am]
sigma = 35e6;       % [S/m] conductivity
alpha = 1/(mu0*sigma);
F = alpha*dt/dr^2;
disp("F="), disp(F) % convergence factor (must be <0.5 in the Descartes-case!)

M = matrix_for_rotrot_cyl(r);

H_init = zeros(n-2,1);

fun = @(t,H) odefun_circularwire_Hphi_FD(t, H, a, M, t_pulse); 

nStep = ceil(Tmax/dt); % no. of steps

H_all = zeros(n, nStep);
t_all = zeros(1, nStep);

H = H_init;
t = 0;
for i = 1:nStep
    Hsurf = current(t, t_pulse)/(2*pi*a); % magnetic field on the surface (r=a) from Amper's law
    H_all(:, i) = [0; H; Hsurf]; % vector of H values along z
    t_all(i) = t;
   
    dHdt = fun(t,H); % dH/dt derivative
     
    H = H + dt*dHdt;
    
    t = t + dt;
end

t_range = [0, dt*nStep];
[t_ode, H_ode] = ode45(fun, t_range, H_init, "RelTol", 0.0001);


%% calculate E_z from H_phi using E = (1/sigma)*rot(H)

s = linspace(dr/2, a-dr/2, n-1); % grid shifted by dr/2

rH_all = diag(r)*H_all; % r*H

E_all = diag(1./s)*diff(rH_all)/dr/sigma; % 1/r * (d(rH)/dr) / sigma -- expanding rot in cyl coords

%% animate the solution

minH = min(min(H_all));
minE = min(min(E_all));
maxH = max(max(H_all));
maxE = max(max(E_all));

if(animate==1)
	figure(1)
	for i = 1:round(nStep/100):nStep
	    
	    subplot(211)
	    plot(r, H_all(:,i))
	    axis([0 a minH maxH])
	    xlabel('r (m)')
	    ylabel('H_{\phi} (A/m)')
	    title(i)
	
	    subplot(212)
	    plot(s, E_all(:,i))
	    axis([0 a minE maxE])
	    xlabel('r (m)')
	    ylabel('E_z (V/m)')
	
	    pause(0.1)
	
	end
end


%%
figure(2)
plot(t_all, H_all(end,:), t_all, H_all(round(2/3*n),:), t_all, H_all(round(1/3*n),:))
xlabel('t (s)')
ylabel('H_{\phi} (A/m)');
legend('r=a', 'r=(2/3)a', 'r=(1/3)a')
title('Forward Euler')

figure(3)
plot(t_all, H_all(end,:), t_ode, H_ode(:,round(2/3*n)), t_ode, H_ode(:,round(1/3*n)))
xlabel('t (s)')
ylabel('H_{\phi} (A/m)');
legend('r=a', 'r=(2/3)a', 'r=(1/3)a')
title('ode45')
