clearvars
close all
clc

% An infinite conductive plate is located at 0 =< z =< a
% Current flows in the x direction. The current per unit length varies
% according to a certain function of t. Herein a half-sine pulse is
% considered.
% The diffusion equation for Hy is solved by the method of lines. The
% spatial discretization is done by finite differences, whereas the time
% evolution is calculated by the fwd Euler scheme.

n = 50;             % discretization of z
a = 0.005;          % [m] thickness of the plate
t_pulse = 1e-3;     % [s] pulse duration
Tmax = 2*t_pulse;   % time window of simulation
dt   = 2e-7;        % [s] time step

z  = linspace(0, a, n); 
dz = z(2)-z(1);

sigma = 35e6;       % [S/m] conductivity
mur   = 1;

mu0   = pi*4e-7;    % [Vs/Am]
mu    = mu0*mur;
alpha = 1/(mu*sigma);
F = alpha*dt/dz^2;
disp("F="), disp(F) % convergence factor (must be <0.5)


M = matrix_for_rotrot_descartes(z);

H_init = zeros(n-2,1); % initial condition

fun = @(t,H) odefun_plate_Hy_FD(t, H, a, M, t_pulse, mu, sigma); 

nStep = ceil(Tmax/dt); % no. of steps

H_all = zeros(n, nStep);
t_all = zeros(1, nStep);

H = H_init;
t = 0;
for i = 1:nStep
    Hsurf = current(t, t_pulse)/a; % magnetic field on the surface (z=a) from Amper's law
    H_all(:, i) = [0; H; Hsurf]; % vector of H values along z
    t_all(i) = t;
   
    dHdt = fun(t,H); % dH/dt derivative
     
    H = H + dt*dHdt;
    
    t = t + dt;
end

% figure(1)
% plot(t_all, H_all(end,:), t_all, H_all(round(2/3*n),:), t_all, H_all(round(1/3*n),:))
% xlabel('t (s)')
% ylabel('H_y (A/m)')
% legend('z=a', 'z=(2/3)a', 'z=(1/3)a')

%%
% reduced order approximation
nSampled = 20;
X = H_all(2:end-1,1:nSampled);
redOrder = 3; % number of base vectors for reduced model

[U,S,V]=svd(X);
figure(2)
semilogy(diag(S)/sum(diag(S)), 'kx')
title('Singular values of snapshot matrix')
xlabel('i (no. of singular value)')
ylabel('\sigma /\Sigma \sigma_i')
S
Uhat = U(:,1:redOrder); % new base, columns are base vectors

%% matrix for the construction on M_red
%T=[zeros(1,redOrder); Uhat; zeros(1,redOrder)];
%T(redOrder,redOrder+2)=1;
%
%% projecting M to the reduced base
%M_red = Uhat'*M*T;
%M_red
%fun_red = @(t,H_red) odefun_plate_Hy_FD_red(t, H_red, a, M_red, t_pulse); 
%
%nStep_red = nStep-nSampled; % no. of steps
%
%H_all_red = [H_all(:,1:nSampled) zeros(n, nStep_red)];
%t_all_red = [t_all(:,1:nSampled) zeros(1, nStep_red)];
%
%H = (H_all(2:end-1,nSampled)'*Uhat)'; % starting from the end of sampling
%t = (nSampled-1)*dt;
%H_red = Uhat'*H;
%
%for i = nSampled+1:nStep
%    
%    Hsurf = current(t, t_pulse)/a; % magnetic field on the surface (z=a) from Amper's law
%    H_all_red = [0;H_red;Hsurf];
%    t_all(i) = t;
%   
%    dHdt_red = fun_red(t,H_red); % dH/dt derivative
%     
%    H_red = H_red + dt*dHdt_red;
%    
%    H_all(:, i) = [0; H; Hsurf]; % vector of H values along z
%    t = t + dt;
%end

