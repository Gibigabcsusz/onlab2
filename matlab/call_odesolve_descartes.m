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

mu0   = pi*4e-7;    % [Vs/Am]
sigma = 35e6;       % [S/m] conductivity
alpha = 1/(mu0*sigma);
F = alpha*dt/dz^2;
disp("F="+F) % convergence factor (must be <0.5)


M = matrix_for_rotrot_descartes(z);

H_init = zeros(n-2,1); % initial condition

fun = @(t,H) odefun_plate_Hy_FD(t, H, a, M, t_pulse); 

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

figure(1)
plot(t_all, H_all(end,:), t_all, H_all(round(2/3*n),:), t_all, H_all(round(1/3*n),:))
xlabel('t (s)')
ylabel('H_y (A/m)')
legend('z=a', 'z=(2/3)a', 'z=(1/3)a')

%%
% reduced order approximation
nSampled = 20;
X = H_all(2:end-1,1:nSampled);
redOrder = 3; % number of base vectors for reduced model

[U,S,V]=svd(X);
% figure(2)
% semilogy(diag(S)/sum(diag(S)), 'kx')
% title('Singular values of snapshot matrix')
% xlabel('i (no. of singular value)')
% ylabel('\sigma /\Sigma \sigma_i')

Uhat = U(:,1:redOrder); % new base, columns are base vectors

H = (H_all(2:end-1,nSampled)'*Uhat)'; % starting from the end of sampling
t = (nSampled-1)*dt;


% the reduced version of M is missing...
% T = [zeros(1,redOrder); Uhat; zeros(1,redOrder)];
% T(1,1)=1;
% T()
% 
% M_red = Uhat'*M*[zeros(1,redOrder); Uhat; zeros(1,redOrder)];