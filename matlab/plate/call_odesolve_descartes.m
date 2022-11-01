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
nSampled = 800;     % number of samples from the full model
redOrder = 5;       % number of base vectors for reduced model
nPlot = 100;        % number of time sample points to plot to keep plot file size small



nStep = ceil(Tmax/dt); % no. of steps
z  = linspace(0, a, n); 
dz = z(2)-z(1);
vPlot=1:round(nStep/nPlot):nStep; % time sample points for plotting

sigma = 35e6;       % [S/m] conductivity
mur   = 1;

mu0   = pi*4e-7;    % [Vs/Am]
mu    = mu0*mur;
alpha = 1/(mu*sigma);
F = alpha*dt/dz^2;
disp("F="), disp(F) % convergence factor (must be <0.5)

M = matrix_for_rotrot_descartes(n,dz);

H_init = zeros(n-2,1); % initial condition

fun = @(t,H) odefun_plate_Hy_FD(t, H, a, M, t_pulse, mu, sigma); 


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


figure(1), hold on
plot(t_all(vPlot), H_all(end,vPlot), 'b-', 'linewidth', 1)
plot(t_all(vPlot), H_all(round(2/3*n),vPlot), 'r-', 'linewidth', 1)
plot(t_all(vPlot), H_all(round(1/3*n),vPlot), 'g-', 'linewidth', 1)
xlabel('t (s)')
ylabel('H_y (A/m)')
legend('z=a', 'z=(2/3)a', 'z=(1/3)a')

%%
% reduced order approximation
X = H_all(2:end-1,1:nSampled);

[U,S,V]=svd(X);
%figure(2)
%semilogy(diag(S)/sum(diag(S)), 'kx')
%title('Singular values of snapshot matrix')
%xlabel('i (no. of singular value)')
%ylabel('\sigma /\Sigma \sigma_i')
%S
Uhat = U(:,1:redOrder); % new base, columns are base vectors

% matrix for the construction on M_red
T=zeros(n,redOrder+2);
T(2:end-1,2:end-1)=Uhat;
T(end,end)=1;
% other construction for the same T
%T2=[zeros(n,1), [zeros(1,redOrder); Uhat; zeros(1,redOrder)], zeros(n,1)];
%T2(end,end)=1;


% projecting M to the reduced base
M_red = Uhat'*M*T;
fun_red = @(t,H_red) odefun_plate_Hy_FD_red(t, H_red, a, M_red, t_pulse, mu, sigma); 

nStep_red = nStep-nSampled; % no. of new steps



H_all_red = [T'*H_all(:,1:nSampled) zeros(redOrder+2, nStep_red)];
t_all_red = [t_all(:,1:nSampled) zeros(1, nStep_red)];

H = H_all(2:end-1,nSampled); % starting from the end of sampling
H_red = Uhat'*H;
t = nSampled*dt;

for i = nSampled+1:nStep
    
    Hsurf = current(t, t_pulse)/a; % magnetic field on the surface (z=a) from Amper's law
    H_all_red(:,i) = [0;H_red;Hsurf];
    t_all(i) = t;
   
    dHdt_red = fun_red(t,H_red); % dH/dt derivative
     
    H_red = H_red + dt*dHdt_red;
    
    t = t + dt;
end


% transforming back from the reduced base to the original base
H_aa = zeros(n,nStep);
H_aa(1,:) = H_all_red(1,:);
H_aa(end,:) = H_all_red(end,:);
H_aa(2:end-1,:) = Uhat*H_all_red(2:end-1,:);


figure(1), hold on
plot(t_all(vPlot), H_aa(end,vPlot), 'b--', 'linewidth', 2)
plot(t_all(vPlot), H_aa(round(2/3*n),vPlot), 'r--', 'linewidth', 2)
plot(t_all(vPlot), H_aa(round(1/3*n),vPlot), 'g--', 'linewidth', 2)

ylim([0,200])
xlabel('t (s)')
ylabel('H_y (A/m)')
legend('z=a', 'z=(2/3)a', 'z=(1/3)a', 'z=a, approx', 'z=(2/3)a, approx', 'z=(1/3)a, approx')

