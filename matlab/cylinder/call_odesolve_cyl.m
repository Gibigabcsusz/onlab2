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
save = 1;           % 1: save pictures, otherwise don't
% sampled time length relative to Tmax
relative_t_sampled = 0.15;
% number of samples from the full Forward Euler model
t_sampled_end_fe = relative_t_sampled*Tmax;
% number of samples from the full ode45 model
t_sampled_end_ode = relative_t_sampled*Tmax;
redOrder = 4;       % number of base vectors for reduced model
n_plot = 100;        % number of time sample points to plot to keep plot file size small
lw = 1;
mur = 1;            % relative permeability
sigma = 35e6;       % [S/m] conductivity

r  = linspace(0, a, n);
dr = r(2)-r(1);

mu0   = pi*4e-7;    % [Vs/Am]
mu = mu0*mur;
alpha = 1/(mu*sigma);
F = alpha*dt/dr^2;
disp("F="), disp(F) % convergence factor (must be <0.5 in the Descartes-case!)

M = matrix_for_rotrot_cyl(r);

H_init = zeros(n-2,1);

fun = @(t,H,M) odefun_circularwire_Hphi_FD(t, H, a, M, t_pulse, mu, sigma); 
funM = @(t,H) fun(t,H,M); 


nStep = ceil(Tmax/dt); % no. of steps

H_all = zeros(n, nStep);
t_all = zeros(1, nStep);

H = H_init;
t = 0;

%% full solution using forward euler method
for i = 1:nStep
    Hsurf = current(t, t_pulse)/(2*pi*a); % magnetic field on the surface (r=a) from Amper's law
    H_all(:, i) = [0; H; Hsurf]; % vector of H values along z
    t_all(i) = t;
    dHdt = fun(t,H,M); % dH/dt derivative
    H = H + dt*dHdt;
    t = t + dt;
end

%% reduced order preparation for Forward Euler method
nSampled_fe = round(t_sampled_end_fe/dt);
[U_fe, S_fe, V_fe] = svd(H_all(2:end-1,1:nSampled_fe));
U_hat_fe = U_fe(:,1:redOrder); % reduced base

% matrix for the construction on M_red
C_fe = zeros(n,redOrder+2);
C_fe(2:end-1,2:end-1) = U_hat_fe;
C_fe(end,end) = 1;

% projecting M to the reduced base
M_red_fe = U_hat_fe'*M*C_fe;

%% reduced order solution using Forward Euler method
nStep_red = nStep-nSampled_fe; % no. of new steps
H_all_red_fe = [C_fe'*H_all(:,1:nSampled_fe) zeros(redOrder+2, nStep_red)];

H = H_all(2:end-1,nSampled_fe); % starting from the end of sampling
H_red = U_hat_fe'*H;
t = nSampled_fe*dt;

for i = nSampled_fe+1:nStep
    Hsurf = current(t, t_pulse)/(2*pi*a); % magnetic field on the surface (r=a) from Amper's law
    H_all_red_fe(:,i) = [0; H_red; Hsurf];
    dHdt_red = fun(t, H_red, M_red_fe); % dH/dt derivative
    H_red = H_red + dt*dHdt_red;
    t = t + dt;
end

% converting back from reduced base to original base
H_aa_fe = zeros(n,nStep);
H_aa_fe(1,:) = H_all_red_fe(1,:);
H_aa_fe(end,:) = H_all_red_fe(end,:);
H_aa_fe(2:end-1,:) = U_hat_fe*H_all_red_fe(2:end-1,:);

%% full solution using the built-in ode45 solver
t_range = [0, Tmax];
[t_ode, H_ode] = ode45(funM, t_range, H_init);
t_ode = t_ode';
H_ode = H_ode';
H_all_ode = [zeros(1, size(t_ode,2)); H_ode; current(t_ode, t_pulse)/(2*pi*a)];

%% reduced order preparation for ode45
% extracting sampled part of full soliution for reduced order solution
nSampled_ode=find(t_ode<t_sampled_end_ode,1,'last');
t_ode_sampled = t_ode(1:nSampled_ode);
H_ode_sampled = H_ode(:,1:nSampled_ode);

% reduced base
[U_ode, S_ode, V_ode] = svd(H_ode_sampled);
U_hat_ode = U_ode(:,1:redOrder);

% matrix for the construction on M_red
C_ode = zeros(n,redOrder+2);
C_ode(2:end-1,2:end-1) = U_hat_ode;
C_ode(end,end) = 1;

% projecting M to the reduced base
M_red_ode = U_hat_ode'*M*C_ode;

%% reduced order solution using the built-in ode45 solver
fun_red_ode = @(t, H_red) fun(t, H_red, M_red_ode);

t_range = [t_ode_sampled(end), Tmax];
H_init = H_ode_sampled(:,end);
H_init_red = U_hat_ode'*H_init;

[t_ode_red, H_ode_red] = ode45(fun_red_ode, t_range, H_init_red);

t_ode_red = t_ode_red(2:end)';
H_ode_red = H_ode_red(2:end,:)';

t_all_ode_red = [t_ode_sampled, t_ode_red];

H_aa_ode = zeros(n,size(t_ode_red,2)+size(t_ode_sampled,2));
H_aa_ode(end,:) = current(t_all_ode_red,t_pulse)/(2*pi*a);
H_aa_ode(2:end-1,:) = [H_ode_sampled, U_hat_ode*H_ode_red];



%% calculate E_z from H_phi using E = (1/sigma)*rot(H)
s = linspace(dr/2, a-dr/2, n-1); % grid shifted by dr/2
rH_all = diag(r)*H_all; % r*H
E_all = diag(1./s)*diff(rH_all)/dr/sigma; % 1/r * (d(rH)/dr) / sigma -- expanding rot in cyl coords

%% animate the solution
minH = min(min(H_all));
minE = min(min(E_all));
maxH = max(max(H_all));
maxE = max(max(E_all));

%% animating the solution
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

%% comparing the two Forward Euler solutions
ix = 1:round(size(t_all,2)/n_plot):size(t_all,2); % index vector for plotting
figure(2)
hold on
plot(t_all(ix), H_all(end,ix), '--', 'linewidth', lw)
plot(t_all(ix), H_all(round(0.7*n),ix), '--', 'linewidth', lw)
plot(t_all(ix), H_all(round(0.3*n),ix), '--', 'linewidth', lw)
plot(t_all(ix), H_aa_fe(round(0.7*n),ix), 'linewidth', lw)
plot(t_all(ix), H_aa_fe(round(0.3*n),ix), 'linewidth', lw)
axis([0 Tmax minH maxH])
xlabel('t [s]')
ylabel('H_{\phi} [A/m]');
xline((nSampled_fe-1)*dt,'linewidth',lw);
legend('r=a', 'r=0.7a', 'r=0.3a', 'r=0.7a (redukált)', 'r=0.3a (redukált)', 'POD bemenet vége')
title('Teljes (n='+string(n)+') és Redukált bázisú (n='+string(redOrder)+') Előrelépő Euler')
hold off

%% comparing the two ode45 solutions
ix_ode = 1:round(size(t_ode,2)/n_plot):size(t_ode,2);
ix_all_ode_red = 1:round(size(t_all_ode_red,2)/n_plot):size(t_all_ode_red,2);
figure(3)
hold on
plot(t_all(ix), H_all(end,ix), '--', 'linewidth', lw)
plot(t_ode(ix_ode), H_all_ode(round(0.7*n),ix_ode), '--', 'linewidth', lw)
plot(t_ode(ix_ode), H_all_ode(round(0.3*n),ix_ode), '--', 'linewidth', lw)
plot(t_all_ode_red(ix_all_ode_red), H_aa_ode(round(0.7*n),ix_all_ode_red), 'linewidth', lw)
plot(t_all_ode_red(ix_all_ode_red), H_aa_ode(round(0.3*n),ix_all_ode_red), 'linewidth', lw)
axis([0 Tmax minH maxH])
xlabel('t (s)')
ylabel('H_{\phi} [A/m]');
xline(t_ode_sampled(end),'linewidth',lw);
legend('r=a', 'r=0.7a', 'r=0.3a', 'r=0.7a (redukált)', 'r=0.3a (redukált)', 'POD bemenet vége')
title('Teljes (n='+string(n)+') és Redukált bázisú (n='+string(redOrder)+') ode45')
hold off

%% plottong the importance of the vectors in the reduced base
% forward Euler case
figure(4)
semilogy(diag(S_fe),'k.','linewidth', lw) %(1:10,1:10)
hold on
axis([0 49 1e-18 1e5])
xticks(0:5:50)
grid on
xline(8.5,'linewidth', lw)
xline(35.5,'linewidth', lw)
ylabel('\Sigma_{i,i}')
xlabel('i')
title('Bázisvektorok fontossága')
annotation("textbox",[0.122 0.1 0.15 0.1],"String","Fontos",'FitBoxToText','on',"HorizontalAlignment","center")
annotation("textbox",[0.395 0.1 0.15 0.1],"String","Nem fontos",'FitBoxToText','on',"HorizontalAlignment","center")
annotation("textbox",[0.723 0.1 0.15 0.1],"String","Numerikus hiba",'FitBoxToText','on',"HorizontalAlignment","center")
hold off

%% plottong some base vectors
% forward Euler case
figure(5)
hold on
plot(U_ode(:,1), 'linewidth', lw)
plot(U_ode(:,2), 'linewidth', lw)
plot(U_ode(:,3), 'linewidth', lw)
plot(U_ode(:,4), 'linewidth', lw)
plot(U_ode(:,5), 'linewidth', lw)
legend('\Psi_1', '\Psi_2', '\Psi_3', '\Psi_4', '\Psi_5', 'Location', 'northwest')
xlabel('Mintavételi pontok a sugár mentén')
ylabel('Vektor értéke az adott pontban')
title('A redukált bázis legfontosabb vektorai')
ylim([-0.5 0.5]);
hold off

figure(6)
hold on
plot(U_ode(:,11), 'linewidth', lw)
plot(U_ode(:,12), 'linewidth', lw)
plot(U_ode(:,13), 'linewidth', lw)
plot(U_ode(:,14), 'linewidth', lw)
legend('\Psi_{11}', '\Psi_{12}', '\Psi_{13}', '\Psi_{14}', 'Location', 'northwest')
xlabel('Mintavételi pontok a sugár mentén')
ylabel('Vektor értéke az adott pontban')
title('A redukált bázis kevésbé fontos vektorai')
ylim([-0.5 0.5]);
hold off

%% calculating error compared to the full model
error_fe = abs(H_all(2:end-1,nSampled_fe+1:nStep)-H_aa_fe(2:end-1,nSampled_fe+1:nStep));
rel_error_fe_all = error_fe'./H_all(2:end-1,nSampled_fe+1:nStep)';
rel_max_error_fe_all = error_fe'/max(max(H_all));

rel_error_fe = zeros(n-2,n_plot);
rel_max_error_fe = zeros(n-2,n_plot);
t_plot = t_sampled_end_fe:((Tmax-t_sampled_end_fe)/(n_plot-1)):Tmax;
for i=1:n-2
    rel_error_fe(i,:) = interp1(t_all(nSampled_fe+1:nStep),rel_error_fe_all(:,i),t_plot)';
    rel_max_error_fe(i,:) = interp1(t_all(nSampled_fe+1:nStep),rel_max_error_fe_all(:,i),t_plot)';
end

figure(7)
imagesc(rel_error_fe,"XData",t_plot,"YData",dr:dr:(a-dr))
hold on
xlabel('t [s]')
ylabel('r [m]')
colorbar;
title('Hiba az összetartozó ponthoz képest');
hold off

figure(8)
imagesc(rel_max_error_fe,"XData",t_plot,"YData",dr:dr:(a-dr))
hold on
xlabel('t [s]')
ylabel('r [m]')
title('Hiba a globális maximumhoz képest');
colorbar;

hold off

%% plotting the solution's waterfall diagram
figure(9)
imagesc(H_all,"XData",t_plot,"YData",dr:dr:(a-dr))
hold on
colorbar;
title('H_{\phi} vízesésdiagramja');
xlabel('t [s]')
ylabel('r [m]')
hold off
figure(10)

imagesc(E_all,"XData",t_plot,"YData",dr:dr:(a-dr))
hold on
colorbar;
title('E_{z} vízesésdiagramja');
xlabel('t [s]')
ylabel('r [m]')
hold off


%% Export images
if(save==1)
    figure(2) % select Euler scheme figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_td'+'.eps','ContentType','vector'); % save
    
    figure(3) % select Euler scheme figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'ode45_'+string(relative_t_sampled)+'_'+string(redOrder)+'_td'+'.eps','ContentType','vector'); % save
    
    figure(4) % select vector importance figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_sv'+'.eps','ContentType','vector'); % save
    
    figure(5) % select vector importance figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_base_1_4'+'.eps','ContentType','vector'); % save
    
    figure(6) % select vector importance figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_base_11_14'+'.eps','ContentType','vector'); % save
    
    figure(7) % select vector importance figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_rel_error'+'.eps','ContentType','vector'); % save
    
    figure(8) % select vector importance figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_rel_max_error'+'.eps','ContentType','vector'); % save
    
    figure(9) % select vector importance figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_hphi_waterfall'+'.eps','ContentType','vector'); % save
    
    figure(10) % select vector importance figure
    ax=gca; % get currently selected figure
    exportgraphics(ax,'euler_'+string(relative_t_sampled)+'_'+string(redOrder)+'_ez_waterfall'+'.eps','ContentType','vector'); % save
end


