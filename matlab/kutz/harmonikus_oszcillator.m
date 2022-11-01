% harmonikus oszcillátoros példa
% J. Nathan Kutz 2. POD-s videója alapján

clearvars; close all; clc

L=30;
n=512;
x2=linspace(-L/2,L/2,n+1);
x=x2(1:n);
t=0:0.2:20; % időlépések
V=(x.^2).'; % potenciál

% hullámszámok az fft-hez idomítva
k=(2*pi/L)*[0:n/2-1 -n/2:-1].';

u=exp(-0.2*x.^2);
ut=fft(u);

[t,utsol]=ode45('harm_rhs',t,ut,[],k,V);

for j=1:length(t)
    usol(j,:)=ifft(utsol(j,:));
end

figure(1)
surfl(x,t,abs(usol)), shading interp

[bal,sv,jobb]=svd(usol);

X=sv';

figure(2)
plot(abs(diag(X)/sum(diag(X)))', 'ko', 'linewidth',2)
hold on
title('Szinguláris értékek abszolút értékei')

