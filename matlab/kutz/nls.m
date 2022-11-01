% Nonlinear Schrödinger Equation a 3. Kutz videóból

clear all
close all

L=30;
n=512;
x2=linspace(-L/2,L/2,n+1);
x=x2(1:n);

% hullámszámok a Fourier-transzformáltakhoz
k=(2*pi/L)*[0:n/2-1 -n/2:-1]';
t=linspace(0,2*pi*3,41*3);


% kezdeti értékek hiperbolikus szekáns függvény...
%u=sech(x); ut=fft(u);
u=sech(x); ut=fft(u);

% időléptetés Fourier-tartományban
[t,utsol] = ode45('nls_rhs',t,ut,[],k);

for j=1:length(t)
    usol(j,:)=ifft(utsol(j,:));
end

%figure(1), surfl(x,t,abs(usol))
X=usol';

[U,s,V]=svd(X);
%figure(2), semilogy(diag(s)/sum(diag(s)),'ko','linewidth',2)

figure(3), plot(real(U(:,1)),'linewidth',2)
figure(4), plot(real(V(:,1)),'linewidth',2)

% U oszlopai az egyes módusok
% V oszlopai (V^* sorai) az egyes módusok együtthatói az időlépések szerint sorban
