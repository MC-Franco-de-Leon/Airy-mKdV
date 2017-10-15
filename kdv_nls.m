function xxx = kdv_nls
%KdV_nls.m:Korteveg de Vries & Nonlinear Schrodinger IVP
%Pseudo-spectral integration of
%u_t + uu_x + u_xxx = 0 and
%iu_t+ u_xx + nu |u|^2 u =0
clear;
global KDVeqn NLSeqn TheEqn
global nu disp a b m k dt;
global Cu B_hat G_u;
KDVeqn = 1; NLSeqn = 2;
TheEqn = input('Enter either 1 for kdv or 2 for nls ');
%initialization of variables for the two cases.
switch TheEqn;
case NLSeqn
N=1024; period=2*pi; nu = 2; b = 3;
Top = b*sqrt(2/nu)*1.1; Bottom = - Top;
Left = - period/2; Right = period/2;
TheEquationName = 'Non-Linear Schrodinger';
case KDVeqn
N=512; period = 20; disp = 0.05;
Top = 1.75; Bottom = -0.1;
Left = - period/2; Right = period/2;
TheEquationName = 'Korteweg de Vries';
end;

h = 2*pi/N; % spatial increment (before scaling).
x = (-pi:h:pi - h); % unscaled space lattice
k = -i*[(0:N/2) (1-N/2:-1)]; % Fourier transform of d/dx.
a=period/(2*pi); % spatial scale-factor
t=0; % initial time
dt=0.01; % time step
y=a*x; % scaled space lattice
u = initial_condition(x); % initial condition
plothandle = plot(y,a*real(u));
set(plothandle,'erasemode','background');
axis([ Left Right Bottom Top]);
title(TheEquationName);
m = mfcn;
B_hat=bhatfcn; % Linear factor of nonlinear part of RHS
C_hat=(1+m)./(1-m); % Cayley transform of linear part of RHS
while 1;
Cu = C_hat.*fft(u);
G_u = G(u); %The nonlinearity of the RHS applied to u
w = u;
for n=1:3;
w = Iterator(w); % =ifft(Cu + B_hat.*fft(G_u + G(w)));
end
u=w;
set(plothandle,'ydata',a*real(u)); drawnow;
t=t+dt;
end
%Subsidiary functions used in this program
function nlnr = G(q);
global KDVeqn NLSeqn TheEqn;
switch TheEqn
case NLSeqn
nlnr = q.*abs(q).^2;
case KDVeqn
nlnr = q.^2;
end;

function ic = initial_condition(x);
global KDVeqn NLSeqn TheEqn nu a b;
switch TheEqn
case NLSeqn
ic = sqrt(2/nu)*b*sech(b*x);
case KDVeqn
ic = exp(-1.2*a*x.^2)/a;
end;
function bhat = bhatfcn;
global KDVeqn NLSeqn TheEqn nu disp a b m k dt;
switch TheEqn
case NLSeqn
bhat = 0.5*i*nu*dt./(1-m);
case KDVeqn
bhat = 0.5*dt*k./(1 - m);
end;
function mfn = mfcn;
global KDVeqn NLSeqn TheEqn nu disp a k dt;
switch TheEqn
case NLSeqn
mfn = 0.5*i*dt*k.^2;
case KDVeqn
mfn = disp*a^(-3)*0.5*dt*k.^3;
end;
function iter = Iterator(w);
global Cu B_hat G_u TheEqn;
iter = ifft(Cu + B_hat.*fft(G_u + G(w)));