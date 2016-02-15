%% Neutron Diffusion Equation K-eigenvalue Criticality Calculation
% Solves neutron diffusion equation (NDE) in slab geometry.
% Finds width of critical slab using one-speed diffusion theory
% with zero flux boundary conditions on the edges. 

%clc, clear

%% Neutron Diffusion Equation in Slab with Fission Source
% The NDE in a slab is given by
%
% $$ -\frac{d}{dx}D(x)\frac{d\phi(x)}{dx} + \Sigma_a \phi(x) = \frac{1}{k}\nu
% \Sigma_f \phi(x) $$
%
% where $D(x)$ is the diffusion coefficient, $\Sigma_a$ and $\Sigma_f$ are
% the absorption and fission macroscopic cross sections, $\nu$ is the
% average number of neutrons emitted in fission, and $k$ is k-effective.

%% Material Properties
%D = 2;
%nusigf = 3*pi;
%siga = pi;
D = 0.9; nusigf = 0.070; siga = 0.09;%0.066;

%% Slab Geometry Width and Discretization
%Lx = pi*((nusigf-siga)/D)^(-0.5);
%Lx = 47.1239;
Lx = 30;

%Discretization of geometry
N = 100;
h = Lx/(N-1);

x = 0:h:Lx;

%% Generation of Leakage and Absorption Matrices
L = (-D/h^2).*full(gallery('tridiag',N,1,-2,1));
A = siga*eye(N,N);

M = L+A;

%% Boundary Conditions $(\phi(0) = \phi(L) = 0)$
M(1,1) = 1; M(1,2) = 0;
M(end,end) = 1; M(end,end-1) = 0;

phi0 = ones(N,1);
phi0(1) = 0;
phi0(N) = 0;

%% Power Iteration Scheme for k-eigenvalue and Flux
% Algorithm:
% We input an initial flux $\phi^{(0)}(x)$ and k-effective value $k_0$ and solve
% the equation:
%
% $$ M \psi^{(0)}(x) = \frac{1}{k} F \phi^{(0)}(x) $$
%
% for $\psi^{(0)}(x)$. Using this function, we calculate the next k-effective
% iterate using
%
% $$ k^{n+1} = \frac{\sum \nu \Sigma_f \psi^{(n)}(x)}{\sum \nu \Sigma_f
% \phi^{(n)}(x)} $$
%
% The new flux $\phi^{(n+1)}(x)$ is calculated
%
% $$ \phi^{(n+1)}(x) = \frac{1}{k} \psi^{(n)}(x) $$.
% 
% This is done until the two-norm difference between k-effective iterations
% is less than some tolerance.
%

beta = 0.0;

%M - 0.9.*eye(N,N);
%M(1,1) = 1; M(1,2) = 0;
%M(end,end) = 1; M(end,end-1) = 0;

% k-effective guess
k = 1.50;

%Tolerance Criterion
tol = 1e-15;

for i = 1:100
    
    kold = k;
    
    psi = (M-beta.*eye(N,N))\(nusigf.*phi0);
    k = sum(nusigf.*psi)/sum(nusigf*phi0);
    phi0 = (1/k).*psi;
    phi0(1) = 0;
    phi0(N) = 0;
    
    residual = norm(abs(k-kold));
    
    if residual <= tol
        
        break
        
    end
    
end

plot(x,phi0);
xlabel('Position (cm)');
ylabel('Flux');
title('Neutron Flux in Slab Geometry for Critical System')
grid on

fprintf('Calculated k-effective value: %f \n',k)
fprintf('K-effective converged in %i iterations \n',i)