%Neutron Diffusion Equation K-eigenvalue Code

clc, clear

%Material Properties
D = 2;
nusigf = 3*pi;
siga = pi;

%Critical size of slab using geometric buckling equal to material bucking.
%Only for a slab. Easily extended to other coordinate systems.
Lx = pi*((nusigf-siga)/D)^(-0.5);

%Discretization of geometry
N = 50;
h = Lx/(N-1);

%Generation of Leakage and Absorption Matrices
L = (-D/h^2).*full(gallery('tridiag',N,1,-2,1));
A = siga*eye(N,N);

M = L+A;

%Boundary Conditions ( phi(0) = phi(Lx) = 0 )
M(1,1) = 1; M(1,2) = 0;
M(end,end) = 1; M(end,end-1) = 0;

phi0 = ones(N,1);
phi0(1) = 0;
phi0(N) = 0;

%K-effective guess
k = 1.00;

tol = 1e-10;

for i = 1:100
    
    kold = k;
    
    psi = M\(nusigf.*phi0);
    k = sum(nusigf.*psi)/sum(nusigf*phi0);
    phi0 = (1/k).*psi;
    phi0(1) = 0;
    phi0(N) = 0;
    
    residual = norm(abs(k-kold));
    
    if residual <= tol
        
        break
        
    end
    
end

fprintf('Calculated k-effective value: %f \n',k)
fprintf('K-effective converged with %i iterations \n',i)