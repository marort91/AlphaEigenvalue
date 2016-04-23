function [S,f] = DiffusionSyntheticAccel_FiniteDiff(Nx,dx,sigt,sigs0,siga,phi_prev,phi_new)

D = 1/(3*sigt);

M = (-D/dx^2).*full(gallery('tridiag',Nx,1,-2,1)) + siga.*eye(Nx,Nx);

S = sigs0.*(phi_new-phi_prev);

f = M\S';

return