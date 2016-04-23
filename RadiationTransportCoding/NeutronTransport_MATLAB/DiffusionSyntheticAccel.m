function f = DiffusionSyntheticAccel(Nx,dx,sigt,sigs0,siga,phi_prev,phi_new)

D = 1/(3*sigt);

M = (-D/dx^2).*full(gallery('tridiag',Nx,1,-2,1)) + siga.*eye(Nx,Nx);
% M = full(gallery('tridiag',Nx+2,1,-2,1));
% M(1,1) = 0;
% M(end) = 0;
% 
% M = (-D/dx^2).*M;
% 
% for i = 2:Nx+1
%     
%     M(i,i) = M(i,i) + siga;
%     
% end

% phi_n(2:length(phi_new)+1) = phi_new;
% phi_n(1) = phi_new(1);
% phi_n(end+1) = phi_new(end);
% 
% phi_o(2:length(phi_new)+1) = phi_prev;
% phi_o(1) = phi_prev(1);
% phi_o(end+1) = phi_prev(end);

%S = sigs0.*(phi_n - phi_o);
S = sigs0.*(phi_new - phi_prev);

cond(M)

f = M\S';

%f = f(2:end-1);

return