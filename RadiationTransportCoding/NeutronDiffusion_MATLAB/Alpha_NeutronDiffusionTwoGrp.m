%Alpha-eigenvalue Search for Two-Group Neutron Diffusion Equation
%
%Use TwoGrpNDE function to calculate k and then uses Hill's Algorithm to
%calculated alpha. Sensitive to initial guesses. Bad guesses causes
%algorithm to either fail or converge on another eigenvalue. Needs to be
%investigated. Sensitive to N and to the tolerance.

clc, clear

%Slab width, discretization, and alpha-eigenvalue guess and modifier
N = 100;
alpha = 0.00;
evm = 0.1;
Lx = 25.4272;
tol = 1e-10;

for i = 1:100
    
    k = TwoGrpNDE(Lx,N,alpha);
    
    if i == 2
        
        alpha_prev = alpha;
        alpha = alpha + evm;
        
    elseif i > 2
        
        if abs(alpha - alpha_prev) < tol
            
            break
            
        end
        
        alpha_new = alpha_prev + (1-kprev)/(k-kprev)*(alpha-alpha_prev);
        alpha_prev = alpha;
        alpha = alpha_new;
        
    else
        
        continue
        
    end
        
    kprev = k;
    
end

fprintf('Alpha-eigenvalue converged: %f in %i iterations \n',alpha,i)