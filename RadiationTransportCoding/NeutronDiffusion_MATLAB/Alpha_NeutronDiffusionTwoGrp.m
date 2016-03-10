%Alpha-eigenvalue Search for Two-Group Neutron Diffusion Equation
%
%Use TwoGrpNDE function to calculate k and then uses Hill's Algorithm to
%calculated alpha. Sensitive to initial guesses. Bad guesses causes
%algorithm to either fail or converge on another eigenvalue. Needs to be
%investigated. Sensitive to N and to the tolerance.

clc, clear

%Slab width, discretization, and alpha-eigenvalue guess and modifier
N = 100;
%alpha = 0.00;
evm = 0.1;
Lx = 25.4272;
tol = 1e-10;

alpha = -0.5:0.05:1.0;
alpha_guess = alpha;

for iter = 1:length(alpha)

for i = 1:100
    
    k = TwoGrpNDE(Lx,N,alpha(iter));
    
    if i == 2
        
        alpha_prev = alpha(iter);
        alpha(iter) = alpha(iter) + evm;
        
    elseif i > 2
        
        if abs(alpha(iter) - alpha_prev) < tol
            
            fprintf('Converged in %i iterations\n',i)
            break
            
        end
        
        alpha_new = alpha_prev + (1-kprev)/(k-kprev)*(alpha(iter)-alpha_prev);
        alpha_prev = alpha(iter);
        alpha(iter) = alpha_new;
        
    else
        
        continue
        
    end
        
    kprev = k;
    
end

alpha_calc(iter) = alpha(iter);

end




%fprintf('Alpha-eigenvalue converged: %f in %i iterations \n',alpha,i)