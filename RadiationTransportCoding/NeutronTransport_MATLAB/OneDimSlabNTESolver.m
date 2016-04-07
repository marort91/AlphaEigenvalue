%% One-dimension, One-speed Discrete Ordinates Neutron Transport Equation Solver
% Solves the 1D neutron transport code using the discrete ordinates method.
% Current work in progress.

%% Neutron Transport Equation
% The general form of the one-speed neutron transport equation is given by 
%
% $$ \mathbf{\hat{\Omega}} \cdot \nabla \psi(\vec{r},\mathbf{\hat{\Omega}})
% + \Sigma_{t}\psi(\vec{r},\mathbf{\hat{\Omega}}) = \int_{4 \pi} \Sigma_{s}
% (\vec{r},\mathbf{\hat{\Omega'}} \rightarrow \mathbf{\hat{\Omega}})
% \psi(\vec{r},\mathbf{\hat{\Omega}}) d \mathbf{\hat{\Omega'}} +
% S(\vec{r},\mathbf{\hat{\Omega}}).
% $$
%
% For the one-dimensional slab geometry, isotropic scattering case the
% equation reduces to
%
% $$ \mu \frac{\partial}{\partial x} \psi(x,\mu) + \Sigma_{t}(x)
% \psi(x,\mu) = \frac{\Sigma_{s}(x)}{2} \int_{-1}^{1} \psi(x,\mu')d\mu' +
% S(x,\mu) $$
%
% To solve the equation, we define N discrete directions and corresponding
% weighting coefficients such that we can approximate the intergral as 
%
% $$ \int_{-1}^{1} d\mu' \psi(\vec{r},\mu') = \phi(\vec{r}) \approx
% \sum_{n=1}^{N} w_{n} \psi(x,\mu_{n}). $$
%
% Using this discretization, the equation is transformed into the
% following:
%
% $$ \mu_{m} \frac{\partial}{\partial x} \psi(x,\mu_{m}) + \Sigma_{t}
% \psi(x,\mu_{m}) = \frac{\Sigma_{s}(x)}{2} \sum_{n=1}^{N} w_{n}
% \psi(x,\mu_{n}) + S(x,\mu_{m}) $$

%% Material Properties

clc, clear, clf, close all
sigt = 1;
sigs0 = 0.2;

if ( sigs0 == 0 )
    
    fprintf('Scattering cross section is set to zero!\n');
    fprintf('\n');
    
end

%% Geometry and Angular Discretization
xL = 0; xR = 1;
Nx = 1000;
Nang = 2;

[mui,wi] = lgwt(Nang,-1,1);

angflux = zeros(Nx,Nang);
half_angular_flux = zeros(Nx+1,Nang);

dx = (xR - xL)/Nx;
x = xL+dx/2:dx:xR-dx/2;

%% Source (constant through problem) 
qx = 1.*ones(1,length(x));


%% Transport Sweep

% Check calculation of scattering source. Potential source of scattering
% errors.

for iter = 1:1000
    
    angfluxprev = angflux;
    
    for l = 1:Nang
        
        if ( mui(l) > 0 )
            
            for i = 1:Nx
                
                q(i) = qx(i);
            
                angflux(i,l) = ( 1 + sigt*dx/(2*abs(mui(l))) )^(-1)*...
                    ( half_angular_flux(i,l) + dx*q(i)/((2*abs(mui(l)))));
                
                half_angular_flux(i+1,l) = 2*angflux(i,l) - half_angular_flux(i,l);
                
            end
            
        else
            
            for i = Nx:-1:1
                
                %q(i) = 0.5*sigs0*angflux(i,:)*wi + qx(i);
            
                angflux(i,l) = ( 1 + sigt*dx/(2*abs(mui(l))) )^(-1)*...
                    ( half_angular_flux(i+1,l) + dx*q(i)/((2*abs(mui(l)))));
                
                half_angular_flux(i,l) = 2*angflux(i,l) - half_angular_flux(i+1,l);
                
            end
            
        end
        
    end
    
    residual = norm(angflux-angfluxprev);
    
    if residual < 1e-10;
        
        fprintf('Scalar flux converged in %i iterations with residual %f\n',iter,residual);
        
        break
        
    end
    
end

phi = zeros(1,length(x));

for i = 1:Nang
    
    phi = phi + 0.5*wi(i).*angflux(:,i)';
    
end

%plot(x,phi);
%max(phi)

hold on

psi_mu_pos = @(x,mu) (q/sigt).*(1-exp(-sigt.*x./mu));
psi_mu_neg = @(x,mu) (q/sigt).*(1-exp((sigt/mu).*(xR-x)));

% for i = 1:Nang
%     
%     plot(x,angflux(:,i),'o')
%     
%     if mui(i) > 0
%         
%         plot(x,psi_mu_pos(x,mui(i)))
%         err(i) = max(abs(angflux(:,i)-psi_mu_pos(x,mui(i))'));
%         
%     else
%         
%         plot(x,psi_mu_neg(x,mui(i)))
%         err(i) = max(abs(angflux(:,i)-psi_mu_neg(x,mui(i))'));
%         
%     end
%     
% end
% % 
% %maxloc = find(max(err) == err);
% %err(maxloc);
% 
% phi_anal = zeros(1,length(x));
% 
% for i = 1:Nang
%     
%     if mui(i) > 0
%         
%         phi_anal = phi_anal + wi(i).*psi_mu_pos(x,mui(i));
%         
%     else
%         
%         phi_anal = phi_anal + wi(i).*psi_mu_neg(x,mui(i));
%         
%     end
%     
% end
% 
% phi_anal = 0.5.*phi_anal;
% 
% max(phi)
% max(phi_anal)