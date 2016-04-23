function [x,scalar_flux] = OneDNeutronTransportSolver(Nx,Nang,sigt,sigs0,bc,accel)

if ( sigs0 == 0 )
    
    fprintf('Scattering cross section is set to zero!\n');
    fprintf('\n');
    
end

xL = 0; xR = 1;

[mui,wi] = lgwt(Nang,-1,1);

angflux = zeros(Nx,Nang);
half_angular_flux = zeros(Nx+1,Nang);
half_angular_flux(:,1) = 0;

dx = (xR - xL)/Nx;
x = xL+dx/2:dx:xR-dx/2;

%% Source (constant through problem) 
qx = 1.*ones(1,length(x));

%% Boundary Conditions
if ( bc == 0 )
    
    fprintf('Vacuum boundary conditions \n');
    fprintf('\n');
    
elseif ( bc == 1 )
    
    fprintf('Reflective boundary condition on right boundary \n')
    fprintf('\n');
    
elseif ( bc == 2 )
    
    fprintf('Reflective boundary condition on both boundaries \n')
    fprintf('\n');
    
end

%% Acceleration Initialization
if ( accel == 1 )
    
    fprintf('Diffusion synthetic acceleration turned on \n');
    fprintf('\n');
    
elseif ( accel == 2 )
    
    fprintf('Transport synthetic acceleration turned on \n');
    fprintf('\n')
    
end

%% Transport Sweep
for iter = 1:1e6
    
    angfluxprev = angflux;
    
    scalar_flux = zeros(1,length(x));
    
    for k = 1:Nang
        
        scalar_flux = scalar_flux + 0.5.*wi(k).*angflux(:,k)';
        
    end
    
    if ( accel == 1 && iter > 1 )
        
        scalar_flux = scalar_flux + f';
        
    end
    
    for l = 1:Nang
        
        if ( mui(l) > 0 )
            
            for i = 1:Nx
                
                if ( bc == 2 )
                    
                    half_angular_flux(1,l) = half_angular_flux(1,Nang+1-l);
                    
                end
                
                q = qx + sigs0.*scalar_flux;
            
                angflux(i,l) = ( 1 + sigt*dx/(2*abs(mui(l))) )^(-1)*...
                    ( half_angular_flux(i,l) + dx*q(i)/((2*abs(mui(l)))));
                
                half_angular_flux(i+1,l) = 2*angflux(i,l) - half_angular_flux(i,l);
                
            end
            
        else
            
            for i = Nx:-1:1
                
                if ( bc == 1 )
                    
                    half_angular_flux(Nx+1,l) = half_angular_flux(Nx+1,Nang+1-l);
                    
                end
                
                if ( bc == 2 )
                    
                    half_angular_flux(Nx+1,l) = half_angular_flux(Nx+1,Nang+1-l);
                    
                end
                
                q = qx + sigs0.*scalar_flux;
            
                angflux(i,l) = ( 1 + sigt*dx/(2*abs(mui(l))) )^(-1)*...
                    ( half_angular_flux(i+1,l) + dx*q(i)/((2*abs(mui(l)))));
                
                half_angular_flux(i,l) = 2*angflux(i,l) - half_angular_flux(i+1,l);
                
            end
            
        end
        
    end
    
    residual = norm(angflux-angfluxprev);
    
    fprintf('Residual: %e     Iteration: %i \n',residual,iter-1);
    
    if residual < 1e-10;
        
        fprintf('Scalar flux converged in %i iterations with residual %e\n',iter-1,residual);
        
        break
        
    end
    
    phi_prev = zeros(1,Nx);
    phi_new = zeros(1,Nx);
    
    for k = 1:Nang
        
        phi_prev = phi_prev + 0.5*wi(k).*angfluxprev(:,k)';
        phi_new = phi_new + 0.5*wi(k).*angflux(:,k)';
        
    end
    
    if ( accel == 1 )
        
        [S,f] = DiffusionSyntheticAccel_FiniteDiff(Nx,dx,sigt,sigs0,sigt-sigs0,phi_prev,phi_new);
        %f = DiffusionSyntheticAccel(Nx,dx,sigt,sigs0,sigt-sigs0,phi_prev,phi_new);
        
    end
    
end

scalar_flux = phi_new;

return