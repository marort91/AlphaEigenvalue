%% Two-dimension, One-speed Discrete Ordinates Neutron Transport Equation Solver
% Solves the 2D neutron transport code using the discrete ordinates method.

clc, clear, clf

%% Material Properies
sigt = 0.1;
sigs0 = 0.0;
siga = sigt - sigs0;

%% Geometry and Angular Discretization
xL = 0; xR = 1;
yB = 0; yT = 1;

Nx = 100; Ny = 100;
Nang = 2;

[mu,eta,wi] = level_sym_table(Nang);

dx = (xR - xL)/Nx;
dy = (yT - yB)/Ny;

%% Generation of Angular Flux Arrays
angular_flux = zeros(Nx,Ny,Nang*(Nang+2)/2);
angular_flux_half_x = zeros(Nx,Nang*(Nang+2)/2);
angular_flux_half_y = zeros(Ny,Nang*(Nang+2)/2);

%% Source Generation
S = 1.*ones(Nx,Ny);
Q = zeros(Nx,Ny);

%% Calculation Parameters
maxiter = 1e3; %Maximum number of sweeps allowed
tol = 1e-8;

%% Transport Sweep
for iter = 1:maxiter
    
    angular_flux_prev = angular_flux;
    
    scalar_flux = zeros(Nx,Ny);
    
    for k = 1:Nang*(Nang+2)/2
        
        scalar_flux = scalar_flux + 0.25.*wi(k).*angular_flux(:,:,k);
        
    end
    
    Q = ( S + sigs0.*scalar_flux );
    
    for l = 1:Nang*(Nang+2)/2
        
        if ( mu(l) > 0 && eta(l) > 0 )
            
            for j = 1:Ny
                
                for i = 1:Nx
                    
                    angular_flux(i,j,l) = ( 2*mu(l)*angular_flux_half_x(i,l)/dx + ...
                        2*eta(l)*angular_flux_half_y(j,l)/dy + Q(i,j) )/...
                        ( 2*mu(l)/dx + 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(i+1,l) = 2*angular_flux(i,j,l) - angular_flux_half_x(i,l);
                    
                end
                
                angular_flux_half_y(j+1,l) = 2*angular_flux(i,j,l) - angular_flux_half_y(j,l);
                
            end
            
        elseif ( mu(l) < 0 && eta(l) > 0 )
            
            for j = 1:Ny
                
                for i = Nx:-1:1
                    
                    angular_flux(i,j,l) = ( -2*mu(l)*angular_flux_half_x(i+1,l)/dx + ...
                        2*eta(l)*angular_flux_half_y(j,l)/dy + Q(i,j) )/...
                        ( -2*mu(l)/dx + 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(i,l) = 2*angular_flux(i,j,l) - angular_flux_half_x(i+1,l);
                    
                end
                
                angular_flux_half_y(j+1,l) = 2*angular_flux(i,j,l) - angular_flux_half_y(j,l);
                
            end
                    
        elseif ( mu(l) > 0 && eta(l) < 0 )
            
            for j = Ny:-1:1
                
                for i = 1:Nx
                    
                    angular_flux(i,j,l) = ( 2*mu(l)*angular_flux_half_x(i,l)/dx - ...
                        2*eta(l)*angular_flux_half_y(j+1,l)/dy + Q(i,j) )/...
                        ( 2*mu(l)/dx - 2*eta(l)/dy + sigt );
                    
                    angular_flux_half_x(i+1,l) = 2*angular_flux(i,j,l) - angular_flux_half_x(i,l);
                    
                end
                
                angular_flux_half_y(j,l) = 2*angular_flux(i,j,l) - angular_flux_half_y(j+1,l);
                
            end
            
        elseif ( mu(l) < 0 && eta(l) < 0 )
            
            for j = Ny:-1:1
                
                for i = Nx:-1:1
                    
                    angular_flux(i,j,l) = ( -2*mu(l)*angular_flux_half_x(i+1,l)/dx - ...
                        2*eta(l)*angular_flux_half_y(j+1,l)/dy + Q(i,j) )/...
                        ( -2*mu(l)/dx - 2*eta(l)/dy - sigt );
                    
                    angular_flux_half_x(i,l) = 2*angular_flux(i,j,l) - angular_flux_half_x(i+1,l);
                    
                end
                
                angular_flux_half_y(j,l) = 2*angular_flux(i,j,l) - angular_flux_half_y(j+1,l);
                
            end
            
        else
            
            error('Angular discretization incorrect \n')
            
        end
        
    end
    
    scalar_flux_prev = zeros(Nx,Ny);
    scalar_flux_new = zeros(Nx,Ny);
    
    for k = 1:Nang*(Nang+2)/2
        
        scalar_flux_prev = scalar_flux_prev + 0.25.*wi(k).*angular_flux_prev(:,:,k);
        scalar_flux_new = scalar_flux_new + 0.25.*wi(k).*angular_flux(:,:,k);
        
    end
    
    residual = norm(scalar_flux_new-scalar_flux_prev);
    fprintf('Residual: %f     Iteration: %i \n',residual,iter-1);
    
    if ( residual < tol )
        
        break
        
    end
       
end