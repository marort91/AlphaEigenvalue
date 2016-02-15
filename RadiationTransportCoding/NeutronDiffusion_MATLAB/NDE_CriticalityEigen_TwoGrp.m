clc, clear, clf, close all

nusigf = [ 0.080 0.080;
           0.000 0.000];
       
sigs = [ 0.000 0.000;
         0.050 0.000];
     
sigR = [ 0.09 0.09 ];
D = [0.9 0.9];

%nusigf = [ 0.070 0.010 ];
%sigtot = [0.066 0.066];
%sigtot = [ 0.1 0.1 ];
%sigss = [0.010 0.010];
%sigR = sigtot - diag(sigs)';
%sigR = [ 0.066 1 ];
%D = [ 0.9 0.9 ];
%sigs = [ 0.0 0.05 ];

% nusigf = [ 0.008476 0.18514 ];
% sigR = [ 0.02619 0.1210 ];
% D = [ 1.2627 0.3543 ];
% siga = [ 0.01207 0.1210];
% 
% sigs = [ sigR(2)-siga(2) sigR(1) - siga(1) ];

%B = @(x) 1 - nusigf(1)/(sigR(1)+D(1)*x) + sigs(2)/(sigR(1)+D(1)*x)*(nusigf(2)/(sigR(2)+D(2)*x));
%B2 = fsolve(B,50);

%Lx = pi/(sqrt(B2))

%This makes slab critical
%Lx = 33.401;
%Lx = 6.79405;
nEgrps = 2;

Lx = 30.0;
%Lx = pi*((nusigf(1)-sigR(1))/D(1))^(-0.5);
N = 100;

h = Lx/(N-1);

x = 0:h:Lx;

L = zeros(N,N,nEgrps);

for i = 1:nEgrps

    L(:,:,i) = (-D(i)/h^2).*full(gallery('tridiag',N,1,-2,1));
    
end

for i = 1:nEgrps
    
    T(:,:,i) = sigR(i).*eye(N,N);
    
end

M = L + T;

for i = 1:nEgrps
    
    M(1,1,i) = 1; M(1,2,i) = 0;
    M(end,end,i) = 1; M(end,end-1,i) = 0;
    
end

phi = ones(N,1,nEgrps);
psi = ones(N,1,nEgrps);
q = zeros(N,1,nEgrps);
f = zeros(N,1,nEgrps);

for i = 1:nEgrps
    
    phi(1,1,i) = 0;
    phi(end,1,i) = 0;
    
end

tol = 1e-15;
keff = 1.0;

for i = 1:1000
    
    k_i = keff;
    f(:,:,:) = 0;
    q(:,:,:) = 0;
    %q(:,:,:) = 0;
    
    for j = 1:nEgrps
        
        % Need to separate fission source and scattering source for my
        % sake.
        
        for k = 1:nEgrps
            
            f(:,1,j) = f(:,1,j) + nusigf(j,k)*phi(:,1,k);
            
        end
        
        for k = 1:nEgrps
            
            %if k ~= j
                
                q(:,1,j) = q(:,1,j) + sigs(j,k)*phi(:,1,k);
                
            %end
            
        end
            
        psi(:,1,j) = M(:,:,j)\(f(:,1,j)+q(:,1,j));
            
    end
    
    Sprev = 0;
    Snew = 0;
    
    for l = 1:nEgrps
        
        for m  = 1:nEgrps
            
            Snew = Snew + sum(nusigf(l,m)*psi(:,1,m));
            Sprev = Sprev + sum(nusigf(l,m)*phi(:,1,m));
        
        %Snew = Snew + sum(nusigf(l).*psi(:,1,l));
        %Sprev = Sprev + sum(nusigf(l).*phi(:,1,l));
        
        end
        
    end
    
    keff = Snew/Sprev;
    
    phi = (1/keff).*psi;
    phi(1,1,:) = 0;
    phi(N,1,:) = 0;
    
    residual = norm(abs(keff-k_i));
    
    if residual <= tol
        
        fprintf('k-effective converged: %f in %i iterations \n',keff,i)
        break
        
    end
    
end

hold on

strs = {'ko-','rx-'};

for i = 1:nEgrps
    
    plot(x,phi(:,1,i),strs{i})%,x,phi(:,1,2),'rx-');
    
end
legend('Fast Flux','Thermal Flux');
grid on