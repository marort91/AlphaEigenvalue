% Implementation of Hansen-Roach Six Group Data.

clc, clear, clf

%Data ordering: fast group -> thermal group
% chi = [0.204 0.344 0.168 0.18 0.09 0.014];
% nusigf = [ 3.557 3.196 3.087 2.988 3.518 5.71];
% sigf = [1.21 1.22 1.22 1.2 1.43 2.34];
% sigc = [ 0.05 0.08 0.11 0.15 0.23 0.6];
% sigtr = [ 4.25 4.5 4.65 5.2 7.9 12];
% 
% sigs = [1.2	0.27 0.37 0.65 0.44 0.06;
%         0.0 1.77 0.24 0.67 0.45 0.07;
%         0.0 0.00 2.30 0.55 0.40 0.07;
%         0.0	0.00 0.00 3.42 0.35 0.08;
%         0.0	0.00 0.00 0.00 6.16	0.08;
%         0.0	0.00 0.00 0.00 0.00	9.06];
%     
% sigst = sum(sigs);
% sigt = sigf + sigc + sigst;
% 
% sigR = sigt - diag(sigs)';
% 
% D = 1./(3.*sigtr);

XSdata = CrossSectionRead(1,{'U235_fast_6grp.txt'});
[chi,nusigf,sigf,D,sigR,sigs] = XSInterpret(XSdata);

nEgrps = length(chi);

alpha = 0;

% Geometry
Lx = 0.553869;
N = 100;
h = Lx/(N-1);

x = 0:h:Lx;

L = zeros(N,N,nEgrps);

for i = 1:nEgrps

    L(:,:,i) = (-D(i)/h^2).*full(gallery('tridiag',N,1,-2,1));
    
end

T = zeros(N,N,nEgrps);

for i = 1:nEgrps
    
    T(:,:,i) = (sigR(i) + alpha).*eye(N,N);
    
end

M = L + T;

for i = 1:nEgrps
    
    M(1,1,i) = 1; M(1,2,i) = 0;
    M(end,end,i) = 1; M(end,end-1,i) = 0;
    
end

F = chi'*nusigf;

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

for i = 1:10000
    
    ki = keff;
    q(:,:,:) = 0;
    f(:,:,:) = 0;
    
    for j = 1:nEgrps
        
        for k = 1:nEgrps
            
            f(:,1,j) = f(:,1,j) + (1/keff).*F(j,k)*phi(:,1,k);
            
        end
        
        for k = 1:nEgrps
            
            if ( k == j )
                
                continue
                
            else
            
            q(:,1,j) = q(:,1,j) + sigs(k,j)*phi(:,1,k);
            
            end
            
        end
        
        psi(:,1,j) = M(:,:,j)\(f(:,1,j)+q(:,1,j));
        
    end
    
    Sprev = 0;
    Snew = 0;
    
    for l = 1:nEgrps
        
        for m  = 1:nEgrps
            
            Snew = Snew + sum(nusigf(m)*psi(:,1,m));
            Sprev = Sprev + sum(nusigf(m)*phi(:,1,m));
        
        end
        
    end
    
    keff = Snew/((1/keff)*Sprev);
    
    phi = psi;
    phi(1,1,:) = 0;
    phi(N,1,:) = 0;
    
    residual = norm(abs(keff-ki));
    
    if residual <= tol
        
        fprintf('k-effective converged: %f in %i iterations \n',keff,i)
        break
        
    end
    
end

hold on
plotstr = {'ko','ro','bo','kx','rx','bx'};

for i = 1:nEgrps
    
    plot(x,phi(:,:,i))%,plotstr{i})
    
end

legend('Grp.1','Grp.2','Grp.3','Grp.4','Grp.5','Grp.6')
xlabel('Distance (cm)')
ylabel('Neutron Flux')