clc, clear, clf, close all

nusigf = [ 0.070 0.070;
           0.000 0.000];
       
sigs = [ 0.000 0.000;
         0.050 0.000];
     
sigR = [ 0.09 0.09 ];
D = [0.9 0.9];

%B = @(x) 1 - nusigf(1)/(sigR(1)+D(1)*x) + sigs(2,1)/(sigR(1)+D(1)*x)*(nusigf(2)/(sigR(2)+D(2)*x));
%B2 = fsolve(B,50);

Lx = 30;

nEgrps = 2;

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
    
    for j = 1:nEgrps
        
        for k = 1:nEgrps
            
            f(:,1,j) = f(:,1,j) + (1/keff).*nusigf(j,k)*phi(:,1,k);
            
        end
        
        for k = 1:nEgrps
                
            q(:,1,j) = q(:,1,j) + sigs(j,k)*phi(:,1,k);
            
        end
            
        psi(:,1,j) = M(:,:,j)\(f(:,1,j)+q(:,1,j));
            
    end
    
    Sprev = 0;
    Snew = 0;
    
    for l = 1:nEgrps
        
        for m  = 1:nEgrps
            
            Snew = Snew + sum(nusigf(l,m)*psi(:,1,m));
            Sprev = Sprev + sum(nusigf(l,m)*phi(:,1,m));
        
        end
        
    end
    
    keff = Snew/((1/keff)*Sprev);
    
    phi = psi;
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
    
    plot(x,phi(:,1,i),strs{i});
    
end
legend('Fast Flux','Thermal Flux');
grid on