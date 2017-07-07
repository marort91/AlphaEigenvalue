function [xi,HzMG,invVzMG,mu] = GenerateMGHzMatrix(xL,xR,M,L,sigt,V)

G = length(diag(sigt));

sigt = diag(sigt);
V = diag(V);

%Spatial Discretization
xi = xL:(xR-xL)/(M):xR;
dxi = max(diff(xi));

%Angular Discretization
[mu,w] = AngularQuad1DSlab(L,-1,1);
mu = sort(mu,'ascend');

%Z Matrix Generation
Z = sparse(M,M);
Zm = Z(M,:)';
Z = [speye(M,M); sparse(1,M)];
Zb = sparse(M+1,1); Zb(end) = 1;

%Derivative Matrix Generation
Dx = sparse(diag(dxi.*ones(M,1)));
Dm = sparse(-[diag(ones(M,1),0),Zm] + [Zm,diag(ones(M,1),0)]);
Cx = Dx\Dm;

%Step Differencing
Splus = 1.0.*sparse([Zm,diag(ones(M,1),0)]);
Sminus = 1.0.*sparse([diag(ones(M,1),0),Zm]);

%Generation of Basis Vectors
e = speye(M+1);
e0m = e(1,:)';
emm = e(M+1,:)';
B1 = e0m';
B2 = emm';

for i = 1:length(mu)
    
    if ( mu(i) < 0 )
        
        S{i} = Sminus;
        B{i} = B2;
        
    else
        
        S{i} = Splus;
        B{i} = B1;
        
    end
    
    Cl{i} = mu(i)*Cx;
    
end

S = blkdiag(S{:});
C = blkdiag(Cl{:});
B = blkdiag(B{:});
Zbbar = kron(speye(L,L),Zb);

for g = 1:G

    Sigma_Bar{g} = kron(speye(L,L),sigt(g).*speye(M,M));
    invV{g} = kron(speye(L,L),1./V(g).*speye(M,M));
    
end

Sigma_Bar = blkdiag(Sigma_Bar{:});
invV_Bar = blkdiag(invV{:});

% Extension to Multigroup
Zbar = kron(eye(L,L),Z);
Z = kron(eye(G,G),Zbar);

C = kron(eye(G,G),C);
B = kron(eye(G,G),B);
Zbbar = kron(eye(G,G),Zbbar);
Sbar = kron(eye(G,G),S);

Cz = C*inv(Z*Sbar + Zbbar*B)*Z;

HzMG = Sigma_Bar + Cz;

invVzMG = Z'*Z*invV_Bar*Sbar*inv(Z*Sbar + Zbbar*B)*Z;

return