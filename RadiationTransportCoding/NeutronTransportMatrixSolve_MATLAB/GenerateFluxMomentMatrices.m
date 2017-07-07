function [L0,L0p,L1,L1p,L2,L2p] = GenerateFluxMomentMatrices(M,L)

[mu,w] = AngularQuad1DSlab(L,-1,1);
mu = sort(mu,'ascend');

Im = speye(M,M);
L0 = 0.5*w(1)*Im;

%Zeroth Moment

for el = 2:L
    
    L0 = [ L0, 0.5*w(el)*Im ];
    
end

L0p = Im;

for el = 2:L
    
    L0p = [ L0p; Im ];
    
end

%First Moment
l1 = mu';
W = diag(w);
L1 = kron((l1*W),Im);

%L1p = 3.*kron(l1',Im);
L1p = 1.*kron(l1',Im);

%Second Moment
l2 = 0.5.*(3.*(mu').^2-1);
W = diag(w);
L2 = kron((l2*W),Im);

%L2p = 5.*kron(l2',Im);
L2p = 1.*kron(l2',Im);

return