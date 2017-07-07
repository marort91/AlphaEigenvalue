function SigmaS = GenerateMGScatteringMomentMatrix(M,sigs0,sigs1,sigs2)

Egrp = length(sigs0);

Gij = numel(sigs0);
sigs0 = sigs0'; sigs0 = sigs0(:)';
sigs1 = sigs1'; sigs1 = sigs1(:)';
sigs2 = sigs2'; sigs2 = sigs2(:)';


for g = 1:Gij
    
    S0{g} = spdiags(sigs0(g).*ones(M,1),0,M,M);
    S1{g} = spdiags(sigs1(g).*ones(M,1),0,M,M);
    S2{g} = spdiags(sigs2(g).*ones(M,1),0,M,M);
    
end

for g = 1:Gij
    
    Sgg{g} = blkdiag(S0{g},S1{g},S2{g});
    
end

%G = sqrt(length(sigs0));
 
for g = 1:Egrp
    
    Sggtmp{g} = horzcat(Sgg(1+(g-1)*Egrp:g*Egrp));
     
end

SigmaS = cell2mat(vertcat(Sggtmp{:}));

return