function F = GenerateMGFissionMatrix(M,F)

Egrp = length(F);

Gij = numel(F);
F = F'; F = F(:)';

for g = 1:Gij
    
    F0{g} = spdiags(F(g).*ones(M,1),0,M,M);
    F1{g} = spdiags(0.*ones(M,1),0,M,M);
    F2{g} = spdiags(0.*ones(M,1),0,M,M);
    
end

for g = 1:Gij
    
    Fgg{g} = blkdiag(F0{g},F1{g},F2{g});
    
end

%G = sqrt(length(sigs0));
 
for g = 1:Egrp
    
    Ftmp{g} = horzcat(Fgg(1+(g-1)*Egrp:g*Egrp));
     
end

F = cell2mat(vertcat(Ftmp{:}));

return