function FzMG = GenerateMGFzMatrix(M,L,F)

G = length(F);

SigmaF = GenerateMGFissionMatrix(M,F);

[L0,L0p,L1,L1p,L2,L2p] = GenerateFluxMomentMatrices(M,L);

LN = kron(eye(G,G),[L0;L1;L2]);
LNplus = kron(eye(G,G),[L0p L1p L2p]);

FzMG = LNplus*SigmaF*LN;

return