function TzMG = GenerateMGTzMatrix(M,L,sigs0,sigs1,sigs2)

G = length(sigs0);

SigmaS = GenerateMGScatteringMomentMatrix(M,sigs0,sigs1,sigs2);

[L0,L0p,L1,L1p,L2,L2p] = GenerateFluxMomentMatrices(M,L);

LN = kron(eye(G,G),[L0;L1;L2]);
LNplus = kron(eye(G,G),[L0p L1p L2p]);

TzMG = LNplus*SigmaS*LN;

return