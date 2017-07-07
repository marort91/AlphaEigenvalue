function XSDatFile_Generate(ZAID,XSDatName)

% XSDatFile_Generate Generates XS DAT File for MATLAB Analysis
% XSDatFile_Generate generates a DAT file containing the cross
%   sections for the material specified by ZAID. Cross sections must be in a
%   directory that ends with a .s extension with a subdirectory with the name
%   of the material ZAID. Your current directory must be the directory
%   containing the cross section directory.
%
% Cross Section Directory: example.s
% Cross Section Subdirectory: ZAID
%
% The directory must contain the following file: enub
%
% The subdirectory must contain the following files:
% chi, p0, p1, p2, p3, p4, p5, p6, p7, siga, sigfnu, sigt

ZAID = num2str(ZAID);
XSdir = dir('*.s');

%CD into Cross Section Directory
cd(XSdir.name)

%Open Energy Boundary File
enub_file = 'enub';
enub = wrev(XS_FileRead(enub_file));
Egrp = length(enub)-1; %Energy groups

%CD in Isotope Directory
cd(ZAID)

%Open Fission Spectra (chi) File
chi_file = 'chi';
chi = wrev(XS_FileRead(chi_file));

%Open p0 Scattering (p0) File
p0_file = 'p0';
p0 = fliplr(flip(reshape((XS_FileRead(p0_file)),[Egrp Egrp])));

%Open p1 Scattering (p1) File
p1_file = 'p1';
p1 = fliplr(flip(reshape((XS_FileRead(p1_file)),[Egrp Egrp])));

%Open p2 Scattering (p2) File
p2_file = 'p2';
p2 = fliplr(flip(reshape((XS_FileRead(p2_file)),[Egrp Egrp])));

%Open p3 Scattering (p3) File
p3_file = 'p3';
p3 = fliplr(flip(reshape((XS_FileRead(p3_file)),[Egrp Egrp])));

%Open p4 Scattering (p4) File
p4_file = 'p4';
p4 = fliplr(flip(reshape((XS_FileRead(p4_file)),[Egrp Egrp])));

%Open p5 Scattering (p5) File
p5_file = 'p5';
p5 = fliplr(flip(reshape((XS_FileRead(p5_file)),[Egrp Egrp])));

%Open p6 Scattering (p6) File
p6_file = 'p6';
p6 = fliplr(flip(reshape((XS_FileRead(p6_file)),[Egrp Egrp])));

%Open p7 Scattering (p7) File
p7_file = 'p7';
p7 = fliplr(flip(reshape((XS_FileRead(p7_file)),[Egrp Egrp])));

%Open Absorption (siga) File
siga_file = 'siga';
siga = diag(wrev(XS_FileRead(siga_file)));

%Open Fission (sigfnu) File
sigfnu_file = 'sigfnu';
sigfnu = wrev(XS_FileRead(sigfnu_file));

%Open Total (sigt) File
sigt_file = 'sigt';
Sigma = diag(wrev(XS_FileRead(sigt_file)));

%Generate Velocity Matrix
mN = 939.5654133; %Mass of Neutron (MeV/c2)
c = 2.99792458e10; %Speed of light (cm/s)

V = zeros(Egrp,Egrp);

for e = 1:Egrp
    
    dE = abs(enub(e+1)-enub(e));
    Emid = mean([enub(e) enub(e+1)]);
    ep = Emid + 0.5*dE;
    em = Emid - 0.5*dE;
    V(e,e) = (1e-6*c*dE)/(sqrt(ep*(ep+2*mN)) - sqrt(em*(em+2*mN)));
    
end 

cd('../..')

savefile = strcat(XSDatName,'.mat');

%Scatter Matrix S
%S = p0;% + p1 + p2 + p3 + p4 + p5 + p6 + p7;
S0 = p0; S1 = p1; S2 = p2; S3 = p3; S4 = p4; S5 = p5; S6 = p6; S7 = p7;
F = chi*sigfnu';

% save(savefile,'S','F','V','Sigma')
save(savefile,'S0','S1','S2','S3','S4','S5','S6','S7','F','V','Sigma')

return