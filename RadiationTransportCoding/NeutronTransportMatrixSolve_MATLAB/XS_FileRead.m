function XS = XS_FileRead(file)

fid = fopen(file);

numvals = textscan(fid,'%d',1);
numvals = numvals{1};

XS = textscan(fid,'%f',numvals);
XS = XS{1};

fclose(fid);

return