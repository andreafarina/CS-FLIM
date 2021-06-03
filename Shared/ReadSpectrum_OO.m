function [data,lambda] = ReadSpectrum_OO(filename,filenamedark)
% read spectral data from Ocean Optics Compact Spectrometer
% A. Ghezzi Polimi
% 2021/05/25
DARK = 0;
if nargin > 1
    DARK = 1;
end

fid=fopen(filename);
cdata=textscan(fid,'%f%f','delimiter','\t', 'HeaderLines', 17 );
fclose(fid);
lambda = cdata{1};
data = cdata{2};

if DARK > 0
    fid=fopen(filenamedark);
    cdata=textscan(fid,'%f%f','delimiter','\t', 'HeaderLines', 17 );
    fclose(fid);
    dataDARK = cdata{2};
    data = data - dataDARK;
end






