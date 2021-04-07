function [data,lambda] = ReadSpectrum_HAM(filename,filenamedark)
% read spectral data from HAMAMATSU COmpact Spectrometer
% A. Farina - CNR
% 2021/04/03
DARK = 0;
if nargin > 1
    DARK = 1;
end

y = readmatrix(filename,'NumHeaderLines',1);
lambda = y(:,1);
data = y(:,2);

if DARK > 0
    y = readmatrix(filenamedark,'NumHeaderLines',1);
    dataDARK = y(:,2);
    data = data - dataDARK;
end
    

