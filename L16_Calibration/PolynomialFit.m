function [lambda,p_baric] = PolynomialFit(lambda0,spectrum,order)
% polynomial fit of the wavelengths based on the baricenter of the
% time-resolved multichannel data
% lambda0:  nominal wavelengths
% spectrum: spectrum obtained by the multichannel detector
% order:    order of the polynomial fit (typical 1).
% lambda:   calbrated lambda
% p_baric:  baricenter position in channel units.

Nchan = length(spectrum);
p = 1:Nchan;
p_baric = (p*spectrum)./sum(spectrum);

%% find polinomial function
c = polyfit(p_baric,lambda0,order);

lambda = polyval(c,p);

%figure(2),hold on, plot(lambda0,p_baric,'o',lambda,p),xlabel('lambda (nm)'),
%ylabel('channel')