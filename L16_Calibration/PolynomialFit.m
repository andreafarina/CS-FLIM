function [lambda,p_baric] = PolynomialFit(lambda0,spectrum,order)

Nchan = length(spectrum);
p = 1:Nchan;
p_baric = (p*spectrum)./sum(spectrum);

%% find polinomial function
c = polyfit(p_baric,lambda0,order);

lambda = polyval(c,p);

%figure(2),hold on, plot(lambda0,p_baric,'o',lambda,p),xlabel('lambda (nm)'),
%ylabel('channel')