function [c] = CountsForDecayPositivityAtT(im,tau,Tbin)

c = tau*(3/2)*sqrt(2*mean(im,'all')*exp(-3)/Tbin)/exp(-3);

end