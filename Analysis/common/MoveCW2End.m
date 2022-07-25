function [had,M] = MoveCW2End(had,M,q)
% We saw that reconstruction with the CW pattern affects the PSNR: there is
% a step in PSNR if, by tuning the CR, the CW pattern is included or not in
% the had vector.

M = M([q(1) 1:(q(1)-1) (q(1)+1):size(M,1)],:); 
had = had([q(1) 1:(q(1)-1) (q(1)+1):length(had)],:);

end