function [had] = PosNeg2Had(spc)
% spc shall be a matrix defined spc(time,lambda,pattern number)
    had = spc(:,:,1:2:end) - spc(:,:,2:2:end);
end

