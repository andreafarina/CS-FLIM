function [p] = patternSin_sequence(L,npatt,nphase)
x = [1:L]; 
y = [1:L]; 
[xx,yy] = meshgrid(x,y); 
p = zeros(L,L,npatt*nphase); 
 
for i = 1:npatt*nphase
    Tx = L/(floor((i-1)/nphase));
    phase = 2*pi/nphase*(i-1);
    p(:,:,i) = 0.5 + 0.5*sin(2*pi/Tx*xx + phase); 
     imshow(p(:,:,i)),pause
end
 
%p = im2uint8(p); 
 
% figure
% for i = 1:npatt*nphase
%     subplot(npatt,nphase,i),imshow(p(:,:,i))
% end
 
end
