function [had] = PosNeg2Pos(spc,N,method)
% spc shall be a matrix defined spc(time,lambda,pattern number)

% N is the number of negative pattern to be use to reconstruct CWs

% method is a way to use CWs: if I have 20 CWs I can do a mean, or
% an interpolation between them.

indexCW = zeros(1,size(spc,3));
indexCW(1:2:end) = 1;
if N == 1
    neg2save = 2;
else
    [~,d] = find((indexCW == 0));
    b = round(linspace(1,numel(d),N));
    neg2save = d(b);
end

CWs = spc(:,:,neg2save-1) + spc(:,:,neg2save);

if strcmp(method,'Mean') || length(neg2save) == 1
    CW = mean(CWs,3);
else
    if size(spc,1) == 1 && size(spc,2) == 1
        CW = permute(interp1(neg2save/2, squeeze(CWs),1:size(spc,3)/2),[3 1 2]);
    elseif size(spc,1) == 1 || size(spc,2) == 1
        p = find(size(spc)~=1);
        [l,k] = meshgrid(neg2save/2,1:size(spc,p(1)));
        [lq,kq] = meshgrid(1:size(spc,3)/2,1:size(spc,p));
        CW = interp2(l,k,squeeze(CWs),lq,kq);
        if size(spc,1) == 1
            CW = permute(CW,[3 1 2]);
        else
            CW = permute(CW,[1 3 2]);
        end
    else
        [t,l,k] = ndgrid(1:size(spc,1),1:size(spc,2),neg2save/2);
        [tq,lq,kq] = ndgrid(1:size(spc,1),1:size(spc,2),1:size(spc,3)/2);
        CW = interpn(t,l,k,CWs,tq,lq,kq);
    end
end

had = 2.*spc(:,:,1:2:end)-CW;

end

