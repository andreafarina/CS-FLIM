function [had] = PosNeg2Pos(spc,N,method)
% spc shall be a matrix defined spc(time,lambda,pattern number)

% N is the number of negative pattern to be use to reconstruct CWs

% method is a way to use CWs: if I have 20 CWs I can do a mean, or
% interpolate between them. 

[~,p] = max(sum(spc,1:2));
CWs = spc;
CWs(:,:,p:(p+1)) = [];
CWs = PosNeg2Had(CWs);
pos_number = round(linspace(0,size(CWs,3),2*N+1));
pos_number = pos_number(2:2:end);

spc_1 = spc(:,:,1:2:end);
CWs = CWs(:,:,pos_number);
CW = Pattern_reconstruction(CWs,spc_1,pos_number,method);
had = 2.*spc_1-CW;

end


function [CW] = Pattern_reconstruction(CWs,spc,pos_number,method)
if strcmp(method,'Mean')
    CW = mean(CWs,3);
else
    D = round(diff(pos_number));
    P = [];
    for l = 1:length(D)
        P = [P repmat(l,[1 D(l)])];
    end
    P = [P repmat(length(pos_number),[1 (size(spc,3) - size(P,2))])];
    CW = CWs(:,:,P);
end
end
