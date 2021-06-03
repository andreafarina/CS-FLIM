function [spc] = CWHandling(spc,CW_type,pattern_type,CW)

if strcmp(pattern_type,'Pos-Neg')
    if strcmp(CW_type,'Measured')
        return;
    elseif strcmp(CW_type,'Measured, but exclude')
        [~,b] = max(sum(spc,1:2));
        spc = Change2CW(spc,b);
    elseif strcmp(CW_type,'Not measured')
        [~,b(1)] = min(sum(spc,1:2));
        dato = spc;
        dato(:,:,b(1)) = max(sum(spc,1:2));
        [~,b(2)] = min(sum(dato,1:2));
        b = sort(b,'ascend');
        spc = Change2CW(spc,b(1));
    end
elseif strcmp(pattern_type,'Shifted')
    if strcmp(CW_type,'Measured')
        return;
    elseif strcmp(CW_type,'Not measured')
        [~,notCW] = min(squeeze(sum(spc,1:2)));
        spc(:,:,notCW) = CW(:,:,notCW);
    elseif strcmp(CW_type,'Measured, but exclude')
        [~,notCW] = max(squeeze(sum(spc,1:2)));
        spc(:,:,notCW) = CW(:,:,notCW);        
    end
end

end

function [spc] = Change2CW(spc,b)

if size(spc,3) > b+3
    spc(:,:,b) = spc(:,:,b+3) + spc(:,:,b+2);
else
    spc(:,:,b) = spc(:,:,b-2) + spc(:,:,b-1);
end

end