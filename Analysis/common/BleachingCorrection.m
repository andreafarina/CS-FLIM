function [ref] = BleachingCorrection(spc,ref)
        N = 2;
        CWs = reshape([sum(spc(:,:,3:7),3) sum(spc(:,:,(end-4):end),3)],[size(spc,1) size(spc,2) N]);
        [t,l,k] = ndgrid(1:size(spc,1),1:size(spc,2),linspace(1,size(spc,3),N));
        [tq,lq,kq] = ndgrid(1:size(spc,1),1:size(spc,2),1:size(ref,3));
        CW = interpn(t,l,k,CWs,tq,lq,kq);
        CW = sum(CW,1);
        CW = CW(1,:,:)./CW(1,:,1);
        ref = ref./CW;
end