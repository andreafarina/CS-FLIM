function Data = PrepareDataFusion(im,camera,lambda,t,thresh,imGT)

if size(camera,1) > 128
    sizeC = [256 256];
else
    sizeC = size(camera);
end
        
            Data.CCD = imresize(camera,sizeC,'box');

            Data.PMT = permute(im,[2 1 3 4]);
            
            MaskCCD = Data.CCD>thresh;

            Data.CCD = Data.CCD.*MaskCCD;
            Data.CCD = Data.CCD./sum(Data.CCD(:));
            MaskPMT = imresize(MaskCCD,size(Data.PMT,3:4));
            MaskPMT = permute(repmat(MaskPMT,[1,1,size(Data.PMT,1),size(Data.PMT,2)]),[3,4,1,2]);
            Data.PMTnoMask = Data.PMT;
            Data.PMT = Data.PMT.*MaskPMT;
            Data.PMT = Data.PMT./sum(Data.PMT(:)).*sum(Data.CCD(:));
            Data.MaskCCD = permute(repmat(MaskCCD,[1,1,size(Data.PMT,1),size(Data.PMT,2)]),[3,4,1,2]);
            if nargin>5
                Data.GT = permute(imGT,[2,1,3,4]);
                Data.GT = Data.GT/sum(Data.GT(:))*sum(Data.CCD(:));
                 
            end

            Data.time = t;
            Data.lambda = lambda;
end