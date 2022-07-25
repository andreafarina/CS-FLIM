function [im_TVAL,opts] = InvertwTVAL3(had,m,M,imA,K,mu,beta)

PLOT_ITERATIONS = 1;

% Rescale reference image from CCD to fit scale of SPC dataset

imA = imresize(imA,[K K],'box');
imA = max(sum(had,1:2),[],3).*imA./sum(imA(:));

% Best mu and b searching algorithm:
%alpha = 0;%0.5;
%imA = imA+alpha.*sqrt(mean(imA(:)))*randn(size(imA));
%opts = SearchTVALMuB(sum(had,1:2),m,M,imA,K);

opts.mu = mu; %default value, suggested by authors
opts.beta = beta; %default value, suggested by authors
opts.maxit = 3000;
opts.tol = 1E-8;
opts.isreal = true;
%% Display
if PLOT_ITERATIONS == 1
    [~, pPix] = min(abs(imA(:)-0.5*max(imA(:))));
    [rPix,cPix] = ind2sub([K K],pPix);
    normalize = @(X)X./max((X));
    c = figure;
    subplot(1,3,1), h1 = imagesc(imA); axis image
    subplot(1,3,2), h2 = plot(squeeze(sum(had,[1 3]))); axis square
    subplot(1,3,3), h3 = plot(squeeze(sum(had,2:3))); axis square
end
%% Solve
[~,b] = max(squeeze(sum(had,1:2)));
for ti = 1:size(had,1)
    parfor li = 1:size(had,2)
        [had_rec, M1] = MoveCW2End(squeeze(had(ti,li,:)),M,b);  
        [U, ~] = TVAL3(M1(1:m,:),had_rec(1:m),K,K,opts);
        if sum(isnan(U(:)))>0
        res(ti,li,:,:) = zeros(K,K);
        else   
        res(ti,li,:,:) = U;
        end
    end
    if PLOT_ITERATIONS == 1
      h1.CData = squeeze(sum(res,1:2));
      h2.YData = normalize(squeeze(sum(res(:,:,rPix,cPix),[1 3:4])));
      h3.YData = normalize(squeeze(sum(res(:,:,rPix,cPix),[2 3:4])));
      drawnow;
    end
    disp(['Processed: ',num2str(100*ti/size(had,1)),'%']);
end
im_TVAL = res;
close(c);
end