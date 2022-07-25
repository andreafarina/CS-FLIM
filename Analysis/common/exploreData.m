function exploreData(camera,im,GT)
%im[time,lambda,x,y]
figure
while true
subplot(2,3,1), imshow(squeeze(camera),[],'Colormap',parula); axis square; axis on; colorbar; title('CCD image');
subplot(2,3,2),imshow(squeeze(sum(sum(im,1),2)),[],'Colormap',parula); axis square; axis on; colorbar; title('Intensity image');
subplot(2,3,3),plot(squeeze(sum(im,[1,3,4])))
subplot(2,3,6),plot(squeeze(sum(im,[2,3,4])))
[x,y] = ginput(1);
row = round(y); col = round(x);
if row<0 ||col<0
    break;
end
slicetime = squeeze(sum(im(:,:,row,col),2));
slicelambda = squeeze(sum(im(:,:,row,col),1));
if nargin>2
    GTslicetime = squeeze(sum(GT(:,:,row,col),2));
    GTslicelambda = squeeze(sum(GT(:,:,row,col),1));
    plotTime = [slicetime,GTslicetime];
    plotLambda = [slicelambda',GTslicelambda'];
    strlegend = {'data','GT'};
else
    plotTime = slicetime';
    plotLambda = slicelambda;
    strlegend = {'data'};
end

hold on
subplot(2,3,4),plot(plotTime),legend(strlegend)
subplot(2,3,5),plot(plotLambda),legend(strlegend)

% end
end