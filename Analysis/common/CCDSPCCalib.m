%% Camera to SPC Calibration
function [] = CCDSPCCalib(Points)
Points(1:9,:) = min(min(Points(10:end,:)));
figure,imagesc(Points),axis image, colorbar, title('Select corners. Draw clockwise:')
points = Points;
h = drawpolygon;

coord(1,1,:) = h.Position(1,:);
coord(1,2,:) = h.Position(2,:);
coord(2,1,:) = h.Position(3,:);
coord(2,2,:) = h.Position(4,:);

d(1) = sqrt((coord(1,1,1)-coord(1,2,1))^2+(coord(1,1,2)-coord(1,2,2))^2);
d(2) = sqrt((coord(1,2,1)-coord(2,1,1))^2+(coord(1,2,2)-coord(2,1,2))^2);
d(3) = sqrt((coord(2,1,1)-coord(2,2,1))^2+(coord(2,1,2)-coord(2,2,2))^2);
d(4) = sqrt((coord(2,2,1)-coord(1,1,1))^2+(coord(2,2,2)-coord(1,1,2))^2);

d = floor(mean(d));

A = zeros(d,d);
A(1,1) = 100;
A(1,d) = 100;
A(d,1) = 100;
A(d,d) = 100;

cmos = A;

camera_phone = uint8(255.*points./max(points(:)));

%% Draw control points
selectedFixedPoints  = [1,1;d,1;1,d;d,d];
selectedMovingPoints = [squeeze(coord(1,2,:))';squeeze(coord(1,1,:))';squeeze(coord(2,1,:))';squeeze(coord(2,2,:))'];
[selectedMovingPoints,selectedFixedPoints] = cpselect(camera_phone,cmos,selectedMovingPoints,selectedFixedPoints,'Wait',true);

%% Find transformation

tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'affine');
%Jregistered = imwarp(camera_phone,tform,'OutputView',imref2d(size(cmos)));
%figure,imshowpair(cmos,Jregistered)
pixels = [d d];

date_save = char(datetime('now','TimeZone','local','Format','yyyyMMdd'));
save(['../FOVmicroscope_',date_save,'.mat'],'tform','pixels');
%%

J = imwarp(squeeze(Points),tform,'OutputView',imref2d(pixels));
J = fliplr(J);
figure(10),imagesc(J),axis image
end
