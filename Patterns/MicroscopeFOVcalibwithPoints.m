
figure,imagesc(Points),axis image, colorbar, title('Draw clockwise')
points = Points;
h = drawpolygon;
%h = drawpolygon('Position',[61.1240455950077,247.836343913881;269.083012836479,51.4520774895263;464.256249473950,259.236469994108;256.172837268903,454.278225910346]);
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
camera_phone = zeros(512,512,3);
%points = imbinarize(points,10000);
camera_phone(:,:,1) = points;
camera_phone(:,:,1) = round(255.*camera_phone(:,:,1));

%% Draw control points
selectedFixedPoints  = [1,1;d,1;1,d;d,d];
selectedMovingPoints = [squeeze(coord(1,2,:))';squeeze(coord(1,1,:))';squeeze(coord(2,1,:))';squeeze(coord(2,2,:))'];
[selectedMovingPoints,selectedFixedPoints] = cpselect(camera_phone,cmos,selectedMovingPoints,selectedFixedPoints,'Wait',true);

%% Find transformation
%camera_phone = sum(double(camera_phone),3);
%cmos = sum(double(cmos),3);

tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'affine');
Jregistered = imwarp(camera_phone,tform,'OutputView',imref2d(size(cmos)));
figure,imshowpair(cmos,Jregistered)
pixels = [d d];
save('520_microscope20x_20220614.mat','tform','pixels');
%%

J = imwarp(squeeze(Points),tform,'OutputView',imref2d(pixels));
J = fliplr(J);
figure(10),imagesc(J),axis image


% Plot FOV
normalize = @(X)X./max(X); 
figure(10),line([round(size(J,1)/2),round(size(J,1)/2)],[1,d], 'Color', 'b');
line([1,d], [round(size(J,1)/2),round(size(J,1)/2)],'Color', 'r');
line([1,d], [1,d],'Color', 'y')
figure(11),plot(normalize(J(:,round(size(J,1)/2)))),hold on
plot(normalize(J(round(size(J,1)/2),:)))
I = eye(d,d);
A = diag(J.*I);
plot(normalize(A))
grid minor