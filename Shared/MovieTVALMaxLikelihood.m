function [lifetime,concentrations] = MovieTVALMaxLikelihood(filenameSPC, filenameCMOS)
% Time Tagging TCPSC + SiPM CS-SPC Data Reconstruction
% Alberto Ghezzi, Politecnico di Milano

% Default dataset
if nargin == 0
    pathfile = '../../../../FLIM Setup/Data/20230622/';
    prefix = 'BEADS_11_folded';
    filenameSPC = [pathfile,prefix,'.mat'];
    filenameCMOS = '../../../../FLIM Setup/Data/20230622/cmos11';
end

addpath ../Analysis/
run('../Analysis/TVAL/warm_up.m')

% Load matrix used for generating patterns 
load('FLIM_Scrambled-Hadamard_1024.mat','M');
M = M(1:2:end,:);

% Measurement settings
intTime = 0.195; %ms, projection time per pattern 
Patt = 256; %number of used patterns
CR = 1-Patt/size(M,2); %compression ratio 
Npx = sqrt(size(M,2));
RepRate = 3.99e7; %laser repetition rate
DMDdeadTime = 0;
disp(filenameSPC);

%% Data Processing
dset = load(filenameSPC);
measurement = dset.measurement;
measurement = permute(measurement,[1 2 4 3]);
lifetime = zeros(size(measurement,1),Npx,Npx);
concentrations = lifetime;
nRep = size(measurement,1);

for i = 1:nRep

    %% Pre-Processing
    spc = permute(squeeze(measurement(i,:,:,:)),[1 3 2]);
    % Pile-up correction:
    for ii = 1:size(spc,3)
        spc(:,1,ii) = coates_correction_EA(squeeze(spc(:,1,ii)),2520,RepRate,intTime*1e-3);
    end
    % Bin and Noise BKG
    [t,spc] = BinData(spc,TimeVectorTCSPC(size(spc,1)),0.2);
    spc = spc(:,:,2:end)-mean(spc(4:10,:,2:end),1);
    % Pick decay from peak
    [~,b] = max(sum(spc,2:3));
    t = t(b:end)-t(b);

    %% Image Reconstruction    
    % Reconstruct flat pattern:
    spc(:,:,(size(spc,3)+1):(Patt/(1-CR)+1)) = zeros(size(spc,1),size(spc,2),(Patt/(1-CR)+1)-size(spc,3));
    notCW = 546;
    spc(:,:,notCW) = spc(:,:,1)+spc(:,:,2);
    had = spc(:,:,[1 3:end]);
    K = sqrt(size(had,3));
    
    % Cut empty temporal bins
    maxTWinLength = 9;
    [~,d] = min(abs(t-maxTWinLength));
    t = t(1:d);
    
    had = had(b:(d+b-1),:,:);
    m = round((1-CR)*K^2);

    % TVAL3 Settings:
    optsTV.mu = 2^8; 
    optsTV.beta = 2^5;

    optsTV.maxit = 3000;
    optsTV.tol = 1E-8;
    optsTV.isreal = true;
    
    % Starting condition:
    if exist("Cimage")
        Cimage_prev = Cimage;
    end
   

    % Tval reconstruction:
    for ti = 1:size(had,1)
        if exist("Cimage_prev")
            optsTV.init = squeeze(Cimage_prev(ti,:,:));
        end
        [Cmap, M1] = MoveCW2End(squeeze(had(ti,1,:)),M,notCW-1);
        Cimage(ti,:,:) = TVAL3(M1(1:m,:),Cmap(1:m,1),K,K,optsTV);
    end
% Maximum Likelihood Data fit
x0 = [1 1]; %[A, tau] Starting values
Ntau = 1;
opts = optimoptions('fminunc','TolFun',1e-6,'TolX',1e-6,...
        'Display','Off','StepTolerance',1e-6,...
        'OptimalityTolerance',1e-6);

for xi = 1:size(Cimage,2)
    for yi = 1:size(Cimage,3)
            ydata = squeeze(Cimage(:,xi,yi))';
            ydata = ydata./max(ydata);
            % Fit:
            x = fminunc(@objectivePOISS,x0,opts);
     
            c = (exp(-t./x(Ntau+1:end)'))'\squeeze(Cimage(:,xi,yi)); 
            plot(squeeze(Cimage(:,xi,yi))),hold on 
            plot(sum(c).*exp(-t./((c'*x(Ntau+1:end)')/sum(c)))),hold off
            drawnow
            lifetime(i,xi,yi) = (c'*x(Ntau+1:end)')/sum(c);
            concentrations(i,xi,yi) = sum(c);  
    end
end
raw_data(i,:,:,:) = Cimage;
end

%% Load camera file:
cmos = load(filenameCMOS);
camera = cmos.camera;

filename = [winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'), filesep,'Result.mp4'];
figure(100)
v = VideoWriter(filename,'MPEG-4');
v.FrameRate = 1/((intTime+DMDdeadTime)*1e-3*Patt);
open(v);
for i = 1:size(camera,1)
    ax(1) = subplot(1,2,1); 
    imagesc(rot90(squeeze(camera(i,:,:)),2)); axis image,colormap(ax(1),'gray'),caxis([0 1300]);
    cb = colorbar;
    cb.Label.String = 'Amplitude [a.u.]';
    ax(2) = subplot(1,2,2);
    DrawFLIMMap(lifetime(i,:,:),concentrations(i,:,:),1,[],mean(X)-3*std(X),mean(X)+3*std(X),ax(2));
    % if following a moving sample, uncomment this:
    %[x,y] = ginput(1);
    %xx(u) = round(x);
    %yy(u) = round(y); u = u +1;
    frame = getframe(gcf);
    writeVideo(v,frame);
    B = rot90(squeeze(camera(i,:,:)),2); 
    A = imresize(squeeze(concentrations(i,:,:)),size(B),'box');
    PSNR(i) = Compute_PSNR(B,A);
end
close(v)


disp(mean(PSNR(not(isnan(PSNR)))))
disp(std(PSNR(not(isnan(PSNR)))))

%% objective and models
function y = forward(x,~)
    y = x(1:Ntau) * (exp(-t./x(Ntau+1:end)'));
    %figure(1),semilogy(t,ydata,'.',t,y,'r-'),ylim([1e-6 1]),xlabel('time (ns)');
end

function phi = objectivePOISS(x) 
    yfwd = forward(x);
    phi = sum(yfwd-ydata.*log(yfwd));
    %figure(1),semilogy(t,ydata,'.',t,yfwd,'r-'),ylim([1e-6 1]),xlabel('time (ns)');
end

end

