function [im_fused] = data_fusion(Data)
% Data fusion codes for the manuscript:
% Soldevila et al. Giga-voxel multidimensional fluorescence imaging combining single-pixel detection and data fusion
% published in OSA Optics Letters (2021)
%
% If the code turns to be useful for you, please consider citing it as:
% Fernando Soldevila, Armin Lenz, Alberto Ghezzi, Andrea Farina, Cosimo D'Andrea, and Enrique Tajahuerce, 
% “Giga-voxel multidimensional fluorescence imaging combining single-pixel detection and data fusion”, 
% Opt. Lett. Doc. ID 434127, Posted 30 Jul 2021 (In press), DOI: 10.1364/OL.434127
%% Pathing
%addpath('.\routines');

%% Cast datasets in 4D lambda-t-x-y format
CCD = permute(Data.CCD,[3 4 1 2]);
PMT = permute(Data.PMT,[4 3 1 2]);
SPEC = permute(Data.L16,[3 4 1 2]);
% Resolutions/channels:
Res.spatLow = size(PMT,3);
Res.spatHigh = size(CCD,3);
Res.tempLow = 1;
Res.tempHigh = size(PMT,2);
Res.specLow = 1;
Res.specHigh = size(SPEC,1);

%% Base functions definitions (used to build objective function and its gradient)

bf.S = @(x)specInt(x);                          %Integrate all lambdas
bf.St = @(x)specDeInt(x, Res.specHigh);         %Generate new spectral channels replicating available one
bf.T = @(x)timeInt(x);                          %Integrate all times
bf.Tt = @(x)timeDeInt(x, Res.tempHigh);         %Generate new time channels replicating available one
bf.K = @(x) spaceResample(x,Res.spatLow);       %Downsample in spatial domain
bf.Kt = @(x) spaceResample(x,Res.spatHigh);     %Upsample in spatial domain

%% Create initial estimation
%%%%%Random initial estimation%%%%%
init = ones(Res.specHigh,Res.tempHigh,Res.spatHigh,Res.spatHigh)...
    + 0.1*randn(Res.specHigh,Res.tempHigh,Res.spatHigh,Res.spatHigh);
Xnew = init/norm(init(:)); %Normalize proposed solution
clear aux_pmt aux_spec init

%%Show initialization
f1 = figure(1);
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
subplot(1,2,1)
imshow(squeeze(sum(sum(Xnew,1),2)),[],'Colormap',parula); axis square; colorbar; title('Intensity image');
subplot(1,2,2)
imshow(squeeze(sum(sum(permute(Xnew,[2 1 3 4]),3),4)),[],'YData',Data.time,'XData',Data.lambda,'Colormap',parula);
axis square; axis on; colorbar; ylabel('Time [ns]'); xlabel('\lambda [nm]')
drawnow

%% Normalize meaasurements
PMT = PMT*norm(tens2vec(bf.K(bf.S(Xnew))))/norm(PMT(:));
SPEC = SPEC*norm(tens2vec(bf.K(bf.T(Xnew))))/norm(SPEC(:));
CCD = CCD*norm(tens2vec(bf.S(bf.T(Xnew))))/norm(CCD(:));

%% Prepare regularization, define parameters and functions
%%%%%Regularization parameters%%%%%
reg.beta = 1;            % Comparison with high-res temporal data term
reg.gamma = 1;           % Comparison with high-res spectral data term
reg.epsilon = 2;         % Comparison with high-res spatial data term

reg.iter = 200;          % Number of iterations (maximum)
reg.plotFreq = 1;        % Show results every plotFreq iterations
reg.initStepSize = 0.1;  % Initial stepsize 
reg.stepSize = [];       % Stepsize
reg.btParam = 0.5;       % Backtracking Parameter: should be between 0.1 and 0.8

%%%%%Define regularization terms%%%%%
%Comparison with high-res temporal data term
F.F1 = @(x) reg.beta*0.5*norm(tens2vec(bf.K(bf.S(x))-PMT))^2;
%Comparison with high-res spectral data term
F.F2 = @(x) reg.gamma*0.5*norm(tens2vec(bf.K(bf.T(x))-SPEC))^2;
%Comparison with high-res spatial data term
F.F3 = @(x) reg.epsilon*0.5*norm(tens2vec(bf.S(bf.T(x))-CCD))^2;

%Complete objective function as sum of all the terms
F.F = @(x) F.F1(x) + F.F2(x) + F.F3(x);% 

%%%%%Define gradient of objective function%%%%%
%Comparison with high-res temporal data term
dF.dF1 = @(x) reg.beta*bf.St(bf.Kt(bf.K(bf.S(x))-PMT));
%Comparison with high-res spectral data term
dF.dF2 = @(x) reg.gamma*bf.Tt(bf.Kt(bf.K(bf.T(x))-SPEC));
%Comparison with high-res spatial data term
dF.dF3 = @(x) reg.epsilon*bf.Tt(bf.St(bf.S(bf.T(x))-CCD));

%Complete gradient of objective function as sum of all the terms
dF.dF = @(x) dF.dF1(x) + dF.dF2(x) + dF.dF3(x);

%% Solve gradient descent

for k=1:reg.iter    
    Xold = Xnew;                            %update solution
    gradient = dF.dF(Xold);                 %calculate gradient
    gradient = gradient*(norm(Xold(:))/norm(gradient(:)));
    [reg.stepSize, breakCond] = backTrackingLineSearch(F.F,Xold,gradient,reg.initStepSize,reg.btParam); %calculate step size
    if breakCond
        disp('**************************************************************************************************************')
        disp('Iteration stopped. Backtracking doesn''t seem to find a good step size. Estimation may be already good enough.')
        break;
    end
    Xnew = Xold - reg.stepSize*gradient;    %descend to find next solution
    %Show gradient descent info
    if k == 1 || mod(k,reg.plotFreq) == 0
        disp('******************************************')
        fprintf('Iteration number %d \n',k)
        fprintf('Step size %d \n',reg.stepSize)
        fprintf('Objective function %d \n',F.F(Xnew))
    end
    %Show results
    if k == 1 || mod(k,reg.plotFreq) == 0
        f2 = figure(2);
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        subplot(1,2,1)
        imshow(squeeze(sum(sum(Xnew,1),2)),[],'Colormap',parula); axis square; axis on; colorbar; title('Intensity image');
        subplot(1,2,2)
        imshow(squeeze(sum(sum(permute(Xnew,[2 1 3 4]),3),4)),[],'YData',Data.time,'XData',Data.lambda,'Colormap',parula);
        axis square; axis on; colorbar; ylabel('Time [ns]'); xlabel('\lambda [nm]')
        drawnow
    end
end


im_fused = permute(Xnew,[2 1 3 4]);
close([f1 f2])
end