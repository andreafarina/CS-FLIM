function [im_fused] = data_fusion(Data,param)

%% Pathing
path = mfilename('fullpath');
addpath([path(1:(end-length(mfilename()))) 'routines']);

%% Cast datasets in 4D lambda-t-x-y format
CCD = permute(Data.CCD,[3 4 1 2]);
PMT = Data.PMT;
PMTlambdaTime = sum(PMT,3:4);

if isfield(Data,'GT')
GT = Data.GT;
end
% Resolutions/channels:
Res.spatLow = size(PMT,3);
Res.spatHigh = size(CCD,3);
Res.tempLow = 1;
Res.tempHigh = size(PMT,2);
Res.specLow = 1;
Res.specHigh = size(PMT,1);

%% Base functions definitions (used to build objective function and its gradient)
global M Mt M_int_time M_int_time_t M_int_spec M_int_spec_t Msqueeze Msqueeze_t
M = SpaceResampleMatrix([Res.spatLow,Res.spatLow],[Res.spatHigh,Res.spatHigh]);
M = kron(M,speye(Res.tempHigh*Res.specHigh));
Mt = M';
Msqueeze = SpaceResampleMatrix([1,1],[Res.spatHigh,Res.spatHigh]);
Msqueeze = kron(Msqueeze,speye(Res.tempHigh*Res.specHigh));
Msqueeze_t = Msqueeze';

% % M_int_time = TimeIntMatrix([Res.specHigh,Res.tempHigh,Res.spatHigh,Res.spatHigh]);
% % M_int_time_t = M_int_time';
% % M_int_spec = SpecIntMatrix([Res.specHigh,1,Res.spatHigh,Res.spatHigh]);
% % M_int_spec_t = M_int_spec';

bf.S = @(x)specInt(x);                          %Integrate all lambdas
bf.St = @(x)specDeInt(x, Res.specHigh);         %Generate new spectral channels replicating available one
%  bf.S = @(x)reshape(M_int_spec*tens2vec(x),[1,1,Res.spatHigh,Res.spatHigh]); %Integrate all lambdas
% bf.St = @(x)reshape(M_int_spec_t*tens2vec(x),[Res.specHigh,1,Res.spatHigh,Res.spatHigh]);         %Generate new spectral channels replicating available one

bf.T = @(x)timeInt(x);                          %Integrate all times
bf.Tt = @(x)timeDeInt(x, Res.tempHigh);         %Generate new time channels replicating available one
% bf.T = @(x)reshape(M_int_time*tens2vec(x),[Res.specHigh,1,Res.spatHigh,Res.spatHigh]);
%  bf.Tt = @(x)reshape(M_int_time_t*tens2vec(x),[Res.specHigh,Res.tempHigh,Res.spatHigh,Res.spatHigh]);

%bf.K = @(x) spaceResample_AF(x,Res.spatLow);       %Downsample in spatial domain
%bf.Kt = @(x) spaceResample_AF(x,Res.spatHigh);     %Upsample in spatial domain
bf.K = @(x) reshape(Mt*x(:),size(PMT));       %Downsample in spatial domain
bf.Kt = @(x) reshape(M*x(:),[Res.specHigh,Res.tempHigh,Res.spatHigh,Res.spatHigh]);     %Upsample in spatial domain

%% Create initial estimation
%%%%%Random initial estimation%%%%%
for li = 1:size(PMT,1)
    for ti = 1:size(PMT,2)
        start(li,ti,:,:) = imresize(squeeze(PMT(li,ti,:,:)),[size(CCD,3) size(CCD,4)],'box');
    end
end
start = start.*Data.MaskCCD;
start = start./sum(start(:)).*sum(CCD(:));
Xnew = start;

% %%Show initialization
% f1 = figure(1);
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% subplot(1,2,1)
% imshow(squeeze(sum(sum(Xnew,1),2)),[],'Colormap',parula); axis square; colorbar; title('Intensity image');
% subplot(1,2,2)
% imshow(squeeze(sum(sum(permute(Xnew,[2 1 3 4]),3),4)),[],'YData',Data.time,'XData',Data.lambda,'Colormap',parula);
% axis square; axis on; colorbar; ylabel('Time [ns]'); xlabel('\lambda [nm]')
% drawnow

%% Prepare regularization, define parameters and functions
%%%%%Regularization parameters%%%%%
if nargin<2
reg.beta = 1;            % Comparison with high-res temporal data term
reg.epsilon = 1;         % Comparison with high-res spatial data term 
reg.gamma = 0; 
reg.delta = 0;             
else
reg.beta = param.beta;            
reg.epsilon = param.epsilon;         
reg.gamma = param.gamma; 
reg.delta = param.delta;             
end
reg.iter = 2000;          % Number of iterations (maximum)
reg.plotFreq = 1;        % Show results every plotFreq iterations
reg.initStepSize = 1e-2;  % Initial stepsize 
reg.stepSize = [];       % Stepsize
reg.btParam = 0.5;       % Backtracking Parameter: should be between 0.1 and 0.8

%%%%%Define regularization terms%%%%%
F.F1 = @(x) norm(tens2vec(bf.K(x)-PMT))^2;
%Comparison with high-res spatial data term
F.F2 = @(x) norm(tens2vec(bf.S(bf.T(x))-CCD))^2;

% Tikhonov regularization
F.F3 = @(x) norm(tens2vec(x))^2;%

% Global map lambda-time
F.F4 = @(x) norm(Msqueeze_t*x(:) - PMTlambdaTime(:))^2;

%Complete objective function as sum of all the terms
F.F = @(x) (reg.beta*0.5*F.F1(x) + reg.epsilon*0.5*F.F2(x) ...
    + reg.gamma*0.5*F.F3(x)) + reg.delta*0.5*F.F4(x); 

%%%%%Define gradient of objective function%%%%%
%Comparison with high-res temporal data term
dF.dF1 = @(x) bf.Kt(bf.K(x)-PMT);
%Comparison with high-res spatial data term
dF.dF2 = @(x) bf.Tt(bf.St(bf.S(bf.T(x))-CCD));
%dF.dF2 = @(x) reg.epsilon*bf.Tt(bf.T(x)-CCD);

dF.dF3 = @(x) (x);
dF.dF4 = @(x) reshape(Msqueeze*(Msqueeze_t*x(:) - PMTlambdaTime(:)),[Res.specHigh,Res.tempHigh,Res.spatHigh,Res.spatHigh]);

%Complete gradient of objective function as sum of all the terms
dF.dF = @(x) reg.beta*dF.dF1(x) + reg.epsilon*dF.dF2(x) + reg.gamma*dF.dF3(x) + reg.delta*dF.dF4(x);

%% Solve gradient descent
if isfield(Data,'GT')
fprintf('Objective function with GT %d \n',F.F(GT));   
end
iobj = 1;
data_fid(1) = (reg.beta>0) * F.F1(Xnew);
regul_CCD(1) = (reg.epsilon > 0)  * F.F2(Xnew);
regul_TK(1) = (reg.gamma > 0) * F.F3(Xnew);
regul_lambbdaT(1) = (reg.delta > 0) * F.F4(Xnew);

obj = F.F(Xnew);
objvec(1) = obj;
first_order_opt = 1e10;
tol = 1e-4;
tol = tol * obj;
k = 1;
obj_start = Inf;

t1 = cputime;

while (iobj < reg.iter) && (obj_start - obj > tol)
    Xold = Xnew;                            %update solution
    gradient = dF.dF(Xold).*Data.MaskCCD;   %calculate gradient
    gradient(isnan(gradient)) = 0;
    gradient = gradient*(norm(Xold(:))/norm(gradient(:)));
    r_new = gradient;
     first_order_opt = norm(r_new(:),inf);
     disp(['first order opt:', num2str(first_order_opt)]);
    if k == 1
         p_new = r_new;
    end
    if k > 1
        beta = norm(r_new(:))^2./norm(r_old(:))^2;
        p_new = r_new + beta*p_old;
    end
    % backtracking
    [reg.stepSize, breakCond,obj_start,obj] = backTrackingLineSearch(F.F,Xold,p_new,reg.initStepSize,reg.btParam); %calculate step size
    
    %if breakCond
       % break;
    if (obj_start - obj) < tol
         disp('reset CG');
         p_new = r_new;
         [reg.stepSize, breakCond,obj_start,obj] = backTrackingLineSearch(F.F,Xold,p_new,reg.initStepSize,reg.btParam); %calculate step size
     end
     % toast linesearch bracketng
     % [reg.stepSize, obj] = toastLineSearch (Xold, p_new, 10, F.F(Xold), F.F);
       
     Xnew = Xold - reg.stepSize*p_new;    %descend to find next solution
     
    %obj = F.F(Xnew);
    p_old = p_new;
    r_old = r_new;

    %Show gradient descent info
    if k == 1 || mod(k,reg.plotFreq) == 0
        disp('******************************************')
        fprintf('Iteration number %d \n',k)
        fprintf('Step size %d \n',reg.stepSize)
        fprintf('Objective function %d \n',obj)
        
    end
    %Show results
    if k == 1 || mod(k,reg.plotFreq) == 0
        f2 = figure(200);
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        subplot(3,2,1)
        imshow(squeeze(sum(sum(Xnew,1),2)),[],'Colormap',parula); axis square; axis on; colorbar; title('Intensity image');
        subplot(3,2,2)
        imshow(squeeze(CCD),[],'Colormap',parula); axis square; axis on; colorbar; title('CCD image');
        subplot(3,2,3)
        imshow(squeeze(sum(sum(permute(Xnew,[2 1 3 4]),3),4)),[],'YData',Data.time,'XData',Data.lambda,'Colormap',parula);
        axis square; axis on; colorbar; ylabel('Time [ns]'); xlabel('\lambda [nm]'),title('lambda - time map')
        Xnew_Tint = sum(sum(sum(Xnew,1),3),4);
        PMT_Tint = sum(sum(sum(PMT,1),3),4);
        Xnew_Sint = sum(Xnew,2:4);
        PMT_Sint = sum(PMT,[2,3,4]);
        subplot(3,2,4),semilogy([Xnew_Tint',PMT_Tint']),title('time')
        subplot(3,2,5),plot([Xnew_Sint(:),PMT_Sint(:)]),title('spectrum')
%         subplot(3,2,2)
%         imshow(squeeze(sum(sum(Xnew(1:2,20:80,:,:),1),2)),[],'Colormap',parula); axis square; axis on; colorbar; title('Intensity image');
%         subplot(3,2,3)
%         imshow(squeeze(sum(sum(Xnew(8:13,6:9,:,:),1),2)),[],'Colormap',parula); axis square; axis on; colorbar; title('Intensity image');
        
        iobj = iobj +1;
        objvec(iobj) = obj;
        data_fid(iobj) = (reg.beta>0) *F.F1(Xnew);
        regul_CCD(iobj) = (reg.epsilon>0) *F.F2(Xnew);
        regul_TK(iobj) = (reg.gamma>0) *F.F3(Xnew);
        regul_lambdaT(iobj) = (reg.delta>0) *F.F4(Xnew);

         subplot(3,2,6),semilogy([objvec',data_fid',regul_CCD',regul_TK',regul_lambdaT']),
         legend('obj','data fid','regul CCD','regul TK','regul lambdaT'),
         xlim([1,100])
%         
%         %figure(201),semilogy([objvec',data_fid',regul_CCD',regul_TK']),%legend('obj','data fid','regul CCD','regul TK')
%         
%         drawnow
    end
    k = k + 1;
end
t2 = cputime;
disp(['time: ',num2str(t2-t1),' s'])
Xnew = reshape(Xnew,size(start));
%disp(['Final Objective function: ',num2str(obj)])
disp(['F1 = ',num2str(F.F1(Xnew))]);
disp(['F2 = ',num2str(F.F2(Xnew))]);
disp(['F3 = ',num2str(F.F3(Xnew))]);
disp(['F4 = ',num2str(F.F4(Xnew))]);

if isfield(Data,'GT')
    err = norm(GT(:) - Xnew(:) );
    disp(['errorGT = ',num2str(err)]);
end
im_fused = permute(Xnew,[2 1 3 4]);
%close([f1 f2])

end