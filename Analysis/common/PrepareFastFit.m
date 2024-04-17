function [app] = PrepareFastFit(app)

            if isfield(app,'visplot')==1
               visplot = app.visplot;
            else    
               visplot = 1;
            end
            
            spc_raw = app.spc;
            app.spc = app.spc(app.t_start:app.t_stop,app.l_start:app.l_stop,:);
           
            if not(strcmp(app.BasisDropDown.Value,'Raster Scan'))
            switch app.PatterntypeDropDown_2.Value
                case 'Pos-Neg'
                    spc = CWHandling(app.spc,app.CWDropDown.Value,'Pos-Neg');
                    had = PosNeg2Had(spc);
                case 'Shifted Simulation'
                    spc = CWHandling(app.spc,app.CWDropDown.Value,'Pos-Neg'); %for CW handl. Pos simulation implies a pos-neg meas.
                    had = PosNeg2Pos(spc,app.NegEditField.Value,app.MethodDropDown.Value);
                case 'Shifted'
                    K = sqrt(2^floor(log2(size(app.spc,3)-app.NegEditField.Value)));
                    had = Shifted2Had(app.spc(:,:,1:(K^2+app.NegEditField.Value)),K,app.NegEditField.Value,app.MethodDropDown.Value,app.CWDropDown.Value);
                case 'Pos'
                    if mod(log2(sqrt(size(app.spc,3))),1)==0
                       % 1024 measurements, with CW
                       had = app.spc;
                    elseif mod(log2(sqrt(size(app.spc,3)/2)),1)==0
                       % Original posneg measurement
                       spc = CWHandling(app.spc,app.CWDropDown.Value,'Pos-Neg');
                       had = spc(:,:,1:2:end);
                    else
                       K = sqrt(2^floor(log2(size(app.spc,3)-app.NegEditField.Value)));
                       if app.NegEditField.Value == 1
                       CW = zeros(size(app.spc));
                       [~,notCW] = min(squeeze(sum(app.spc(:,:,1:(K^2+app.NegEditField.Value)),1:2)));
                       CW(:,:,notCW) = app.spc(:,:,1)+app.spc(:,:,2);
                       spc = CWHandling(app.spc(:,:,1:(K^2+app.NegEditField.Value)),app.CWDropDown.Value,'Shifted',CW);    
                       had = spc(:,:,[1 3:end]);
                       else
                       [~,had] = Shifted2Had(app.spc(:,:,1:(K^2+app.NegEditField.Value)),K,app.NegEditField.Value,app.MethodDropDown.Value,app.CWDropDown.Value);
                       end
                     end
            end
            else
                had = [];
            end
            app.had_spectrum.YData = squeeze(sum(app.spc,1:2));
            
            if strcmp(app.CompressivesensingSwitch.Value,'On') && strcmp(app.CompressionTypeDropDown.Value,'LSQR')
                had = PutCofficients2Zero(had,app.RatioEditField.Value);    
            end
            
            app.had = had;
            
            if size(app.loaded,1) == 0 && not(strcmp(app.BasisDropDown.Value,'Raster Scan'))
                try
                    load(['../../Calibrations/FLIM_',app.BasisDropDown.Value,'_',num2str(size(had,3))],'M','Pc','Pr')
                catch
                    [file, path] = uigetfile(['FLIM_',app.BasisDropDown.Value,'_',num2str(size(had,3))]);
                    figure(app.UIFigure)
                    if file == 0
                        return;
                    end
                    load([path file],'M');
                    app.loaded = 1;
                    app.M = M;
                end
            elseif size(app.loaded,1) == 1 && not(strcmp(app.BasisDropDown.Value,'Raster Scan'))
                    M = app.M;
            end
            
            if not(strcmp(app.BasisDropDown.Value,'Raster Scan')) && not(strcmp(app.PatterntypeDropDown_2.Value,'Pos'))
                M = M(1:2:end,:) - M(2:2:end,:);
            elseif strcmp(app.PatterntypeDropDown_2.Value,'Pos')
                M = M(1:2:end,:);
            end
            
            run('./TVAL/warm_up.m')

            t = app.t; lambda = app.lambda(app.l_start:app.l_stop); spc = app.spc; maxTWinLength = app.t_stop; ratio = app.RatioEditField.Value/100;
            %% Global Fit
            [~,a] = max(squeeze(sum(spc,1:2)));
            [~,b] = max(sum(spc,2:3));
            t = t(b:end)-t(b);
            [~,d] = min(abs(t-maxTWinLength));
            t = t(1:d);
            STmap = spc(b:(d+b-1),:,a)';

            %Number of components based on SVD
            [~,B,~] = svd(STmap);
            [~,D,~] = svd(real(sqrt(STmap)));
            D = diag(D);
            N = sum(diag(B(1:(size(B,1)-1),1:(size(B,1)-1))+diag(B(2:size(B,1),2:size(B,1))))>=sqrt(2)*D(1));                        
            %N = sum(diag(B)>1);%
            figure,loglog(diag(B),'-o'),hold on,line([1 size(spc,2)],0.8.*[D(1) D(1)],'color','r'),grid minor
            ylabel('Amplitude [a.u.]'), title('SVD'),xlabel('Value')
            prompt = {'Number of coefficients:'};
            dlgtitle = 'Input';
            dims = [1 35];
            definput = {num2str(N)};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            if size(answer,1) == 0
              answer = definput;
            end
            t0max = 4;
            N = str2double(answer);
            tau_0 = rand(N,1)+(((1:N)-1)*t0max/N)';
            NON_NEG = 0; showplots = 1;
            [tau,DAS,T,spcFit] = GlobalFit(t,lambda,STmap,tau_0',NON_NEG,showplots);

            had = had(b:(d+b-1),:,:);

            K = sqrt(size(had,3));
            had = had(:,:,1:round(K^2*(1-ratio))); % da poco dopo il picco

            Cmap = FastFit(had,DAS./sum(DAS,2),T);
            %Cmap = FastFit(had,DAS,T);
            %% TVAL Reconstruction
            opts.mu = eval(app.TVALmu.Value); %default value, suggested by authors
            opts.beta = eval(app.TVALbeta.Value); %default value, suggested by authors
            opts.maxit = 3000;
            opts.tol = 1E-8;
            opts.isreal = true;
            
            [~,e] = max(squeeze(sum(had,1:2)));
            [Cmap1, M1] = MoveCW2End(Cmap,M,e);
            m = round((1-ratio)*K^2);
            
            tic
            for i = 1:length(tau)
                Cimage(i,:,:) = (TVAL3(M1(1:m,:),Cmap1(1:m,i),K,K,opts));
                %Cimage(i,:,:) = SHReconstruction(permute(squeeze(Cmap1(1:m,i)),[2 3 1]),M1(1:m,:),0);
            end
            toc
            assignin('base', 'tval_opts', opts);
            assignin('base', 'ratio', ratio);
            %% Display results
            CreateFLIMMaps(Cimage,tau,DAS,lambda,N);
            app.im(1,1,:,:) = squeeze(sum(sum(T)'.*Cimage,1));
            assignin('base', 'DAS', DAS);
            assignin('base', 'tau', tau);
           
end
%% AUXILIARY FUNCTIONS  %%
function [lifetime,concentrations] = CreateFLIMMaps(Cimage,tau,DAS,lambda,N)
for i = 1:length(tau)
  CimagePerc(i,:,:) = squeeze(Cimage(i,:,:));
  hfig = figure;
  wl = abs(DAS(i,:))*lambda'/sum(abs(DAS(i,:)));
  S = Hyperspectral2RGB(wl,permute(CimagePerc(i,:,:),[4 1 2 3]));
  imagesc(S),title(['\tau = ',num2str(round(tau(i),1)),' ns']),axis image, axis off
  cbh = colorbar;
  get(get(cbh, 'Children'))
  [r,g,b] = Wavelength2color(wl);
  colorMap = [linspace(0,r,256)', linspace(0,g,256)',linspace(0,b,256)'];
  colormap(colorMap);
  II = CimagePerc(i,:,:);
  set(cbh,'XTickLabel',num2cell(round(100.*linspace(0,max(II(:)),length(cbh.Ticks)))./100))
end
% DRAW TAU AVE LIFETIME MAP
SIM = 0; %flag for simulation
if not(SIM==1)
Cimage(Cimage<0) = 0;
concentrations = Cimage;
Csum = squeeze(sum(concentrations,1));
mask = find(Csum>mean(Csum(:)));
mask = reshape(ismember(1:length(Csum(:)),mask),[size(Csum,1) size(Csum,2)]);
else
    Cimage(Cimage<0) = 0;
    concentrations = Cimage;
    Csum = squeeze(sum(concentrations,1));
    mask = find(Csum>0.083.*mean(Csum(:)));
    mask = reshape(ismember(1:length(Csum(:)),mask),[size(Csum,1) size(Csum,2)]);
end
A = repmat(permute(mask,[3 1 2]),[N 1 1]).*concentrations./repmat(permute(Csum,[3 1 2]),[N 1 1]);
lifetime(1,:,:) = reshape(tau*reshape(A,N,[]),[1 size(Csum,1) size(Csum,2)]);
concentrations = permute(Csum.*mask,[3 1 2]);
lifetime = lifetime(1,:,:);
FLIMSetColorScaleAndPlot(lifetime,concentrations);
assignin('base', 'lifetime', lifetime);
assignin('base', 'concentrations', concentrations);
end

