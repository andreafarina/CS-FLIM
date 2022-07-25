function [] = ImageReconstruction(app)

            
            spc_raw = app.spc;
            app.spc = app.spc(app.t_start:app.t_stop,app.l_start:app.l_stop,:);
            
            if app.IntegralovertimeCheckBox.Value == 1
                app.spc = sum(app.spc,1);
            end
            
            if app.IntegraloverwavelengthCheckBox.Value == 1
                app.spc = sum(app.spc,2);
            end
            
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
            
            flag_flip = app.FlipCheckBox.Value;
            if strcmp(app.CompressivesensingSwitch.Value,'Off') || (strcmp(app.CompressivesensingSwitch.Value,'On')&&strcmp(app.CompressionTypeDropDown.Value,'LSQR'))
                switch app.BasisDropDown.Value
                    case 'Walsh-Hadamard'
                        im = WHReconstruction(had,flag_flip);
                    case 'Scrambled-Hadamard'
                        im = SHReconstruction(had,M,flag_flip);
                    case 'Raster Scan'
                        im = reshape(app.spc(:,:,1:app.XEditField_2.Value*app.XEditField.Value),[size(app.spc,1) size(app.spc,2) app.XEditField_2.Value app.XEditField.Value]);
                end
            else
                run('./TVAL/warm_up.m')
                ratio = app.RatioEditField.Value/100;
                K = size(app.had,3);
                m = round((1-ratio)*K);
                [im,opts] = InvertwTVAL3(app.had,m,M,app.camera,sqrt(K),eval(app.TVALmu.Value),eval(app.TVALbeta.Value));
                assignin('base', 'tval_opts', opts);
            end
            app.im = im;
            app.spc = spc_raw;

end