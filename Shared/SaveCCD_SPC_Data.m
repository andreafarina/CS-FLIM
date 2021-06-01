function SaveCCD_SPC_Data
% ================2048=====================================================
% Save raw Data into a .mat file
% Andrea Farina - CNR - 09/09/2015
% Alberto Ghezzi - Polimi - 26/05/2021
% =========================================================================
%% Select if CCD file or SPC file
CCD = 0;
SPC = 1;
MANTA = 0;
%% Path
data_folder = '\\nas-dfis.fisica.polimi.it\work\DATA\CS_FLIM\UJI2';
out_folder = 'C:\Users\ghezz\OneDrive - Politecnico di Milano\';
ccd_ext='SPE';
tr_ext='sdt'; %if manta -> sdtm

%% Measurement info
day = '20210430';
prefix = 'UJI_20210430';
N_Pattern_IN = 2049;%16*16*2+1;                          % Input pattern number                                          % Angle number
N_Pattern_OUT = 1;%2*8*8;%16*16*2+1;               % Output pattern number

CF = 0; % Select if the SPC has been acquired with Continuous Flow (1 = Yes)

% SPC parmaters
if SPC == 1
    N_lambda = 16;
    Num_chan = 16384/4;  % Manta 16384 Becker 4096
elseif MANTA == 1
    N_lambda = 32;       % 32 MANTA; 16 L 16; CCD 1
    Num_chan = 16384;    %Manta 16384 Becker 4096
elseif CCD == 1
    N_lambda = 1;
end

% CCD parameters
Npixel = 512;                                         % 512x512 pixwl sensor
B = 1;                                                % Hardware binning


%% ========================================================================
%% Construct path
data_path = [data_folder,filesep,day,filesep];
out_path = [out_folder,filesep,day,filesep];

%% otutup filename
if CCD == 1
    output_file = [prefix,'_CCD_raw'];
elseif SPC == 1
    output_file = [prefix,'_SPC_raw'];
elseif MANTA == 1
    output_file = [prefix,'_SPCm_raw'];
end
%% Data format
filepath = [data_path prefix];

%% CCD file
if CCD == 1
    camera = zeros(Npixel/B,Npixel/B,N_lambda,N_Pattern_IN,N_Pattern_OUT);
    tic;
    i = 1;
    for j = 1:N_Pattern_IN
        for k = 1:N_Pattern_OUT
            camera(:,:,i,j,k) = loadSPE([filepath '_' num2str(0) '_' num2str(j-1) '_' num2str(k-1) '.',ccd_ext],[Npixel/B Npixel/B]);
        end
    end
    toc;
    
end
%% SPC file
if SPC == 1
    spc = zeros(Num_chan,N_lambda,N_Pattern_IN,N_Pattern_OUT);
    tic;
    if CF == 0
        for j = 1:N_Pattern_IN
            for k = 1:N_Pattern_OUT
                all_data = f_read_sdt_01([filepath '_' num2str(0) '_' num2str(j-1) '_' num2str(k-1) '.',tr_ext]);%-im_b;        % carico la curva dello SPAD
                spc(:,:,j,k) = reshape(all_data,[Num_chan N_lambda]);
            end
        end
        spc = squeeze(spc);
    else
        data = ReadSDT_ContFlow(filepath,Num_chan,N_lambda,(N_Pattern_OUT-1)/(2^17/Num_chan),2^17/Num_chan);
        [~,b] = min(sum(data,1:2));
        spc(:,:,1) = data(:,:,b); % Replicate the "off pattern" measurement and put it at the beginning, to be used as a background
        spc(:,:,2:end) = data;
        spc = squeeze(spc);
    end
    toc;
    
    if length(size(spc)) == 3 %Case Pattern In AND Out to be developed
        data = sum(spc,1:2);
        figure,plot(data(:)), title('Acquired patterns (sum over t and \lambda)')
        xlim([-5 size(spc,3)])
        grid minor
    end
end

%% SPC file
if MANTA == 1
    spc=zeros(Num_chan,N_lambda,N_Pattern_OUT);
    tic;
    for k = 1:N_Pattern_OUT
        all_data = readBinaryMANTA([filepath '_' num2str(0) '_' num2str(0) '_' num2str(k-1) '.' tr_ext]);
        spc(:,:,k) = reshape(all_data,[Num_chan N_lambda]);
    end
    toc;
end

%% Save to data file
if exist(out_path,'file') == 0
    mkdir(out_path)
end
out_filepath = [out_path output_file];
if CCD == 1
    save(out_filepath,'camera','-v7.3');
elseif SPC == 1 || MANTA == 1
    save(out_filepath,'spc','-v7.3');
end
end
