function [handle] = Read_PTU(tInt,R,N,filename,pathname)  
%function Read_PTU_AG % Read PicoQuant Unified TTTR Files
% This is demo code. Use at your own risk. No warranties.
% Marcus Sackrow, PicoQuant GmbH, December 2013
% Peter Kapusta, PicoQuant GmbH, November 2016
% Edited script: text output formatting changed by KAP.

% Note that marker events have a lower time resolution and may therefore appear
% in the file slightly out of order with respect to regular (photon) event records.
% This is by design. Markers are designed only for relatively coarse
% synchronization requirements such as image scanning.

% T Mode data are written to an output file [filename].out
% We do not keep it in memory because of the huge amout of memory
% this would take in case of large files. Of course you can change this,
% e.g. if your files are not too big.
% Otherwise it is best process the data on the fly and keep only the results.

% All HeaderData are introduced as Variable to Matlab and can directly be
% used for further analysis


% some constants
tyEmpty8      = hex2dec('FFFF0008');
tyBool8       = hex2dec('00000008');
tyInt8        = hex2dec('10000008');
tyBitSet64    = hex2dec('11000008');
tyColor8      = hex2dec('12000008');
tyFloat8      = hex2dec('20000008');
tyTDateTime   = hex2dec('21000008');
tyFloat8Array = hex2dec('2001FFFF');
tyAnsiString  = hex2dec('4001FFFF');
tyWideString  = hex2dec('4002FFFF');
tyBinaryBlob  = hex2dec('FFFFFFFF');
% RecordTypes
rtMultiHarpT3    = hex2dec('00010307');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $07 (MultiHarp)
rtMultiHarpT2    = hex2dec('00010207');% (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $07 (MultiHarp)

TTResultFormat_TTTRRecType = 0;
TTResult_NumberOfRecords = 0;
MeasDesc_Resolution = 0;
MeasDesc_GlobalResolution = 0;

% start Main program
%[filename, pathname]=uigetfile('*.ptu', 'T-Mode data:');
fid=fopen([pathname filename]);

fprintf(1,'\n');
Magic = fread(fid, 8, '*char');
if not(strcmp(Magic(Magic~=0)','PQTTTR'))
    error('Magic invalid, this is not an PTU file.');
end
Version = fread(fid, 8, '*char');
fprintf(1,'Tag Version: %s\n', Version);

%% READ HEADER
while 1
    % read Tag Head
    TagIdent = fread(fid, 32, '*char'); % TagHead.Ident
    TagIdent = (TagIdent(TagIdent ~= 0))'; % remove #0 and more more readable
    TagIdx = fread(fid, 1, 'int32');    % TagHead.Idx
    TagTyp = fread(fid, 1, 'uint32');   % TagHead.Typ
    % TagHead.Value will be read in the
    % right type function
    TagIdent = genvarname(TagIdent);    % remove all illegal characters
    if TagIdx > -1
        EvalName = [TagIdent '(' int2str(TagIdx + 1) ')'];
    else
        EvalName = TagIdent;
    end
    fprintf(1,'\n   %-40s', EvalName);
    % check Typ of Header
    switch TagTyp
        case tyEmpty8
            fread(fid, 1, 'int64');
            fprintf(1,'<Empty>');
        case tyBool8
            TagInt = fread(fid, 1, 'int64');
            if TagInt==0
                fprintf(1,'FALSE');
                eval([EvalName '=false;']);
            else
                fprintf(1,'TRUE');
                eval([EvalName '=true;']);
            end
        case tyInt8
            TagInt = fread(fid, 1, 'int64');
            fprintf(1,'%d', TagInt);
            eval([EvalName '=TagInt;']);
        case tyBitSet64
            TagInt = fread(fid, 1, 'int64');
            fprintf(1,'%X', TagInt);
            eval([EvalName '=TagInt;']);
        case tyColor8
            TagInt = fread(fid, 1, 'int64');
            fprintf(1,'%X', TagInt);
            eval([EvalName '=TagInt;']);
        case tyFloat8
            TagFloat = fread(fid, 1, 'double');
            fprintf(1, '%e', TagFloat);
            eval([EvalName '=TagFloat;']);
        case tyFloat8Array
            TagInt = fread(fid, 1, 'int64');
            fprintf(1,'<Float array with %d Entries>', TagInt / 8);
            fseek(fid, TagInt, 'cof');
        case tyTDateTime
            TagFloat = fread(fid, 1, 'double');
            fprintf(1, '%s', datestr(datenum(1899,12,30)+TagFloat)); % display as Matlab Date String
            eval([EvalName '=datenum(1899,12,30)+TagFloat;']); % but keep in memory as Matlab Date Number
        case tyAnsiString
            TagInt = fread(fid, 1, 'int64');
            TagString = fread(fid, TagInt, '*char');
            TagString = (TagString(TagString ~= 0))';
            fprintf(1, '%s', TagString);
            if TagIdx > -1
                EvalName = [TagIdent '{' int2str(TagIdx + 1) '}'];
            end
            eval([EvalName '=[TagString];']);
        case tyWideString
            % Matlab does not support Widestrings at all, just read and
            % remove the 0's (up to current (2012))
            TagInt = fread(fid, 1, 'int64');
            TagString = fread(fid, TagInt, '*char');
            TagString = (TagString(TagString ~= 0))';
            fprintf(1, '%s', TagString);
            if TagIdx > -1
                EvalName = [TagIdent '{' int2str(TagIdx + 1) '}'];
            end
            eval([EvalName '=[TagString];']);
        case tyBinaryBlob
            TagInt = fread(fid, 1, 'int64');
            fprintf(1,'<Binary Blob with %d Bytes>', TagInt);
            fseek(fid, TagInt, 'cof');
        otherwise
            error('Illegal Type identifier found! Broken file?');
    end
    if strcmp(TagIdent, 'Header_End')
        break
    end
end
fprintf(1, '\n----------------------\n');
%% READ AND DECODE DATA
endOfHeader = ftell(fid);
% % Constant provided by PicoQuant:
T3WRAPAROUND = 1024;
OverflowCorrection = 0; OverflowCorrectionReminder = 0;
% Init
b = 0;
fseek(fid,0, 'eof');
channel = zeros((ftell(fid)-endOfHeader)/4,1);
fseek(fid,endOfHeader, 'bof');
disp('Looking for markers...')
% Read 
while 1 %not(channel == 1)
        T3Record = fread(fid, 1e7, 'ubit32');     % all 32 bits:
        if size(T3Record,1)==0
            break;
        end
        channel((b+1):(b+length(T3Record)),1) = bitand(bitshift(T3Record,-25),63);    % the next 6 bits:
        special((b+1):(b+length(T3Record)),1) = bitand(bitshift(T3Record,-31),1);    
        b = b+length(T3Record);
end
posMarkers = find((channel==1)&(special==1));
posMarkers(end+1) = length(channel);
markIntervals = diff(posMarkers);
markIntervals(end) = markIntervals(end-1)*3;
%Rewind dataset and skip header
fseek(fid, endOfHeader+4, -1); %after read marker go back to previous position
disp('Markers found.')
bar = waitbar(0,'Reading data');
if length(markIntervals) < N*R
    disp('Missing markers. Measurement incomplete.')
    R = ceil((length(markIntervals)/N));
    %handle = [];
    %return;
end
T3Record = fread(fid, posMarkers(1)-2, 'ubit32'); % first void measurements: keep them for OVF
nsync = bitand(T3Record,1023);                  % the lowest 10 bits:
channel = bitand(bitshift(T3Record,-25),63);    % the next 6 bits:
OverflowCorrection = cumsum((channel==63).*(T3WRAPAROUND.*nsync))+OverflowCorrectionReminder; %Conversion with overflow correction. Channel 63 means overflow
OverflowCorrectionReminder = OverflowCorrection(end);
%spc = zeros(2520,1,N); measurement = zeros(R,2520,1,N+1);
k = 1; r = 1;
for i = 1:length(markIntervals)

    T3Record = fread(fid, markIntervals(i)+1, 'ubit32');
    fseek(fid, -4, 0);
    
    nsync = bitand(T3Record,1023);                  % the lowest 10 bits:
    dtime = bitand(bitshift(T3Record,-10),32767);   % the next 15 bits:
    channel = bitand(bitshift(T3Record,-25),63);    % the next 6 bits:
    special = bitand(bitshift(T3Record,-31),1);     % the last bit:
    %
    OverflowCorrection = cumsum((channel==63).*(T3WRAPAROUND.*nsync))+OverflowCorrectionReminder; %Conversion with overflow correction. Channel 63 means overflow
    timeTag = (nsync + OverflowCorrection)';
    rawData = [special channel (nsync + OverflowCorrection) 1e9*timeTag'./TTResult_SyncRate dtime];
    rawData = rawData(not(rawData(:,2)==63),[1 2 4 5]);   % remove overflows
    
    if i == length(markIntervals)
       rawData(end+1,:) = [1 1 rawData(end,2) 0];
    end
    spc(:,:,k) = CreateHistogram(tInt,rawData);
    
    if k == N
       A = zeros(2520,size(spc,2),N+1);
       A(:,:,2:(N+1)) = spc;
       measurement(r,:,:,:) = A;
       spc = zeros(2520,size(spc,2),N);
       k = 1; r = r + 1;
    else
        k = k + 1;
    end
    waitbar(round(i/length(markIntervals),2),bar,'Reading data');

end
close(bar);
if exist([pathname, 'Converted\'],'file') == 0
    mkdir([pathname, 'Converted\'])
end

save([pathname, 'Converted\',filename(1:end-4),'_folded.mat'],'measurement','-v7.3');
handle = [pathname, 'Converted\',filename(1:end-4),'_folded.mat'];
end

function [spc] = CreateHistogram(intTime,photons)
intTime = intTime*1e6;
markers = find((photons(:,1)==1)&(photons(:,2)==1));
N = (length(markers)-1);
spc = zeros(2^16,1,N);
for k = 1:N
    P = photons((markers(k)+1):(markers(k+1)-1),:);
    P(:,3) = P(:,3)-photons(markers(k),3);
    P(P(:,2)==63,:) = [];
    P(P(:,3)>intTime,:) = []; %remove counts arrived after the projection time
    for ch = 1:16
    counts = P(:,2)==(ch-1);
    if sum(counts) > 0
        spc(:,ch,k) = histcounts(P(counts,4),0:2^16);
    end    
    
    end
end
spc = spc(1:2520,:,:);
end