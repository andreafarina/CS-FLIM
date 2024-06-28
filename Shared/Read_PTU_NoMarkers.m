function [t,dtof] = Read_PTU_NoMarkers()  
%function Read_PTU_AG % Read PicoQuant Unified TTTR Files
% This is demo code. Use at your own risk. No warranties.
% Marcus Sackrow, PicoQuant GmbH, December 2013
% Peter Kapusta, PicoQuant GmbH, November 2016
% Edited script: text output formatting changed by KAP.

% Edited by Alberto Ghezzi (Politecnico di Milano) for increasing
% conversion speed. 2023

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

% start Main program
[filename, pathname]=uigetfile('*.ptu', 'T-Mode data:');
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
    
    if strcmp(EvalName,'MeasDesc_Resolution')
        dt = TagFloat;
    end
     
    if strcmp(EvalName,'MeasDesc_GlobalResolution')
        T = TagFloat;
    end

    if strcmp(TagIdent, 'Header_End')
        break
    end
end
fprintf(1, '\n----------------------\n');
%% READ DAND DECODE DATA
endOfHeader = ftell(fid);

T3WRAPAROUND = 1024;% % Constant provided by PicoQuan
OverflowCorrection = 0;
chunkSize = 1e7;

fseek(fid,0, 'eof');
sizefile = (ftell(fid)-endOfHeader)/4;
fseek(fid,endOfHeader, 'bof');
bar = waitbar(0,'Reading data');
dtof = 0; i = 1;
while 1   
    T3Record = fread(fid, chunkSize, 'ubit32');
    if size(T3Record,1) == 0
        close(bar);
        break;
    end
    waitbar(round(i*chunkSize/sizefile,2),bar,'Reading data');
    nsync = bitand(T3Record,1023);                  % the lowest 10 bits:
    dtime = bitand(bitshift(T3Record,-10),32767);   % the next 15 bits:
    channel = bitand(bitshift(T3Record,-25),63);    % the next 6 bits:
    %special = bitand(bitshift(T3Record,-31),1);    % the last bit:
    
    dtime(channel==1) = 0;
    
    true_nSync = cumsum((channel==63).*(T3WRAPAROUND.*nsync)); %Conversion with overflow correction. Channel 63 means overflow
    OverflowCorrection = true_nSync(end) + OverflowCorrection;
    timeTag = (nsync + true_nSync)';
    dataset = [(1:length(nsync))' channel (nsync + true_nSync) 1e9*timeTag'./TTResult_SyncRate dtime];
    dataset = dataset(not(dataset(:,2)==63),[2 4 5]);
    dataset(dataset(:,1)==63,:) = [];  
    hist = histcounts(dataset(:,3),0:2^16);
    dtof = hist+dtof;
    i = i+1;
end
if exist([pathname, '/Converted/'],'file') == 0
    mkdir([pathname, '/Converted/'])
end
t = dt:dt:(T+100*dt);
save([pathname, '/Converted/',filename(1:end-4),'_refolded.mat'],'dtof','t','-v7.3');

end
