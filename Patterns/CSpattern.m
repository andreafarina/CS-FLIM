classdef CSpattern
    %CSPATTERN this class creates/allocate/save 2D
    % patterns
    %   Detailed explanation goes here
    
    properties
        type
        order
        Npatt
        kind        % posneg % shifted % shiftedCW
        cw
        dim
        stack       % image stacks
        Tmatrix     % matrix of true patterns
        Mmatrix     % measurement matrix
        seed = 122; % random seed
        % scrambledH dedicated
        Pr
        Pc
        % shiftedCW dedicated
        indexCW
    end
    
    methods
        function obj = CSpattern(type,order,seed)
            %init Construct an instance of this class
            %   create M
            if nargin > 2
                obj.seed = seed;
                rnd(see);
            end
            switch lower(type)
                case 'hadamard'
                    obj.dim = [order,order];
                    obj.Mmatrix = hadamard(order.^2);
                    obj.Tmatrix = obj.Mmatrix;
                    obj.stack = reshape(obj.Tmatrix,[order,order,order^2]);
                    obj.order = order;
                    obj.cw = 1;
                    
                case 'wh'
                    for i = 1:order^2
                        d = zeros(order,order);
                        d(i) = 1;
                        hk = fwht2D(d,order);
                        obj.stack(:,:,i) = hk;
                        obj.Mmatrix(i,:) = hk(:);
                        obj.Tmatrix = obj.Mmatrix;
                    end
                    obj.dim = [order,order];
                    obj.order = order;
                    obj.cw = 1;
                    
                case 'scrambledh'
                    obj.dim = [order,order];
                    K = order^2;
                    H = hadamard(K);
                    I = eye(K);
                    obj.Pr = I(randperm(K),:);
                    obj.Pc = I(randperm(K),:);
                    obj.Mmatrix = obj.Pr * H * obj.Pc;
                    obj.Tmatrix = obj.Mmatrix;
                    obj.stack = reshape(obj.Tmatrix,[order,order,order^2]);
                    obj.order = order;
                    obj.cw = 1;
                    
                case 'bernoulli'
                    obj.dim = [order,order];
                    obj.Mmatrix = random('bino',1,1/2,[order^2,order^2])*2 - 1;
                    obj.Tmatrix = obj.Mmatrix;
                    obj.stack = reshape(obj.Mmatrix,[order,order,order^2]);
                    obj.order = order;
                    obj.cw = 0;
                    
                case 'fourier'
                    % to be implemented
                otherwise
                    disp('Pattern not found');
                    return;
            end
            obj.Npatt = size(obj.Mmatrix,1);
            obj.type = lower(type);
        end
        
        function obj = restore(obj)
            obj.Mmatrix = obj.Tmatrix;
            obj.stack = reshape(obj.Tmatrix,[obj.order,obj.order,obj.order^2]);
            obj.dim = [obj.order,obj.order];
            obj.Npatt = size(obj.Tmatrix,1);
        end
        
        function obj = resample(obj,dim)
            rstack = zeros(dim(1),dim(2),obj.Npatt);
            for i = 1:obj.Npatt
                rstack(:,:,i) = imresize(obj.stack(:,:,i),dim,'box');
            end
            obj.stack = rstack;
            %obj.Mmatrix = (reshape(rstack,[prod(dim),obj.Npatt]))';Not needed to resample meas-matrix
            obj.dim = dim;
        end
        
        function obj = ConvertToPosNeg(obj)
            if obj.CheckKind
                return;
            end
            nrow = obj.Npatt;
            Mpos = obj.Mmatrix;
            Mpos(Mpos<0) = 0;
            Mneg = -obj.Mmatrix;
            Mneg(Mneg<0) = 0;
            M = [Mpos;Mneg];
            % permutation matrix
            Prr = zeros(2*nrow);
            Prr(1:nrow,1:2:end) = eye(nrow);
            Prr((1:nrow)+nrow,2:2:end) = eye(nrow);
            obj.Mmatrix = Prr' * M;
            obj.stack = reshape(obj.Mmatrix',[obj.dim(1),obj.dim(2),2*obj.Npatt]);
            obj.Npatt = 2 * obj.Npatt;
            obj.kind = 'posneg';
        end
        
        function obj = ConvertToShifted(obj)
            if obj.CheckKind
                return;
            end
            obj.Mmatrix = obj.Mmatrix + abs(min(obj.Mmatrix,[],2));
            obj.Mmatrix = obj.Mmatrix./max(obj.Mmatrix,[],2);
            obj.stack = reshape(obj.Mmatrix',[obj.dim(1),obj.dim(2),obj.Npatt]);
            obj.kind = 'shifted';
        end
        
        function obj = padToFitDim(obj,dim2)
            dd = round((dim2-obj.dim)/2);
            rstack = padarray(obj.stack,[dd(1),dd(2),0],'both');
            obj.dim = dim2;
            obj.stack = rstack;
        end
            
        function obj = TransposeStack(obj)
            obj.stack = permute(obj.stack,[2,1,3]);
        end
        
        function obj = DeleteCW(obj)
            if obj.cw == 0
                disp('this pattern does not have a CW');
            end
            a = sum(obj.Tmatrix,2);
            [~,imax] = max(a);
            obj.Mmatrix(imax,:) = [];
            obj.stack(:,:,imax) = [];
            obj.Npatt = obj.Npatt - 1;
            obj.cw = 0;
        end
        
        function obj = ConvertToCWShifted(obj,num)
            %ADDPOSNEGCW add posneg patterns evenly separated in the stack
            % only for hadamard, wh, scrambled
            if strcmpi(obj.type,'fourier')
                disp('This does not work for Fourier patterns');
            end
            if obj.cw == 1
                disp('Please Delete the CW!');
                return;
            end
            if ~strcmpi(obj.kind,'posneg')
                disp('Please convert to posneg!');
                return;
            end
            obj.indexCW = zeros(1,obj.Npatt);
            obj.indexCW(1:2:end) = 1;
            if num == 1
                obj.indexCW(2) = 1;
            else
                [~,d] = find((obj.indexCW == 0));
                b = round(linspace(1,numel(d),num));
                obj.indexCW(d(b)) = 1;
            end
            obj.Npatt = sum(obj.indexCW);
            obj.Mmatrix(~obj.indexCW,:) = [];
            obj.stack(:,:,~obj.indexCW) = [];
            obj.kind = 'shiftedCW';
        end
        
        function obj = AddZeroPattern(obj)
            d = zeros(1,size(obj.Mmatrix,2));
            obj.Mmatrix = [d;obj.Mmatrix];
            obj.Npatt = obj.Npatt + 1;
            obj.stack = reshape(obj.Mmatrix',[obj.dim(1),obj.dim(2),obj.Npatt]);
        end
        
        function [] = ShowPatterns(obj,space)
            %SHOWPATTERNS show patterns as mosaic
            %   Detailed explanation goes here
            if nargin<2
                space = 1;
            end
            Nx = size(obj.stack,1);
            Ny = size(obj.stack,2);
            stride = Nx + space;
            A = mean(obj.stack(:))*ones(stride*ceil(sqrt(obj.Npatt)),stride*ceil(sqrt(obj.Npatt)));
            figure;%ax = gca;
            for ip = 1:obj.Npatt
                d = zeros(stride,stride);
                d(space + 1:stride,1:Ny) = obj.stack(:,:,ip);
                id = zeros(ceil(sqrt(obj.Npatt)));
                id(ip) = 1;
                A = A + kron(id',d);
                imagesc(A)
            end
        end
        
        function [] = ShowMmatrix(obj)
            %SHOWMMATRIX show the measurement matrix
            %   Detailed explanation goes here
            figure,imagesc(obj.Mmatrix);
        end
        
        function [] = SaveStack(obj,bitdepth,fmt)
            if sum(obj.Mmatrix(:)<0) > 0
                disp('Convert to positive before save!');
                return;
            end
            switch bitdepth
                case 1
                    s = logical(obj.stack);
                case 8
                    s = im2uint8(obj.stack);
                case 16
                    s = im2uint16(obj.stack);
            end
            selpath = uigetdir(pwd);
            prefix = [selpath,filesep,obj.type,'_',num2str(obj.order)];
            switch lower(fmt(1:3))
                case 'png'
                    for i = 1:obj.Npatt
                        imwrite(s(:,:,i),[prefix,'_',num2str(i-1),'.',fmt]);
                    end
                case 'tif'
                    filename = [prefix,'.tif'];
                    imwrite(s(:,:,1),filename);
                    for i = 2:obj.Npatt
                        imwrite(s(:,:,i),filename,'WriteMode','append');
                    end
                case 'mat'
                    pattern = obj.stack;
                    M = obj.Mmatrix;
                    save(prefix,'pattern','M');
                case 'hdf'
                    filename = [prefix,'.hdf5'];
                    dataset = '/Sequence';
                    %s = permute(s,[3,2,1]);
                    H5WriteStack(filename,dataset,s,['uint',num2str(bitdepth)]);
                    h5writeatt(filename,dataset,'type',obj.type);
                    h5writeatt(filename,dataset,'order',uint8(obj.order));
                    h5writeatt(filename,dataset,'kind',obj.kind);
                    h5writeatt(filename,dataset,'cw',uint8(obj.cw));
            end
        end
    end
    
    methods (Access = protected)
        function setYN = CheckKind(obj)
            if ~isempty(obj.kind)
                disp(['Kind already set to ',obj.kind]);
                setYN = 1;
            else
                setYN = 0;
            end
        end
    end
end