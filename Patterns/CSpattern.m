classdef CSpattern
    %PATTTERN this class creates/allocate/save 2D
    % squared patterns
    %   Detailed explanation goes here
    
    properties
        type
        order
        Npatt
        posneg = false;
        dim
        stack
        Mmatrix
        % random seed
        seed = 122;
        % scrambled dedicated
        Pr
        Pc
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
                    obj.stack = reshape(obj.Mmatrix,[order,order,order^2]);
                    obj.order = order;
                
                case 'wh'
                    for i = 1:order^2
                        d = zeros(order,order);
                        d(i) = 1;
                        hk = fwht2D(d,order);
                        obj.stack(:,:,i) = hk;
                        obj.Mmatrix(i,:) = hk(:);
                    end
                    obj.dim = [order,order];
                    obj.order = order;
                
                case 'scrambled'
                    obj.dim = [order,order];
                    K = order^2;
                    H = hadamard(K);
                    I = eye(K);
                    obj.Pr = I(randperm(K),:);
                    obj.Pc = I(randperm(K),:);
                    obj.Mmatrix = obj.Pr * H * obj.Pc;
                    obj.stack = reshape(obj.Mmatrix,[order,order,order^2]);
                    obj.order = order;
                    
                case 'bernoulli'
                    obj.dim = [order,order];
                    obj.Mmatrix = random('bino',1,1/2,[order^2,order^2])*2 - 1;
                    obj.stack = reshape(obj.Mmatrix,[order,order,order^2]);
                    obj.order = order;
                
                case 'fourier'
                    % to be implemented
            end
            obj.Npatt = size(obj.Mmatrix,1);
            obj.type = type;
        end
        
        function obj = resample(obj,dim)
            rstack = zeros(dim(1),dim(2),obj.Npatt);
            for i = 1:obj.Npatt
                rstack(:,:,i) = imresize(obj.stack(:,:,i),dim,'box');
            end
            obj.stack = rstack;
            obj.Mmatrix = (reshape(rstack,[prod(dim),obj.Npatt]))';
            obj.dim = dim;
        end
        
        function obj = TrueToPosNeg(obj)
            nrow = obj.Npatt;
            Mpos = obj.Mmatrix;
            Mpos(Mpos<0) = 0;
            Mneg = -obj.Mmatrix;
            Mneg(Mneg<0) = 0;
            M = [Mpos;Mneg];
            % permutatio matrix
            Prr = zeros(2*nrow);
            Prr(1:nrow,1:2:end) = eye(nrow);
            Prr((1:nrow)+nrow,2:2:end) = eye(nrow);
            obj.Mmatrix = Prr' * M;
            obj.stack = reshape(obj.Mmatrix',[obj.dim(1),obj.dim(2),2*obj.Npatt]);
            obj.Npatt = 2 * obj.Npatt;
            obj.posneg = true;
        end
        
        function obj = TrueToShifted(obj)
            obj.Mmatrix = obj.Mmatrix + abs(min(obj.Mmatrix,[],2));
            obj.Mmatrix = obj.Mmatrix./max(obj.Mmatrix,[],2);
            obj.stack = reshape(obj.Mmatrix',[obj.dim(1),obj.dim(2),obj.Npatt]);
        end
        
        function [] = showPatterns(obj,space)
            %SHOWPATTERNS show patterns as mosaic
            %   Detailed explanation goes here
            if nargin<2
                space = 1;
            end
            Nx = size(obj.stack,1);
            Ny = size(obj.stack,2);
            stride = Nx + space;
            A = mean(obj.stack(:))*ones(stride*round(sqrt(obj.Npatt)),stride*round(sqrt(obj.Npatt)));
            figure;ax = gca;
            for ip = 1:obj.Npatt
                d = zeros(stride,stride);
                d(space + 1:stride,1:Ny) = obj.stack(:,:,ip);
                id = zeros(round(sqrt(obj.Npatt)));
                id(ip) = 1;
                A = A + kron(id',d);
                imagesc(ax,A);
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
            end
        end
    end
end
    
