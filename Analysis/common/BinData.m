function [t,data,dt] = BinData(data_nobin,t_nobin, dt,bin_type)
if nargin<4
    bin_type = 'Bin';
end
data = data_nobin;
t = t_nobin;
bin = round(length(t)/(dt/(t(2)-t(1))));

if strcmp(bin_type,'Bin') == 1
    if bin < length(t)
        N = length(data(:,1,1,1));
        K = 1:N;
        D = K(rem(N,K)==0);
        [~,p] = min(abs(bin-D));
        bins = D(p);
        subs = discretize(1:length(data(:,1,1,1)),bins)';
        
        for li = 1:size(data,2)
            for xi = 1:size(data,3)
                for yi = 1:size(data,4)
                    data_1(:,li,xi,yi) = accumarray(subs,data(:,li,xi,yi),[],@sum);
                end
            end
        end
        data = data_1(:,:,:,:);
        if (size(discretize(t,bins),1) == size(t,1)) && (size(t,1)>1)
            t = accumarray(discretize(t,bins),t,[],@mean);
        else
            t = accumarray(discretize(t,bins)',t,[],@mean)';
        end
        %t = t(2:end-1);
    end
    if abs((t(2)-t(1))-dt)>(dt/2)
        disp('Some problems determining the desired bin size.')
    end
    dt = t(2)-t(1);
     
else
    
    for li = 1:size(data,2)
        for xi = 1:size(data,3)
            for yi = 1:size(data,4)
                data(:,li,xi,yi) = movmean(data(:,li,xi,yi),round(length(t)/bin));
            end
        end
    end
    dt = t(round(length(t)/bin))-t(1)+t(2)-t(1);
end
end