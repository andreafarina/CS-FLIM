function [t,data,dt] = BinData(data_nobin,t_nobin, dt)
data = data_nobin;
t = t_nobin;      
            bin = round(length(t)/(dt/(t(2)-t(1))));
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
                
                t = accumarray(discretize(t,bins)',t,[],@mean)';
                %t = t(2:end-1);
            end
            dt = t(2)-t(1);
end