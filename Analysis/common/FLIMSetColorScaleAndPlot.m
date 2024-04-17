function [] = FLIMSetColorScaleAndPlot(lifetime,concentration)

N = size(lifetime,1);

CsumFast = squeeze(sum(concentration,1));
maskFast = find(CsumFast>mean(CsumFast(:)));
maskFast = reshape(ismember(1:length(CsumFast(:)),maskFast),[size(CsumFast,1) size(CsumFast,2)]);
CsumFast(1,1) = 0;

prompt = {'LowerBound:','UpperBound:'};
            dlgtitle = 'Coloscale Axis Limits:';
            dims = [1 35];
            definput = {num2str(min(lifetime(:))),num2str(max(lifetime(:)))};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            if size(answer,1) == 0
              answer = definput;
            end
ub = str2double(answer{2});
lb = str2double(answer{1});

lifetime(lifetime(:)>ub)=0;
concentration(lifetime(:)<lb)=0;
lifetime(lifetime(:)<lb)=0;

for li=1:size(lifetime,1)
    figure(size(lifetime,1)+li+50),hold on
    for n = 1:N
        subplot(1,N,n)
        taus = squeeze(lifetime(li,:,:,n));
        h = histogram(taus(taus>0),round(numel(lifetime(li,:,:,n))/70));
        h.BinEdges = (lb-0.05):0.1:(ub+0.05);
        plot(linspace(mean([h.BinEdges(1:2)]),mean([h.BinEdges(end-1:end)]),h.NumBins),h.Values,'o');
    end
end

lifetime(1,1,32) = ub;
lifetime(1,1,31) = lb;
DrawFLIMMap(lifetime(1,:,:),permute(CsumFast.*maskFast,[3 1 2]),N,zeros(1,size(lifetime,2),size(lifetime,3)));
end