function [h] = DrawFLIMMap(lifetime,concentration,N,residuals)

% Plot the results
for li = 1:size(lifetime,1)
    for n = 1:N
        h = figure(li+50); subplot(1,N+1,n),
        % Colormap is from 0 to 1. To avoid 1s and 0s to be both
        % colored in red,  The limit
        % is 0°,260°. Then I set the 0 to the min lifetime
        A = squeeze(lifetime(li,:,:,n));
        B = A-min(min(A(A>0)));
        B(B<0) = 0;
        I(:,:,1) = B./max(max(B));
        % Additionally, I want shot times in blue and long in
        % red:
        %I(:,:,1) =  1-I(:,:,1);
        %I limit the color wheel to a degree
        % lower than 360° (both 0° and 360° are red).
        I(:,:,1) = I(:,:,1).*(260/360);
        I(:,:,2) = ones([size(I,1) size(I,2)]);
        I(:,:,3) = squeeze(concentration(li,:,:,n)./max(concentration(:)));
        
        imagesc((hsv2rgb(I))),axis image,
        hsv = [flipud(linspace(0,260/360,255)); ones(1,255);ones(1,255)]';
        rgb = hsv2rgb(hsv);
        %colormap(flipud(rgb))
        colormap(rgb)
        cb = colorbar;
        cb.Label.String = 'Lifetime [ns]';
        cb.Ticks = linspace(0,1,8);
        cb.TickLabels = round(linspace(min(min(A(A>0))),max(max(A)),length(cb.Ticks)),1);
        %title(['Composite lifetime map for \tau_',num2str(n)])
    end
    if nargin == 4
    resplot = subplot(1,N+1,N+1);
    Z = squeeze(residuals(li,:,:));
    h = imagesc(Z);axis image,colormap(resplot,'parula'),colorbar,caxis([0 1])
    % Set nan values to transparent: 
    set(h,'alphadata',Z>0)
    % Make the background black: 
    set(gca,'color','black')
    title('Residuals')
    end
end


for li=1:size(lifetime,1)
    figure(size(lifetime,1)+li+50)
    for n = 1:N
        subplot(1,N,n)
        taus = squeeze(lifetime(li,:,:,n));
        histogram(taus(taus>0),round(numel(lifetime(li,:,:,n))/70))
    end
end

end
