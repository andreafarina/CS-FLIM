function [tau,DAS,T,spcFit] = GlobalFit(t,lambda,spc,tau_0,NON_NEG,showplots)
global wl time wlTmap flag neg

wl = lambda;
time = t;
wlTmap = spc;
flag = false;
neg = NON_NEG;

if  nargin<4
    % Initial Condition
    disp('NOT ENOUGH INPUTS, RUN SIMULATION')
    time = TimeVectorTCSPC(256);
    wl = 1:16;
    DAS(:,1) = 300.*normpdf(wl,wl(4),2);
    DAS(:,2) = 100.*normpdf(wl,wl(13),2);
    tau_0 = [0.6 3];
    k0 = 1./tau_0;
    ref = poissrnd(exp(-time'*k0)*DAS');
    wlTmap = ref';
    tau_0 = [1 2];
end

k0 = 1./tau_0;

optimfun = @(x,x_data) global_analysis(x);

options = optimset('Display','off',...
    'Jacobian','off',...
    'DerivativeCheck','off',...
    'MaxIter',10000,...
    'MaxFunEvals',10000,...
    'TolFun',1e-10,'TolX',1e-10); %

xdata = 1:length(wl)*length(time);
lb = zeros(size(k0));
ub = 100.*ones(size(k0));
tic
kf = lsqcurvefit(optimfun,k0,xdata,wlTmap(:),lb,ub,options);
[~,b] = unique(round(100.*kf));
if length(b)<length(kf)
kf = kf(b);
end
toc
flag = showplots;
[~,T,DAS] = global_analysis(kf);
spcFit = (T*DAS)';
tau = 1./kf;

if flag
    figure(3), subplot(1,3,1);
    imagesc(wl,time,wlTmap'); xlabel('Wavelength'); ylabel('Time');
    xlim([lambda(1)-(lambda(2)-lambda(1)) lambda(end)+(lambda(2)-lambda(1))])
    box on; axis tight square;
    title('Data'); cb = colorbar;
    A = str2double(cb.TickLabels);
    if 10*round(max(wlTmap(:))/10)-A(end)>mean(diff(A))/2
        upvalue = 10*round(max(wlTmap(:))/10)+mean(diff(A))/2;
    else
        upvalue = 10*round(max(wlTmap(:))/10);
    end
    caxis([0 upvalue]);
    subplot(1,3,2);
    imagesc(wl,time,spcFit');xlabel('Wavelength'); ylabel('Time');
    xlim([lambda(1)-(lambda(2)-lambda(1)) lambda(end)+(lambda(2)-lambda(1))])
    grid on; box on; axis tight square;
    title('Fit'); cb = colorbar;
    caxis([0 upvalue]);
    
    subplot(1,3,3);
    B = wlTmap'-spcFit';
    imagesc(wl,time,B);xlabel('Wavelength'); ylabel('Time');
    xlim([lambda(1)-(lambda(2)-lambda(1)) lambda(end)+(lambda(2)-lambda(1))])
    grid on; box on; axis tight square;
    title('Residuals');cb = colorbar;
    caxis(round([-3*sqrt(upvalue) 3*sqrt(upvalue)]));
    
    figure,
    for li = 1:length(wl)
        subplot(4,4,li), plot(time,wlTmap(li,:),'.'),hold on,plot(time,spcFit(li,:)),title([num2str(round(wl(li))), 'nm'])
    end
end

end

function [diff,T,S] = global_analysis(x)
global wl time wlTmap flag neg
T = exp(-x.*time');
% if neg
% ops = optimset('Display','off');
% A = repmat(wlTmap,[size(T,2) size(T,2)])';
% x = lsqlin(A,T(:),[],[],repmat(eye(size(wlTmap,1)),[1 size(T,2)]),sum(wlTmap,2)./sum(T(:)),zeros(size(A,2),1),1e6.*ones(size(A,2),1),[],ops);
% S = reshape(x,[size(wlTmap,1) size(T,2)])';
% else
S = T\wlTmap';
% end
if neg
S(S<0) = 0;
end
Dfit = (T*S)';

if flag
    figure(2); clf;
    subplot(2,1,1);
    plot(wl,S);
    title('Decay Associated Spectra');
    grid minor; axis tight;
    xlabel('Wavelength [nm]');
    ylabel('Amplitude [a.u.]')
    xlim([wl(1)-(wl(2)-wl(1)) wl(end)+(wl(2)-wl(1))])
    L = repmat({'DAS'},[1 size(S,1)]);
    for i = 1:size(S,1)
       L{i}(4) = num2str(i);
    end
    legend(L)

    subplot(2,1,2);
    plot(time,T);
    title('Temporal evolution');
    grid on; axis tight;
    xlabel('time');
    ylabel('amplitude')

    drawnow;
end
diff = Dfit(:);
end
