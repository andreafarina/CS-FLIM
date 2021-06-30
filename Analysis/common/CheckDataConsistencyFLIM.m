function [lifetime,concentration] = CheckDataConsistencyFLIM(lifetime,concentration,Athreshold,Tthreshold)
        %% Discard values below a threshold
        %Athreshold = Athreshold/100;
        Tthreshold = Tthreshold/100;
        % values below a threshold are 0. to move the zero lifetime in the last page:
        lifetime(concentration(:)<=Athreshold) = 1e5;
        concentration(concentration(:)<=Athreshold) = 0;
    
    %% sort tau1 and tau2 according to tau1<tau2
    for li = 1:size(lifetime,1)
       taus = squeeze(lifetime(li,:,:,:));
       Cons = squeeze(concentration(li,:,:,:));
       [B,I] = sort(taus,3);
       
       r = repmat(1:size(taus,1),[1 size(taus,2)*size(taus,3)]);
       c = repmat(1:size(taus,1),[size(taus,2) size(taus,3)]);
       c = c(:)';
       p = I(:)';
       s = sub2ind(size(taus),r,c,p);
       % Sort concentratinos according to tau sorting
       As = reshape(Cons(s),size(taus));
       
       lifetime(li,:,:,:) = B;
       concentration(li,:,:,:) = As;     
    end
    lifetime(concentration(:)<=Athreshold) = 0;
    %% Check data similarity
    % if tau1 and/or tau2 and/or tau3 differ for less than 1%, consider it
    % as a monoexp. Sum the amplitudes and take the mean for the lifetime
    for li = 1:size(lifetime,1)
       for xi = 1:size(lifetime,2)
       for yi = 1:size(lifetime,3)
         taus = squeeze(lifetime(li,xi,yi,:));
         taus = diff(taus)./taus(1:end-1);
         taus(abs(taus)<Tthreshold) = 0;
         taus(abs(taus)>=Tthreshold) = 1;
         pos = find(taus==0);
         if pos == 1
             taus = [0; taus];
             val_t= [mean(lifetime(li,xi,yi,logical(1-taus))), squeeze(lifetime(li,xi,yi,logical(taus)))];
             val_t = [val_t zeros(1,length(taus)-length(val_t))];
             lifetime(li,xi,yi,:) = val_t;
             val_c = [sum(concentration(li,xi,yi,logical(1-taus))), squeeze(concentration(li,xi,yi,logical(taus)))];
             val_c = [val_c zeros(1,length(taus)-length(val_c))];
             concentration(li,xi,yi,:) = val_c;
         elseif pos == 2
             taus = [taus; 0];
             val_t= [squeeze(lifetime(li,xi,yi,logical(taus))), mean(lifetime(li,xi,yi,logical(1-taus)))];
             val_t = [val_t zeros(1,length(taus)-length(val_t))];
             lifetime(li,xi,yi,:) = val_t;
             val_c = [squeeze(concentration(li,xi,yi,logical(taus))), sum(concentration(li,xi,yi,logical(1-taus)))];
             val_c = [val_c zeros(1,length(taus)-length(val_c))];
             concentration(li,xi,yi,:) = val_c;   
         end
       end
       end
    end


end