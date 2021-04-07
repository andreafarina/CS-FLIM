function y = EqualizeEfficiency(data,range)
% data: [Nlambda,Neff]
% range: [i1,j1;
%        [i2, j2;...]
% where j1 = i2, etc...
% equalize efficiency curves taken with different reference samples
Ncurve = size(data,2);
if Ncurve < 2
    disp('No curves ti equalize, only 1 curve');
    y = data;
    return;
end
%% sort range and data in ascending order
[range,iord] = sort(range);
%% check consistency
for i = 1:Ncurve-1
    if range(i,2)~=range(i+1,1)
        disp('intervals not correctly joint');
        return;
    end
end
%% reorder data
data = data(:,iord(:,1)');
for i = 2:Ncurve
    y = data(:,i-1);
    fac = data(range(i-1,2),i-1)./data(range(i,1),i);
    y(range(i,1):range(i,2)) = data(range(i,1):range(i,2),i).*fac;
end