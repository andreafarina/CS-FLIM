function y = ReadSDT_ContFlow(prefix_file,Ntime,Nlambda,Nbanks,Nsteps)
% read TCSPC data obtained with Continuous Flow operation

if Nbanks<10
    strform = '%01.f';
elseif Nbanks<100
    strform = '%02.f';
elseif Nbanks <1000
    strform = '%03.f';
elseif Nbanks <10000
    strform = '%04.f';
end
y = zeros(Ntime,Nlambda,Nbanks*Nsteps);
for i = 1:Nbanks
    d = f_read_sdt_01([prefix_file,'_c',num2str(i-1,strform),'.sdt']);
    %d = f_read_sdt_01([prefix_file,'.sdt']);
    d = reshape(d,[Ntime,Nlambda,Nsteps]);
    y(:,:,(1:Nsteps) + (i-1)*Nsteps) = d;
end
end