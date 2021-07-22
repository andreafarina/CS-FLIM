function [] = Convert_sdt2ASCII(filename,spc,t,lambda)
% filename must be included without extension
M = [t' spc];
FID=fopen([filename,'.txt'],'w');
fprintf(FID,'SPC\n');
fprintf(FID,'additional comments\n');
fprintf(FID,'Wavelength explicit\n');
fprintf(FID,['intervalnr ',num2str(size(spc,2)),'\n']);
fprintf(FID,'\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f',round(lambda,3));
writematrix(M,[filename,'.txt'],'WriteMode','append','Delimiter','tab');

fclose(FID);
movefile([filename,'.txt'],[filename,'.ascii'])

end