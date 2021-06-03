function y = readBinaryMANTA(filename)
fid = fopen(filename,'r','ieee-be');
nlambda = fread(fid,1,'uint32');
ntime = fread(fid,1,'uint32');
y = uint16(fread(fid,[ntime,nlambda],'double'));
fclose(fid);
end   