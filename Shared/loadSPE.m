
% Read in SPE-file (ROPER SCIENTIFIC). Very simple function, does not
% really now about file format (reads in only UINT16, size of image has to
% be given).
%
% Size of image is given as [w, h]; if not given, [512 512] is assumed.
% offset = (12587012-12582912)/2;
%
% A. Bassi - Dipartamento di Fisica - Politecnico di Milano - ??/??/????
% N. Ducros- Dipartamento di Fisica - Politecnico di Milano - 11/02/2010

function out = loadSPE(file,image_size)

if nargin>=2,   Nx = image_size(1); Ny = image_size(2);
else            Nx = 512; Ny = 512;
end;

fid     = fopen(file,'r');
im      = fread(fid,Nx*Ny,'uint16');
fclose(fid);
out = reshape(im,Nx,Ny)';
return;
