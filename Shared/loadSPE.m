
% Read in SPE-file (ROPER SCIENTIFIC). Very simple function, does not
% really now about file format (reads in only UINT16, size of image has to
% be given).
%
% Size of image is given as [w, h]; if not given, [512 512] is assumed.
% offset = (12587012-12582912)/2;
%
% A. Bassi - Dipartamento di Fisica - Politecnico di Milano - ??/??/????
% N. Ducros- Dipartamento di Fisica - Politecnico di Milano - 11/02/2010
% A. Ghezzi- Dipartamento di Fisica - Politecnico di Milano - 14/09/2021

function out = loadSPE(file,image_size)

if nargin>=2   
    Nx = image_size(1); Ny = image_size(2);
    if length(image_size) == 2
        Nz = 1;
    else
        Nz = image_size(3);
    end
else
    Nx = 512; Ny = 512; Nz = 1;
end

fid     = fopen(file,'r');
im      = fread(fid,Nx*Ny*Nz,'uint16');
fclose(fid);
out = squeeze(permute(reshape(im,Nx,Ny,Nz),[2 1 3]));
return;
