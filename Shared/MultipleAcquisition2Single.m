function [] = MultipleAcquisition2Single(filename,N,K)

% Takes a file .mat acquired with multiple short acquisitions and sums
% N is the number of measurements to be summed
% K is the number of pattern per measurement

% This script must be run from the folder containg the data

disp('Loading file...')
dato = load(filename);
disp('Loading file: DONE')

A = dato.spc(:,:,2:end);
A = reshape(A,[size(A,1) size(A,2) K size(A,3)/K]);

if size(A,4)<N
    N = size(A,4);
    disp('Warning: N is greater than the acquired measurement')
end

sumA = zeros(size(A,1), size(A,2), K);
for i = 1:N
  sumA = sumA + A(:,:,:,i);  
end

B = dato.spc(:,:,1);

spc = zeros(size(A,1), size(A,2),K+1);
spc(:,:,1) = B;
spc(:,:,2:end) = sumA;

out_path = './summed/';
if exist(out_path,'file') == 0
    mkdir(out_path)
end

save([out_path, filename,'_summedover',num2str(N)],'spc')

end