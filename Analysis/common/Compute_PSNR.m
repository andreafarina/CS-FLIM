% Compute_PSNR -- Compute the Peak Signal to Noise Ratio between two

% images, one being the reference image.

%

%  --> Inputs

%       I_ref       reference image

%       I           image to compare

%

%  --> Ouputs

%       PSNR        PSNR in decibels between I_ref and I

%       MSE         Mean Square Error between I_ref and I

%

%  --> Usage

%      [PSNR,MSE] = Compute_PSNR(I_ref,I)

%

%   Author : F. Rousset - 02/12/14



function [PSNR_val,MSE,S] = Compute_PSNR(I_ref,I)



[H,W] = size(I_ref);



if size(I,1) ~= H || size(I,2) ~= W

    error('Images must have the same size !');

end





min_I_ref = min(I_ref(:));

max_I_ref = max(I_ref(:));

min_I = min(I(:));

max_I = max(I(:));



I = (max_I_ref - min_I_ref) * (I - min_I) / (max_I - min_I) + min_I_ref;

vect_diff = I_ref(:) - I(:);



MSE = 1/length(vect_diff) * sum(vect_diff .^ 2);

PSNR_val = 10*log10((max_I_ref - min_I_ref)^2 / MSE); %num: dynamic range
%PSNR_val = 10*log10((max_I_ref)^2 / MSE);

S = ssim(I,I_ref);
end



