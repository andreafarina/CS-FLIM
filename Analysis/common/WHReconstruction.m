function [im] = WHReconstruction(had,flag_flip)
% performs a WH inversion on a had(time,lambda,coeff) dataset.
% flag_flip is a tick that flips left-right the reconstructed image
L = size(had,2);
T = size(had,1);
c2d = reshape(had,[size(had,1) size(had,2) sqrt(size(had,3)) sqrt(size(had,3))]);
parfor ti = 1:T
    for li = 1:L
        A = ifwht2D(squeeze(c2d(ti,li,:,:)),sqrt(size(had,3)));
        if flag_flip == 1
            im(ti,li,:,:) = fliplr(A);
        else
            im(ti,li,:,:) = A;
        end
    end
end
end