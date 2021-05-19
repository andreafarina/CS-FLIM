function [im] = SHReconstruction(had,Pc,Pr,flag_flip)
% performs a ScrambledHadamard inversion on a had(time,lambda,coeff) dataset.
% flag_flip is a tick that flips left-right the reconstructed image
L = size(had,2);
T = size(had,1);
for li = 1:L
    parfor ti = 1:T
        x = Pc'*fwht(eye(sqrt(size(had,3))^2),[],'hadamard')*Pr'*squeeze(had(ti,li,:));
        A = reshape(x,[sqrt(size(had,3)) sqrt(size(had,3))]);
        if flag_flip == 1
            im(ti,li,:,:) = fliplr(A);
        else
            im(ti,li,:,:) = A;
        end
    end
    disp(['Progress...',num2str(round(100*li/L)),'%'])
end
end