function [im] = SHReconstruction(had,M,flag_flip,Pc,Pr)
% performs a ScrambledHadamard inversion on a had(time,lambda,coeff) dataset.
% flag_flip is a tick that flips left-right the reconstructed image
L = size(had,2);
T = size(had,1);

%B = fwht(eye(sqrt(size(had,3))^2),[],'hadamard');
for li = 1:L
    parfor ti = 1:T
        % x = Pc'*B*Pr'*squeeze(had(ti,li,:));
        [x,flag] = lsqr(M,squeeze(had(ti,li,:)));
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

