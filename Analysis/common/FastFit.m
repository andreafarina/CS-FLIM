function [A] = FastFit(meas,DAS, T)

A = zeros(size(meas,3),size(DAS,1));

%wlTmap1 = T(:,1)*DAS(1,:);
%wlTmap2 = T(:,2)*DAS(2,:);
%G = repmat(permute(A,[3 1 2]),[2 1 1]);


tic
% For each pattern, fit the map to the experimental one to retrieve the weight of each pattern
for l = 1:size(meas,3)
    %G = squeeze(meas(:,:,l));
    A(l,:) = diag(pinv(T)*squeeze(meas(:,:,l))*pinv(DAS));
    %A(l,:) = diag([wlTmap1(:)';wlTmap2(:)']'\[reshape(G,1,[]); reshape(G,1,[])]');
end
toc

end