function M = SpaceResampleMatrix(oldSize,newSize)
%% direct operator
Nold = prod(oldSize);
Nnew = prod(newSize);
M = sparse(zeros(Nnew,Nold));
for i = 1:Nold
    d = zeros(oldSize);
    d(i) = 1;
    dL = imresize(d,newSize,'box');
    M(:,i) = dL(:);
end