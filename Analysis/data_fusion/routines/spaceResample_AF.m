function [ hypcube ] = spaceResample_AF( input, newSize )
%Given a 4D (lambda,t,x,y) tensor [input], spaceResample provides a
%spatially downsampled/upsampled 4D tensor with a size of [newSize]. Also
%normalizes the 4D hyperspectral tensor so norm(hypcube(:)) = 1;
global M Mt
%preallocating
hypcube = zeros(size(input,1),size(input,2),newSize,newSize);
oldSize = size(input,3);
if (newSize)>(oldSize)
    Mop = M;
else
    Mop = Mt;
end
%downsample process
for i=1:size(hypcube,1)
    for j=1:size(hypcube,2)
        d = input(i,j,:,:);
        temp = Mop * d(:);
        hypcube(i,j,:,:) = reshape(temp,[newSize,newSize]);
    end
end

end