function [had,position_to_survive] = PutCofficients2Zero(had,ratio)
% According to the amplitude of Hadamard coefficients, put to zero the
% lowest according to a number defined percentually (0,1) by "ratio"
compress_had = had;
had = zeros(size(had));
[~,ordered_coeff] = sort(abs(squeeze(sum(compress_had,1:2))),'descend');
ordered_coeff = ordered_coeff(1:round(size(had,3)*(1-ratio/100)),1);
pos_survive_coeff = ismember(1:size(had,3),ordered_coeff)';
position_to_survive = repmat(pos_survive_coeff,[size(had,1) size(had,2) 1]);
had = reshape(compress_had(:).*position_to_survive(:),[size(had,1) size(had,2) size(had,3)]);
end