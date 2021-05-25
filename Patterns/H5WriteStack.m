function [data] = H5WriteStack(filename,dataset,data,type)
%H5WRITESTACK Summary of this function goes here
%   Detailed explanation goes here
%dataset = 'image';
if ~exist(filename,'file')
    h5create(filename,dataset,size(data),'Datatype',type);
end
try
    h5info(filename,dataset)
catch
    h5create(filename,dataset,size(data),'Datatype',type);
end

h5write(filename,dataset,data)
% data = h5read(filename,[location,dataset]);% occhio a invertire X e Y dell'array ROI
% %,[654,522,1],[370,425,64]);
% labels = h5readatt(filename,[location,dataset],'DIMENSION_LABELS');
% el_size = h5readatt(filename,[location,dataset],'element_size_um');
end

