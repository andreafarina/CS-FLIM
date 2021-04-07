function [img] = PointGrayReadRaw(name)
img = multibandread(name, [1200 1920 3],'uint8=>uint8',0, 'bip', 'ieee-be');
img = flipud(img);
end