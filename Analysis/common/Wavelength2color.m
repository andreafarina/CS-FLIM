function [x,y,z] = Wavelength2color(a)
    lambda = [380 420 440 490 510 580 645 780];
    r = [97 106 0 0 0 255 255 97]./255;
    g = [0 0 0 255 255 255 0 0]./255;
    b = [97 255 255 255 0 0 0 0]./255;
    x = interp1(lambda,r,a);
    y = interp1(lambda,g,a);
    z = interp1(lambda,b,a);
end