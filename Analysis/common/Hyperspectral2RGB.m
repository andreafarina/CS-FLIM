function [S] = Hyperspectral2RGB(lambda,im)

% Tool to render an RGB image from a multispectral data
% Alberto Ghezzi

if lambda(1) < 380 || lambda(end) > 780
    disp('Wavelength range out of visible range')
    return
end

    [r,g,b] = Wavelength2color(lambda);
    
    for li = 1:size(im,2)     
        I = squeeze(sum(im(:,li,:,:),1)); 
        S_r(:,:,li) = r(li).*I;
        S_g(:,:,li) = g(li).*I;
        S_b(:,:,li) = b(li).*I;
    end
    S_r = sum(S_r,3);
    S_g = sum(S_g,3);
    S_b = sum(S_b,3);
    
    S_r_n = 255.*S_r./max([max(S_r,[],1:2) max(S_g,[],1:2) max(S_b,[],1:2)]);
    S_g_n = 255.*S_g./max([max(S_r,[],1:2) max(S_g,[],1:2) max(S_b,[],1:2)]);
    S_b_n = 255.*S_b./max([max(S_r,[],1:2) max(S_g,[],1:2) max(S_b,[],1:2)]);
    
    S = uint8(reshape([S_r_n S_g_n S_b_n],[size(im,3:4) 3]));
end

function [x,y,z] = Wavelength2color(a)
    lambda = [380 420 440 490 510 580 645 780];
    r = [97 106 0 0 0 255 255 97]./255;
    g = [0 0 0 255 255 255 0 0]./255;
    b = [97 255 255 255 0 0 0 0]./255;
    x = interp1(lambda,r,a);
    y = interp1(lambda,g,a);
    z = interp1(lambda,b,a);
end