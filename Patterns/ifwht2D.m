function y = ifwht2D(x,K)
    y = ifwht (ifwht (x)')./(K^2);
end