function y = fwht2D(x,K)
    y = fwht (fwht (x)')*K^2;
end