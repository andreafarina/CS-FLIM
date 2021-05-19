function t = TimeVectorTCSPC(n_chan)
    spc_gain = 2;
    factor = 50e3/n_chan/spc_gain; % 50e3 max intervallo TAC
    t = (1:n_chan)*factor - factor/2;
    t = t./1e3; %from ps to ns
end