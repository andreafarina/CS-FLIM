function y = coates_correction_EA(current_curve_dis,num_chan,f_laser,T_int_s)
% current_curve_dist must have dimension 1 x num_chan
% num_chan is the total number of bins used by timing electronics
% f_laser is the frequency of the laser in Hz
% T_int_s is the acquisition time of the single curve in s

num_sweep=f_laser*T_int_s;
corrected_counts= zeros(num_chan,1);    
for ch_c= 1 : 1 : num_chan
    %calculation probability of an event occurring in channel i in one cycle
    summa_cor=0;
    for i= 1 : 1 : (ch_c-1)
        summa_cor= summa_cor + current_curve_dis(i);
    end
    prob_ch(1,ch_c)= current_curve_dis(ch_c)/(num_sweep-summa_cor);

    true_prob(1,ch_c)= - log(1-prob_ch(1,ch_c));
    true_counts(1,ch_c) = num_sweep * true_prob(1,ch_c);
    corrected_counts(ch_c,1) = true_counts(1,ch_c);
end
y = corrected_counts;
end
