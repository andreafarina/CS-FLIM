function [lifetime,concentration,N] = Convert2TauAverage(lifetime,concentration)
   A = concentration./sum(concentration,4);
   A(isnan(A)) = 0;
   lifetime = sum(lifetime.*A,4);
   concentration = sum(concentration,4);
   N = 1;
end