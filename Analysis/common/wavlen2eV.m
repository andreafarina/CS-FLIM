function [eV_axis] = wavlen2eV(lambda)
% lambda expressed in nm
coulomb = 1.602176634e-19;
h = 6.62607015e-34;
c = physconst('lightspeed');

eV_axis = h.*c./(coulomb.*lambda.*1e-9);
end

