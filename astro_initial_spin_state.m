function [mx_init, my_init, mz_init] = astro_initial_spin_state( ...
        atomtype_, natomW, natomL, dwcalc)
%ASTRO_INITIAL_SPIN_STATE Build ASTRO's deterministic initial spin arrays.

atomtype_ = astro_validate_atom_distribution(atomtype_, natomW, natomL);

mx_init = zeros(natomW, natomL);
my_init = zeros(natomW, natomL);
mz_init = zeros(natomW, natomL);
phi_ = 0;

for ctL = 1:natomL
    for ctW = 1:natomW
        is_tm = atomtype_(ctW, ctL) == 0;
        if dwcalc && ctL >= round(natomL/2)
            theta_degrees = 5 + 180*is_tm;
        else
            theta_degrees = 5 + 180*(~is_tm);
        end
        theta_ = theta_degrees/180*pi;
        mx_init(ctW, ctL) = sin(theta_)*cos(phi_);
        my_init(ctW, ctL) = sin(theta_)*sin(phi_);
        mz_init(ctW, ctL) = cos(theta_);
    end
end
end
