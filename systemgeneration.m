distrib();
atomtype_=astro_validate_atom_distribution(atomtype_,natomW,natomL,compositionn);
[mx_init,my_init,mz_init]=astro_initial_spin_state(atomtype_,natomW,natomL,dwcalc);

if loadstartm
    clear mx_init my_init mz_init
    load(startmname);
    if natomW~=natomWcheck || natomL~=natomLcheck
       error('system not consistent') 
    end
    clear natomxcheck natomycheck
    mx_init=mmxstart;
    my_init=mmystart;
    mz_init=mmzstart;
end
