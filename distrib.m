%% distribute atoms
%1 RE,0 TM
tmp=randperm(natom,round(natom*compositionn));
atomtype_=10*ones(natomW,natomL);
for ctL=1:natomL
    for ctW=1:natomW
        atomtype_(ctW,ctL)=ismember(ctW+(ctL-1)*natomW,tmp);
    end
end
clear tmp ctW ctL ctz
if load_fixed_atom_distrib
    clear atomtype_
    load('atomtypee.mat');
elseif save_fixed_atom_distrib%save debug data
    save('atomtypee.mat','atomtype_');%change this to save(ddebugfilename);
    error('distribution mat file has been saved, run the program again by setting load_fixed_atom_distrib=1')
end