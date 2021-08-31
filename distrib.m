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
