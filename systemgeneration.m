distrib();
mx_init=zeros(natomW,natomL);%initial magnetization
my_init=zeros(natomW,natomL);
mz_init=zeros(natomW,natomL);

if dwcalc
    phi_=0;
    for ctL=1:natomL
        for ctW=1:natomW
            if ctL<round(natomL/2)
                if atomtype_(ctW,ctL)==0%TM
                    thet_=5/180*pi;
                else
                    thet_=(5+180)/180*pi;
                end
                mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
                my_init(ctW,ctL)=sin(thet_)*sin(phi_);
                mz_init(ctW,ctL)=cos(thet_);
            else
                if atomtype_(ctW,ctL)==0%TM
                    thet_=(5+180)/180*pi;
                else
                    thet_=5/180*pi;
                end
                mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
                my_init(ctW,ctL)=sin(thet_)*sin(phi_);
                mz_init(ctW,ctL)=cos(thet_);
            end
        end
    end
else
    phi_=0;
    for ctL=1:natomL
        for ctW=1:natomW
            if atomtype_(ctW,ctL)==0%TM
                thet_=5/180*pi;
            else
                thet_=(5+180)/180*pi;
            end
            mx_init(ctW,ctL)=sin(thet_)*cos(phi_);
            my_init(ctW,ctL)=sin(thet_)*sin(phi_);
            mz_init(ctW,ctL)=cos(thet_);
        end
    end
end
clear ctL ctW
