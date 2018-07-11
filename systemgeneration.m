distrib();
mx_init=zeros(natomx,natomy);%initial magnetization
my_init=zeros(natomx,natomy);
mz_init=zeros(natomx,natomy);
if bc%
    if dwcalc
        phi_=0;
        for cty=1:natomy
            for ctx=1:natomx
                if cty<round(natomy/2)
                    if atomtype_(ctx,cty)==1%RE
                        thet_=5/180*pi;
                    else
                        thet_=(5+180)/180*pi;
                    end
                    mx_init(ctx,cty)=sin(thet_)*cos(phi_);
                    my_init(ctx,cty)=sin(thet_)*sin(phi_);
                    mz_init(ctx,cty)=cos(thet_);
                else
                    if atomtype_(ctx,cty)==1%RE
                        thet_=(5+180)/180*pi;
                    else
                        thet_=5/180*pi;
                    end
                    mx_init(ctx,cty)=sin(thet_)*cos(phi_);
                    my_init(ctx,cty)=sin(thet_)*sin(phi_);
                    mz_init(ctx,cty)=cos(thet_);
                end
            end
        end
    else
        phi_=0;
       for cty=1:natomy
           for ctx=1:natomx
               if atomtype_(ctx,cty)==1%RE
                   thet_=5/180*pi;
               else
                   thet_=(5+180)/180*pi;
               end
               mx_init(ctx,cty)=sin(thet_)*cos(phi_);
               my_init(ctx,cty)=sin(thet_)*sin(phi_);
               mz_init(ctx,cty)=cos(thet_);
           end
       end
    end
    clear cty ctx
else
    %to do
end