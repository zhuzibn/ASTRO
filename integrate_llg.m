mmx_=zeros(natomx,natomy,totstep);
mmy_=zeros(natomx,natomy,totstep);
mmz_=zeros(natomx,natomy,totstep);

BD=zeros(natomx,natomy,'gpuArray');
muigpu=zeros(natomx,natomy,'gpuArray');
scalgpu=zeros(natomx,natomy,'gpuArray');
AsimnextL=zeros(natomx,natomy,'gpuArray');
AsimnextW=zeros(natomx,natomy,'gpuArray');
AsimpreviousL=zeros(natomx,natomy,'gpuArray');
AsimpreviousW=zeros(natomx,natomy,'gpuArray');
for cty=1:natomy
    for ctx=1:natomx
        if atomtype_(ctx,cty)==1%RE
            muigpu(ctx,cty)=musRE;
            scalgpu(ctx,cty)=gamRE/(1+alp^2);%scale parameter
            BD(ctx,cty)=BDRE;
        else
            muigpu(ctx,cty)=musTM;
            scalgpu(ctx,cty)=gamTM/(1+alp^2);%scale parameter
            BD(ctx,cty)=BDTM;
        end
        if atomtype_(ctx,cty)==1%local atom is RE
            if cty==natomy
                AsimnextL(ctx,cty)=0;
            else
                if atomtype_(ctx,cty+1)==1%the other atom is RE
                    AsimnextL(ctx,cty)=Jgdgd;
                else%the other atom is TM
                    AsimnextL(ctx,cty)=Jfegd;
                end
            end
            
            if cty==1
                AsimpreviousL(ctx,cty)=0;
            else
                if atomtype_(ctx,cty-1)==1%the other atom is RE
                    AsimpreviousL(ctx,cty)=Jgdgd;
                else%the other atom is TM
                    AsimpreviousL(ctx,cty)=Jfegd;
                end
            end
            
            if ctx==natomx
                AsimnextW(ctx,cty)=0;
            else
                if atomtype_(ctx+1,cty)==1%the other atom is RE
                    AsimnextW(ctx,cty)=Jgdgd;
                else%the other atom is TM
                    AsimnextW(ctx,cty)=Jfegd;
                end
            end
            
            if ctx==1
                AsimpreviousW(ctx,cty)=0;
            else
                if atomtype_(ctx-1,cty)==1%the other atom is RE
                    AsimpreviousW(ctx,cty)=Jgdgd;
                else%the other atom is TM
                    AsimpreviousW(ctx,cty)=Jfegd;
                end
            end
            
        else%local atom is TM
            if cty==natomy
                AsimnextL(ctx,cty)=0;
            else
                if atomtype_(ctx,cty+1)==1%the other atom is RE
                    AsimnextL(ctx,cty)=Jfegd;
                else%the other atom is TM
                    AsimnextL(ctx,cty)=Jfefe;
                end
            end
            
            if cty==1
                AsimpreviousL(ctx,cty)=0;
            else
                if atomtype_(ctx,cty-1)==1%the other atom is RE
                    AsimpreviousL(ctx,cty)=Jfegd;
                else%the other atom is TM
                    AsimpreviousL(ctx,cty)=Jfefe;
                end
            end
            
            if ctx==natomx
                AsimnextW(ctx,cty)=0;
            else
                if atomtype_(ctx+1,cty)==1%the other atom is RE
                    AsimnextW(ctx,cty)=Jfegd;
                else%the other atom is TM
                    AsimnextW(ctx,cty)=Jfefe;
                end
            end
            
            if ctx==1
                AsimpreviousW(ctx,cty)=0;
            else
                if atomtype_(ctx-1,cty)==1%the other atom is RE
                    AsimpreviousW(ctx,cty)=Jfegd;
                else%the other atom is TM
                    AsimpreviousW(ctx,cty)=Jfefe;
                end
            end
        end
    end
end
clear ctx cty
BF=chi*BD;
ct3run=round((runtime)/gpusave);
ct3=1;
while ~(ct3>ct3run)
    
    mmx=zeros(natomx,natomy,gpusteps,'gpuArray');
    mmy=zeros(natomx,natomy,gpusteps,'gpuArray');
    mmz=zeros(natomx,natomy,gpusteps,'gpuArray');
    
    if ~(ct3==1)
        mmx(:,:,1)=tmp2xn0;mmy(:,:,1)=tmp2yn0;mmz(:,:,1)=tmp2zn0;
    else
        mmx(:,:,1)=mx_init;mmy(:,:,1)=my_init;mmz(:,:,1)=mz_init;
    end
    clear tmpx tmpy tmpz
    ct1=1; %count 1
    while ct1<gpusteps
        mmxtmp=mmx(:,:,ct1);
        mmytmp=mmy(:,:,ct1);
        mmztmp=mmz(:,:,ct1);
        if bc%not periodic condition
            mmxnextL=circshift(mmxtmp,[0,-1]);%length direction
            mmynextL=circshift(mmytmp,[0,-1]);
            mmznextL=circshift(mmztmp,[0,-1]);
            mmxpreviousL=circshift(mmxtmp,[0,1]);
            mmypreviousL=circshift(mmytmp,[0,1]);
            mmzpreviousL=circshift(mmztmp,[0,1]);
            mmxnextL(:,end)=0;mmynextL(:,end)=0;mmznextL(:,end)=0;
            mmxpreviousL(:,1)=0;mmypreviousL(:,1)=0;mmzpreviousL(:,1)=0;
            
            mmxnextW=circshift(mmxtmp,[-1,0]);%width direction
            mmynextW=circshift(mmytmp,[-1,0]);
            mmznextW=circshift(mmztmp,[-1,0]);
            mmxpreviousW=circshift(mmxtmp,[1,0]);
            mmypreviousW=circshift(mmytmp,[1,0]);
            mmzpreviousW=circshift(mmztmp,[1,0]);
            mmxnextW(end,:)=0;mmynextW(end,:)=0;mmznextW(end,:)=0;
            mmxpreviousW(1,:)=0;mmypreviousW(1,:)=0;mmzpreviousW(1,:)=0;
        else%periodic condition
            mmxnextL=circshift(mmxtmp,[0,-1]);
            mmynextL=circshift(mmytmp,[0,-1]);
            mmznextL=circshift(mmztmp,[0,-1]);
            mmxpreviousL=circshift(mmxtmp,[0,1]);
            mmypreviousL=circshift(mmytmp,[0,1]);
            mmzpreviousL=circshift(mmztmp,[0,1]);
        end
        
        hex_x=-(AsimnextL.*mmxnextL+AsimpreviousL.*mmxpreviousL+...
            AsimnextW.*mmxnextW+AsimpreviousW.*mmxpreviousW)./muigpu;%[T]
        hex_y=-(AsimnextL.*mmynextL+AsimpreviousL.*mmypreviousL+...
            AsimnextW.*mmynextW+AsimpreviousW.*mmypreviousW)./muigpu;
        hex_z=-(AsimnextL.*mmznextL+AsimpreviousL.*mmzpreviousL+...
            AsimnextW.*mmznextW+AsimpreviousW.*mmzpreviousW)./muigpu;
        
        hani_x=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');%anisotropy
        hani_y=zeros(size(hex_x,1),size(hex_x,2),'gpuArray');
        hani_z=2*Ksim./muigpu.*mmztmp;%[T]
        
        hdmi_x=Dsim./muigpu.*(-mmznextW+mmzpreviousW);%[T]
        hdmi_y=Dsim./muigpu.*(-mmznextL+mmzpreviousL);
        hdmi_z=Dsim./muigpu.*(mmxnextW-mmxpreviousW+mmynextL-mmypreviousL);
        
        if (0)%cpu calculation of dipole field, used for benchmaking
            cpuhdipolex=zeros(natomx,natomy);
            cpuhdipoley=zeros(natomx,natomy);
            cpuhdipolez=zeros(natomx,natomy);
            mmxttmp=gather(mmxtmp);
            mmyttmp=gather(mmytmp);
            mmzttmp=gather(mmztmp);
            for cty=1:natomy
                for ctx=1:natomx
                    for cty2=1:natomy
                        for ctx2=1:natomx
                            if ctx==ctx2 && cty==cty2
                                tmp=[0,0,0];
                            else
                                dist=sqrt(((ctx2-ctx)*d)^2+((cty2-cty)*d)^2);
                                rij=[(ctx2-ctx)*d,(cty2-cty)*d,0];
                                sj=[mmxttmp(ctx2,cty2),mmyttmp(ctx2,cty2),mmzttmp(ctx2,cty2)];
                                tmp=mu_0/(4*pi)*(-sj*gather(muigpu(ctx2,cty2))/dist^3+...
                                    3*rij*gather(muigpu(ctx2,cty2))*dot(sj,rij)/dist^5);
                            end
                            cpuhdipolex(ctx,cty)=cpuhdipolex(ctx,cty)+tmp(1);
                            cpuhdipoley(ctx,cty)=cpuhdipoley(ctx,cty)+tmp(2);
                            cpuhdipolez(ctx,cty)=cpuhdipolez(ctx,cty)+tmp(3);
                        end
                    end
                end
            end
            clear mmxttmp mmyttmp mmzttmp ctx cty ctx2 cty2 sj tmp
        end
        hdipolex_=zeros(natomx,natomy,'gpuArray');
        hdipoley_=zeros(natomx,natomy,'gpuArray');
        hdipolez_=zeros(natomx,natomy,'gpuArray');
        if dipolee
            for cty=1:natomy%gpu calculation of dipole field, used for production
                for ctx=1:natomx
                    [rijx_tmp,rijy_tmp]=meshgrid(0:natomy-1,0:natomx-1);
                    rijx_tmp=rijx_tmp-cty+1;%note:rijx_tmp corresponds to column
                    rijy_tmp=rijy_tmp-ctx+1;%note:rijy_tmp corresponds to row
                    dist_=0.4*1e-9*sqrt(rijx_tmp.^2+rijy_tmp.^2);
                    rijx_=0.4*1e-9*rijx_tmp;rijy_=0.4*1e-9*rijy_tmp;
                    dot_sr=muigpu.*(mmxtmp.*rijy_+mmytmp.*rijx_);
                    hdipolex=mu_0/(4*pi)*(3*rijy_.*dot_sr./dist_.^5-muigpu.*mmxtmp./dist_.^3);%[T]
                    hdipolex(ctx,cty)=0;
                    hdipoley=mu_0/(4*pi)*(3*rijx_.*dot_sr./dist_.^5-muigpu.*mmytmp./dist_.^3);
                    hdipoley(ctx,cty)=0;
                    hdipolez=mu_0/(4*pi)*(-muigpu.*mmztmp./dist_.^3);
                    hdipolez(ctx,cty)=0;
                    hdipolex_(ctx,cty)=sum(sum(hdipolex));
                    hdipoley_(ctx,cty)=sum(sum(hdipoley));
                    hdipolez_(ctx,cty)=sum(sum(hdipolez));
                end
            end
        end
        clear ctx cty ctx2 cty2
        hhx=hex_x+hani_x+hdmi_x+hdipolex_+Hext(1);
        hhy=hex_y+hani_y+hdmi_y+hdipoley_+Hext(2);
        hhz=hex_z+hani_z+hdmi_z+hdipolez_+Hext(3);
        if rk4==2%4th predictor-corrector
            if ct3==1 && ~(ct1>3)
                [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                    tstep,hhx,hhy,hhz);
            else
                [sxx,syy,szz]=arrayfun(@atomgpupc4,tmpxn0,tmpyn0,tmpzn0,...
                    tmpxn1,tmpyn1,tmpzn1,tmpxn2,tmpyn2,tmpzn2,tmpxn3,tmpyn3,tmpzn3,...
                    scalgpu,alp,tstep,hhx,hhy,hhz);
            end
        elseif rk4==1 %rk4
            [sxx,syy,szz]=arrayfun(@atomgpurk4,mmxtmp,mmytmp,mmztmp,psjSHEx,psjSHEy,psjSHEz,scalgpu,alp,...
                tstep,hhx,hhy,hhz,BD,BF);
        else%heun
            [sxx,syy,szz]=arrayfun(@atomgpu,mmxtmp,mmytmp,mmztmp,scalgpu,alp,...
                tstep,hhx,hhy,hhz);%
        end
        
        mmx(:,:,ct1+1)=sxx; mmy(:,:,ct1+1)=syy; mmz(:,:,ct1+1)=szz;
        ct1=ct1+1;
        if ~(ct3==1 && ~(ct1>3)) && ct1>3
            tmpxn0=mmx(:,:,ct1);tmpyn0=mmy(:,:,ct1);tmpzn0=mmz(:,:,ct1);
            tmpxn1=mmx(:,:,ct1-1);tmpyn1=mmy(:,:,ct1-1);tmpzn1=mmz(:,:,ct1-1);
            tmpxn2=mmx(:,:,ct1-2);tmpyn2=mmy(:,:,ct1-2);tmpzn2=mmz(:,:,ct1-2);
            tmpxn3=mmx(:,:,ct1-3);tmpyn3=mmy(:,:,ct1-3);tmpzn3=mmz(:,:,ct1-3);
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==2
            tmpxn0=mmx(:,:,ct1);tmpyn0=mmy(:,:,ct1);tmpzn0=mmz(:,:,ct1);
            tmpxn1=tmp2xn0;tmpyn1=tmp2yn0;tmpzn1=tmp2zn0;
            tmpxn2=tmp2xn1;tmpyn2=tmp2yn1;tmpzn2=tmp2zn1;
            tmpxn3=tmp2xn2;tmpyn3=tmp2yn2;tmpzn3=tmp2zn2;
        elseif ~(ct3==1 && ~(ct1>3)) && ct1==3
            tmpxn0=mmx(:,:,ct1);tmpyn0=mmy(:,:,ct1);tmpzn0=mmz(:,:,ct1);
            tmpxn1=mmx(:,:,ct1-1);tmpyn1=mmy(:,:,ct1-1);tmpzn1=mmz(:,:,ct1-1);
            tmpxn2=tmp2xn0;tmpyn2=tmp2yn0;tmpzn2=tmp2zn0;
            tmpxn3=tmp2xn1;tmpyn3=tmp2yn1;tmpzn3=tmp2zn1;
        end
    end
    tmp2xn0=mmx(:,:,end);tmp2yn0=mmy(:,:,end);tmp2zn0=mmz(:,:,end);
    tmp2xn1=mmx(:,:,end-1);tmp2yn1=mmy(:,:,end-1);tmp2zn1=mmz(:,:,end-1);
    tmp2xn2=mmx(:,:,end-2);tmp2yn2=mmy(:,:,end-2);tmp2zn2=mmz(:,:,end-2);
    mmx_(:,:,(ct3-1)*gpusteps+1:ct3*gpusteps)=gather(mmx);
    mmy_(:,:,(ct3-1)*gpusteps+1:ct3*gpusteps)=gather(mmy);
    mmz_(:,:,(ct3-1)*gpusteps+1:ct3*gpusteps)=gather(mmz);
    ct3=ct3+1;
end

clear mmx mmy mmz tmp2xn0 tmp2yn0 tmp2zn0 tmp2xn1 tmp2yn1 tmp2zn1
clear tmp2xn2 tmp2yn2 tmp2zn2
mmx=mmx_(:,:,1:savetstep:end);
mmy=mmy_(:,:,1:savetstep:end);
mmz=mmz_(:,:,1:savetstep:end);
clear mmx_ mmy_ mmz_
t=t(1:savetstep:end);
