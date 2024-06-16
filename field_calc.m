mmxnextL=circshift(mmxtmp,[0,-1]);%length direction
mmynextL=circshift(mmytmp,[0,-1]);
mmznextL=circshift(mmztmp,[0,-1]);
mmxpreviousL=circshift(mmxtmp,[0,1]);
mmypreviousL=circshift(mmytmp,[0,1]);
mmzpreviousL=circshift(mmztmp,[0,1]);

mmxnextW=circshift(mmxtmp,[1,0]);%width direction
mmynextW=circshift(mmytmp,[1,0]);
mmznextW=circshift(mmztmp,[1,0]);
mmxpreviousW=circshift(mmxtmp,[-1,0]);
mmypreviousW=circshift(mmytmp,[-1,0]);
mmzpreviousW=circshift(mmztmp,[-1,0]);

if bc%not periodic condition
    mmxnextL(:,end)=0;mmynextL(:,end)=0;mmznextL(:,end)=0;
    mmxpreviousL(:,1)=0;mmypreviousL(:,1)=0;mmzpreviousL(:,1)=0;
    mmxnextW(1,:)=0;mmynextW(1,:)=0;mmznextW(1,:)=0;
    mmxpreviousW(end,:)=0;mmypreviousW(end,:)=0;mmzpreviousW(end,:)=0;
else%periodic condition
    %do nothing
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
% Renjie pointed out the error in calculating hdmi in the previous
% version of code, which has been corrected below.
% 2023.09, During the plot of skyrmion, Zhang Xue find the top-down
% direction is wrong. It is induced by the wrong nextW direction,
% which has been corrected.
hdmi_x=Dsim./muigpu.*(-mmznextL+mmzpreviousL);%[T]
hdmi_y=Dsim./muigpu.*(-mmznextW+mmzpreviousW);
hdmi_z=Dsim./muigpu.*(mmxnextL-mmxpreviousL+mmynextW-mmypreviousW);
if thermalenable
    %equation (15) in Atomistic spin model simulations of magnetic nanomaterials
    %calculate once for one time step
    Hthermal1=sqrt(2*kb*T*alp./(muigpu.*gamatom.*tstep));%[T]
    Hthermalx=normrnd(0,Hthermal1);
    Hthermaly=normrnd(0,Hthermal1);
    Hthermalz=normrnd(0,Hthermal1);
else
    Hthermalx=zeros(natomW,natomL,'gpuArray');
    Hthermaly=zeros(natomW,natomL,'gpuArray');
    Hthermalz=zeros(natomW,natomL,'gpuArray');
end
dipole_calc();
hhx=hex_x+hani_x+hdmi_x+hdipolex_+Hext(1)+Hthermalx;
hhy=hex_y+hani_y+hdmi_y+hdipoley_+Hext(2)+Hthermaly;
hhz=hex_z+hani_z+hdmi_z+hdipolez_+Hext(3)+Hthermalz;