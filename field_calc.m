deterministic_fields=astro_evaluate_field_terms(mmxtmp,mmytmp,mmztmp,...
    AsimnextL,AsimpreviousL,AsimnextW,AsimpreviousW,muigpu,Ksim,Dsim,Hext,bc);
hex_x=deterministic_fields.exchange.x;
hex_y=deterministic_fields.exchange.y;
hex_z=deterministic_fields.exchange.z;

hani_x=deterministic_fields.anisotropy.x;%anisotropy
hani_y=deterministic_fields.anisotropy.y;
hani_z=deterministic_fields.anisotropy.z;%[T]
% Renjie pointed out the error in calculating hdmi in the previous
% version of code, which has been corrected below.
% 2023.09, During the plot of skyrmion, Zhang Xue find the top-down
% direction is wrong. It is induced by the wrong nextW direction,
% which has been corrected.
hdmi_x=deterministic_fields.dmi.x;%[T]
hdmi_y=deterministic_fields.dmi.y;
hdmi_z=deterministic_fields.dmi.z;
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
hhx=deterministic_fields.total.x+hdipolex_+Hthermalx;
hhy=deterministic_fields.total.y+hdipoley_+Hthermaly;
hhz=deterministic_fields.total.z+hdipolez_+Hthermalz;
