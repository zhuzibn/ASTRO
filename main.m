% 2D atomistic model for FiM
% require nvidia GPU
clear all;clc;close all;tic
%% control parameter
ddebug=0;
loadstartm=0;%1:load mat file; 0:direct calculate
%gpuDevice(1)%select GPU device
constantfile;
clear gam
rk4=1;%1:rk4,0:heun Method,2:4th predictor-corrector
bc=1;%0.periodic condition;1,not periodic
dipolee=0;%enable dipole?
DMIenable=0;
dwcalc=0;%1:simulate dw motion 0: no domain wall
thermalenable=0;%enable thermal field?
%% system generation
natomx=10;natomy=10;%no. of cells along x,y direction
compositionn=0.1;%composition percentage (X) of RE element, e.g. GdX(FeCo)1-X
d=0.4e-9;%[m],lattice constant
natom=natomx*natomy;
if ddebug
    load('ddebug.mat');
end
systemgeneration();
%% FiM parameters
Ksim=0.4e-3*ele;%[J], easy-axis anisotropy
Jgdgd=-1.26e-21;Jfefe=-2.835e-21;Jfegd=1.09e-21;%[J/link][1]
gTM=2.2;gRE=2;%g-factor
gamTM=gTM*mub/(hbar*ele);%1/(s.T)refer to "PRL 97, 217202 (2006), Jiang Xin"
gamRE=gRE*mub/(hbar*ele);
tz=d;%[m],thickness of FiM
musRE=7.63*mub;musTM=2.217*mub;%[J/T]magnetic moment [1]
msRE=musRE/d^3;%[A/m], saturation magnetization
msTM=musTM/d^3;
if DMIenable
    Dsim=128e-6*ele;%[J], DMI
else
    Dsim=0;
end
alp=0.02;%Gilbert damping
%% electrical parameters
jcSOT=1e9;%[A/m2]
jcSTT=1.5e9;%[A/m2]
Hext=[0,0,0e-3];% corresponding to runtime2
%% SOT parameters
SOT_DLT=1;%1(0),enable(disable) SOT damping torque
SOT_FLT=0;%1(0),enable(disable) SOT field-like torque
psjSHE=[0,1,0];%spin flux polarization
psjSHEx=psjSHE(1);
psjSHEy=psjSHE(2);
psjSHEz=psjSHE(3);
thetaSH=0.2;%spin hall angle
if SOT_FLT
    chi=0;%ratio of FLT/DLT
else
    chi=0;
end
BDSOTRE=hbar/2*thetaSH*jcSOT/(msRE*tz);%[T]
BDSOTTM=hbar/2*thetaSH*jcSOT/(msTM*tz);
%% STT parameters
STT_DLT=1;%1(0),enable(disable) SOT damping torque
STT_FLT=0;%1(0),enable(disable) SOT field-like torque
psjSTT=[0,0,1];%spin flux polarization
psjSTTx=psjSTT(1);
psjSTTy=psjSTT(2);
psjSTTz=psjSTT(3);
etaSTT=0.8;%spin hall angle
if STT_FLT
    chiSTT=0;%ratio of FLT/DLT
else
    chiSTT=0;
end
BDSTTRE=hbar/2*etaSTT*jcSTT/(msRE*tz);%[T]
BDSTTTM=hbar/2*etaSTT*jcSTT/(msTM*tz);
%% other parameters
T=100;%[K]
%% time control
gpusave=2e-12;%how often saving gpu data
tstep=2e-15;
gpusteps=round(gpusave/tstep);
runtime=3*gpusave;%second run for dw motion
savetstep=100;%to reduce data size
totstep=round(runtime/tstep);
t=linspace(tstep,runtime,totstep);

if (SOT_DLT || SOT_FLT) && ~(rk4==1)
    error('only rk4 is implemented for spin torque driven')
end

if loadstartm
    clear mx_init my_init mz_init atomtype_
    load('startm_x0.1_10x10.mat')
    if natomx~=natomxcheck || natomx~=natomxcheck || compositionn~=compositionncheck
       error('system not consistent') 
    end
    clear natomxcheck natomycheck compositionncheck
    mx_init=mmxstart;
    my_init=mmystart;
    mz_init=mmzstart;
end
if(0)%view initial state
    gridx = 1:natomx;
    gridy = 1:natomy;
    [gridplotx,gridploty] = meshgrid(gridx,gridy);
    gridz=zeros(natomx,natomy);
    figure;hold on%initial magnetization
    for cty=1:natomy
        for ctx=1:natomx
            if atomtype_(ctx,cty)==0%TM
                quiver3(gridplotx(natomx-ctx+1,cty),gridploty(natomx-ctx+1,cty),gridz(ctx,cty),...
                    mx_init(ctx,cty),my_init(ctx,cty),mz_init(ctx,cty),'r');
            else
                quiver3(gridplotx(natomx-ctx+1,cty),gridploty(natomx-ctx+1,cty),gridz(ctx,cty),...
                    mx_init(ctx,cty),my_init(ctx,cty),mz_init(ctx,cty),'b');
            end
        end
    end
end
%% dynamic calc
integrate_llg(); toc
%% save data
save('final.mat')






