% 2D atomistic model for FiM
% require nvidia GPU
clearvars -except astro_output_file;clc;close all;astro_main_timer=tic;
if ~exist('astro_output_file','var') || isempty(astro_output_file)
    astro_output_file='final.mat';
end
astro_output_file=char(astro_output_file);
%% control parameter
constantfile;
clear gam
astro_defaults=astro_default_config();
astro_config=astro_defaults;
astro_default_fields=fieldnames(astro_defaults);
astro_conditional_fields={'natom_mc_W','natom_mc_L','fixededgeL', ...
    'mxleft','myleft','mzleft','mxright','myright','mzright'};
for astro_default_index=1:numel(astro_default_fields)
    astro_default_name=astro_default_fields{astro_default_index};
    if ~ismember(astro_default_name,astro_conditional_fields)
        eval([astro_default_name '=astro_defaults.(astro_default_name);']);
    end
end
%% use fixed atom distribution
if save_fixed_atom_distrib && load_fixed_atom_distrib
   error('only one of load_fixed_atom_distrib or save_fixed_atom_distrib can be enabled'); 
end
if load_fixed_atom_distrib
    display('use fixed atom distribution')
elseif save_fixed_atom_distrib
    display('save fixed atom distribution')
else
   display('use random atom distribution')
   rng('shuffle');
end
%% enable dipole
%0:not calcualte dipole, 1: direct calculation using CPU
%2:direct calculation using CPU, 3:calculation using macrocell
if dipolemode==3
    natom_mc_W=astro_defaults.natom_mc_W;%number of atoms in the macrocell along x direction
    natom_mc_L=astro_defaults.natom_mc_L;%number of atoms in the macrocell along y direction
end
%% fix atoms at the edge
if enablefixedge
    fixededgeL=astro_defaults.fixededgeL;%No of atoms fixed along Length
    mxleft=astro_defaults.mxleft;
    myleft=astro_defaults.myleft;
    mzleft=astro_defaults.mzleft;
    mxright=astro_defaults.mxright;
    myright=astro_defaults.myright;
    mzright=astro_defaults.mzright;
end
clear astro_defaults astro_default_fields astro_default_index astro_default_name astro_conditional_fields

%% optional control
%gpuDevice(1)%select GPU device
%% system generation
%note this is different to the x,y,z in h_ex or hdmi etc
natom=natomW*natomL;
systemgeneration();
%% time control
gpusteps=round(gpusave/tstep);
runtime=gpurun_number*gpusave;%second run for dw motion
totstep=round(runtime/tstep);
t=(0:savetstep:totstep)*tstep;

if ~mod(gpusteps,savetstep)==0
    error('gpusteps should be multiple integer times of savetstep, otherwise there might be errors')
end

final_m_savestep=numel(t);

if rk4~=1
    error('ASTRO:UnsupportedSolver', ...
        ['Unsupported ASTRO solver mode rk4=%g. Only rk4=1 (RK4) is ' ...
        'supported. rk4=0 (Heun) and rk4=2 (predictor-corrector) are ' ...
        'disabled because they are not validated production paths.'], rk4);
end

if(0)%view initial state
    addpath('C:\Users\zzf-m\OneDrive\code_softwares\general\gitcontrol\3D_vector_plot')
    mmx_show=zeros(natomW,natomL,2);%initial magnetization
    mmy_show=zeros(natomW,natomL,2);
    mmz_show=zeros(natomW,natomL,2);
    mmx_show(:,:,1)=mx_init;
    mmy_show(:,:,1)=my_init;
    mmz_show(:,:,1)=mz_init;
    plottime=0e-12;% the time you want to see the magnetization state
    plotWstep=1;%grid size along Width direction
    plotLstep=1;%grid size along length direction
    scale3d=0.5;%scale factor of the arrow size
    arrowwidth=10;
    plotmode=0;
    %0:plot for 2D atomistic model;1:cross-section plot for 3D atomistic
    %model;2:plot for 3D atomistic model
    colorbarplot=0;%for FiM, choose 0 so that Fe and Gd are represented by red
    % and blue, respectively. For FM, choose 1 so that the colorbar indicates
    % the magnitude of mz

    % movie options
    generatemovie=0;%1:generate movie, 0: no movie
    movieduration=10;%seconds
    plotmoviestep=10;%movie step

    threedplott(mmx_show,mmy_show,mmz_show,runtime,atomtype_,plottime,plotWstep,plotLstep,...
        scale3d,arrowwidth,plotmode,colorbarplot,generatemovie,movieduration,plotmoviestep)
end
%% dynamic calc
integrate_llg();
astro_elapsed_seconds=toc(astro_main_timer);
%% save data
astro_config.natom=natom;
astro_config.gpusteps=gpusteps;
astro_config.runtime=runtime;
astro_config.totstep=totstep;
astro_config.final_m_savestep=final_m_savestep;
astro_result=astro_build_production_result(astro_config, atomtype_, ...
    mx_init, my_init, mz_init, mmx, mmy, mmz, t, ...
    astro_elapsed_seconds, astro_output_file);
astro_output_dir=fileparts(astro_output_file);
if ~isempty(astro_output_dir) && ~exist(astro_output_dir,'dir')
    mkdir(astro_output_dir);
end
save(astro_output_file,'astro_result')
fprintf('ASTRO production result saved to %s\n', astro_output_file);
