function config = astro_default_config()
%ASTRO_DEFAULT_CONFIG Production defaults shared by main and benchmarks.

constantfile;
clear gam

config = struct();

config.rk4 = 1; % 1: RK4, 0: Heun, 2: 4th predictor-corrector
config.bc = 1; % 0: periodic condition; 1: not periodic
config.DMIenable = 1;
config.dwcalc = 1; % 1: simulate DW motion; 0: no domain wall
config.thermalenable = 0; % enable thermal field
config.loadstartm = 0; % 1: load magnetization file; 0: direct calculate
config.startmname = 'startm_20x250.mat';
config.load_fixed_atom_distrib = 0; % load fixed atom distribution
config.save_fixed_atom_distrib = 0; % save fixed atom distribution
config.dipolemode = 0; % 0: disabled; 1/2: direct; 3: macrocell
config.natom_mc_W = 2; % atoms in macrocell along width
config.natom_mc_L = 2; % atoms in macrocell along length
config.enablefixedge = 0; % fix atoms at the edge
config.fixededgeL = 3; % atoms fixed along length
config.mxleft = 0;
config.myleft = 0;
config.mzleft = 1;
config.mxright = 0;
config.myright = 0;
config.mzright = -1;

config.natomW = 20; % cells along width
config.natomL = 30; % cells along length
config.compositionn = 0.1; % RE/Gd composition fraction
config.d = 0.4e-9; % [m], lattice constant

config.Ksim = 0.4e-3*ele; % [J], easy-axis anisotropy per atom
config.Jgdgd = -1.26e-21; % [J/link]
config.Jfefe = -2.835e-21; % [J/link]
config.Jfegd = 1.09e-21; % [J/link]
config.gTM = 2.2; % TM g-factor
config.gRE = 2; % RE g-factor
config.gamTM = config.gTM*mub/(hbar*ele); % [1/(s.T)]
config.gamRE = config.gRE*mub/(hbar*ele); % [1/(s.T)]
config.tz = config.d; % [m], FiM thickness
config.musRE = 7.63*mub; % [J/T], RE magnetic moment per atom
config.musTM = 2.217*mub; % [J/T], TM magnetic moment per atom
config.msRE = config.musRE/config.d^3; % [A/m], RE saturation magnetization
config.msTM = config.musTM/config.d^3; % [A/m], TM saturation magnetization
if config.DMIenable
    config.Dsim = 128e-6*ele; % [J], DMI
else
    config.Dsim = 0;
end
config.alp = 0.02; % Gilbert damping

%% electrical parameters
config.jcSOT = 1e9; % [A/m2]
config.jcSTT = 1.5e9; % [A/m2]
config.Hext = [0, 0, 0e-3]; % [T], external field

%% SOT parameters
config.SOT_DLT = 1; % enable SOT damping-like torque
config.SOT_FLT = 0; % disable SOT field-like torque
config.psjSHE = [0, 1, 0]; % spin flux polarization
config.psjSHEx = config.psjSHE(1);
config.psjSHEy = config.psjSHE(2);
config.psjSHEz = config.psjSHE(3);
config.thetaSH = 0.2; % spin Hall angle
config.chi = 0; % FLT/DLT ratio
config.BDSOTRE = config.SOT_DLT*hbar/2*config.thetaSH* ...
    config.jcSOT/(config.msRE*config.tz); % [T]
config.BDSOTTM = config.SOT_DLT*hbar/2*config.thetaSH* ...
    config.jcSOT/(config.msTM*config.tz); % [T]

%% STT parameters
config.STT_DLT = 0; % disable STT damping-like torque
config.STT_FLT = 0; % disable STT field-like torque
config.psjSTT = [0, 0, 1]; % spin flux polarization
config.psjSTTx = config.psjSTT(1);
config.psjSTTy = config.psjSTT(2);
config.psjSTTz = config.psjSTT(3);
config.etaSTT = 0.8; % spin Hall angle
config.chiSTT = 0; % FLT/DLT ratio
config.BDSTTRE = config.STT_DLT*hbar/2*config.etaSTT* ...
    config.jcSTT/(config.msRE*config.tz); % [T]
config.BDSTTTM = config.STT_DLT*hbar/2*config.etaSTT* ...
    config.jcSTT/(config.msTM*config.tz); % [T]

%% other parameters
config.T = 100; % [K]

config.gpusave = 1e-12; % [s], GPU run chunk duration
config.gpurun_number = 2;
config.tstep = 2e-16; % [s]
config.savetstep = 100; % integration steps between saved states
end
