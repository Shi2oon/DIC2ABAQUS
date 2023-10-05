% A code to prepare HR-EBSD and strain data to be inputed to Abaqus. the data is
% fomrated as receieved from the xEBSD v3 software.
% restoredefaultpath; warning off;       
close all;     addpath(genpath(pwd));              
clc;clear;     set(0,'defaultAxesFontSize',20);       set(0,'DefaultLineMarkerSize',12)   
DS = com.mathworks.mde.desk.MLDesktop.getInstance();  DS.closeGroup('Variables');  
tic;
%% Setting the scene
Dir.fullpath     = pwd; % filr Directory
Dir.fillname     = '2D_11I_6II WS\mm_DISP';       % DONT include extension,in .mat or .dat format
Dir.input_unit   = 'mm';        % meter (m) or milmeter (mm) or micrometer(um);
Dir.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
Dir.Operation    = 'DIC';       % calculation mode\n Str = Strain, xED = xEBSD, DIC for Displacement
Dir.unique       = 'Synthetic_Data';   % uniqe name for data set, can be = Dir.fillname;
Dir.stressstat   = 'plane_stress'; % 'plane_strain

%% INPUT MATERIAL PROPERTIES AND DATA
% Poisson's ratio,          Young's Modulus [Pa],      		Material Name     
  Dir.nu    = 0.3;          Dir.E  = 210E9;                 Dir.Mat = 'Ferrite';
% 'E' for Elastic or 'R' for Ramberg-Osgood or 'A' for Elastic-Anisotropic
% or 'P' for elasto-plastic
  Dir.type  = 'E'; 
% if 'Ramberg-Osgood' type of material input                Yield Stress [Pa] 
%   Dir.Exponent = 26.67;     Dir.Yield_offset = 1.24;        Dir.yield = 4E9;
% if 'Elastic-Anisotropic' you need to define the stifness tensor 
%     Dir.Stiffness = NaN(6,6);
% if 'Plastic' you need to define the stifness tensor 
% Dir.type  = 'P'; 
% Dir.Plastic_Strain = [0.0010,0.0023,0.0038,0.0059,0.0085,0.0120,0.0150,0.0205,0.0260,0.0336]'-0.0010;
% Dir.Plastic_Stress = [150,350,547,849,1156,1569,1934,2616,3142,3912]'*1e6;
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of USER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Locate the crack 
% [alldata,Dir] = ThingsWentWrong(Dir,Maps);
[DATA,UnitOffset,Dir, msk,SaveD] = Locate_Crack(0,Dir.input_unit,[Dir.fullpath '\' Dir.fillname],Dir); 
    
%% prepare and run abaqus cae
[Dir.Abaqus,Abaqus.CAE] = PrintRunCode(Dir, ...
        msk,SaveD,ceil(min(size(DATA.X1))*0.5-2),UnitOffset);

%% Post Processing
sizE = diff(unique(DATA.X1(:)));
[Abaqus.J,Abaqus.Keff,Abaqus.KI,Abaqus.KII] = PlotKorJ(Dir.Abaqus,Dir.E,...
    UnitOffset,Dir.input_unit,sizE(2));
% Abaqus.Stress_Maps = readAbaqusOutput(Dir.Abaqus,Dir.unique);
% imagesc(Abaqus.Stress_Maps.X(1,:),Abaqus.Stress_Maps.Y(:,1),Abaqus.Stress_Maps.S12); axis image
save([fileparts(Dir.Abaqus) '\AbaqusOutput.mat'],'Abaqus','Dir');

fprintf('All processing Took %.2f hours\n',toc/3600);
