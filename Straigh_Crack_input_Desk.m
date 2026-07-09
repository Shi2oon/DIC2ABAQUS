% This the input desk for straight cracks where you do not need to exclude
% the crack geomtery yourself. Here that will be done automatically.

% To genereat syntahtic data to validate this the code use
% "InputDesk_Validate".
clc;clear;close all
addpath(genpath([pwd '\Functions']));
addpath(genpath([pwd '\Miscellaneous']));
DataDirect = fullfile(pwd,'\Miscellaneous\1KI-2KII-3KIII_Data.dat'); % file location
Maps.results = erase(DataDirect,'.dat');
% Domain size (square, crack tip at centre).
Maps.input_unit     = 'mm';            % meter (m) or milmeter (mm) or micrometer(um);
Maps.pixel_size   = 1;              % if DIC values are in pixel, 1 if in physical units;
Maps.Operation    = 'DIC';          % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = 'Calibration';
Maps.Mat 	  = 'Ferrite'; %Material Name 

%% INPUT MATERIAL PROPERTIES AND DATA       
% Select material type:
%   'E' = Elastic
%   'R' = Ramberg–Osgood
%   'A' = Elastic–Anisotropic
%   'U' = UMAT (User Material)
Maps.type = 'E';  % <<<--- CHANGE THIS TO SELECT YOUR MATERIAL MODEL


switch Maps.type
    % --------------------------------------------------------------------
    case 'E'  % Elastic Material
        Maps.Mat = 'Ferrite';       % Material Name (no spaces)
        Maps.E  = 220E9;             % Young's Modulus [Pa]
        Maps.nu = 0.3;            % Poisson's ratio

    % --------------------------------------------------------------------
    case 'R'  % Ramberg–Osgood Material
        % The hardening exponent must be > 1
        % Reference: https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/stm/default.htm?startat=ch04s03ath111.html
        Maps.Mat = 'Ferrite';
        Maps.E = 220e9;             % Young's Modulus [Pa]
        Maps.nu = 0.3;              % Poisson's ratio
        Maps.Exponent = 26.67;      % Hardening exponent
        Maps.Yield_offset = 1.24;   % Offset for yield definition
        Maps.yield = 4e9;           % Yield stress [Pa]

    % --------------------------------------------------------------------
    case 'A'  % Elastic–Anisotropic Material
        Maps.Mat = 'Anisotropic_Ferrite';
        Maps.Stiffness = [283 121 121 0 0 0;
                          121 283 121 0 0 0;
                          121 121 283 0 0 0;
                          0   0   0  81 0 0;
                          0   0   0   0 81 0;
                          0   0   0   0  0 81] * 1e9;  % [Pa]
        Maps.nu = 0.30;
        Maps.E  = 210e9; % Reference value for anisotropic case

        % EBSD data required for this model
        Maps.EBSDfilename = 'EBSDsquareExample.mat';  % must be in same directory as this script

    % --------------------------------------------------------------------
    case 'U'  % UMAT (User Material)
        % See OXFORD-UMAT documentation: https://github.com/TarletonGroup/CrystalPlasticity
        % Material properties defined in usermaterials.f
        Maps.Mat = 'UserDefined_Ferrite';
        Maps.depvar = 50;           % Number of user-defined state variables
        Maps.materialID = 1;        % (1) bcc, (2) fcc, (3) hcp
        Maps.PROPS = 0;             % NPROPS value for OXFORD-UMAT

        % File path of main UMAT. Necessary environment files should be included in the same folder
        Maps.UMATfilepath = fullfile(pwd, 'OXFORD-UMAT\OXFORD-UMAT v3.3\OXFORD-UMAT.f');

        % Full path to abaqus command shortcut. Abaqus should be linked with Fortran compiler        
        Maps.abqCmdShortcutPath = ...
            'C:\ProgramData\Microsoft\Windows\Start Menu\Programs\Dassault Systemes SIMULIA Abaqus CAE 2017\Abaqus Command.lnk'; 

        % EBSD data required for this model
        Maps.EBSDfilename = 'EBSDsquareExample.mat';  % must be in same directory as this script
        Maps.E  = 210e9; % Placeholder value

    % --------------------------------------------------------------------
    otherwise
        error('Unknown Maps.type selected. Use ''E'', ''R'', ''A'', or ''U''.');
end
%% %%%%%%%%%%%%%%%%%%%%%%% END of USER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%
Data = importdata(DataDirect);
[~,RawData ] = reshapeData(Data.data*Maps.pixel_size);

% to remove rigid body movement for stereo DIC
% [theta,RawData.Ux,RawData.Uy,RawData.Uz,rotCentre] = ...
%     RotRemoving('true',Maps.X,Maps.Y,RawData.Ux,RawData.Uy,RawData.Uz);
% to remove rigid body movement for 2D DIC
% [theta,RawData.Ux,RawData.Uy,rotCentre] = ...
%     RotRemoving('true',Maps.X,Maps.Y,RawData.Ux,RawData.Uy);

Maps.X  = RawData.X1;       Maps.Y = RawData.Y1;   
Maps.Ux = RawData.Ux;       Maps.Uy = RawData.Uy;
try %stereo
	Maps.Z = RawData.Z1;      Maps.Uz = RawData.Uz;
catch %2D
	.Z = RawData.X1*0;      Maps.Uz = RawData.Uy*0;
end

% [data] = Cropping10(Maps.X,Maps.Y,Maps.Ux, Maps.Uy);
% Maps.X = data.X;       Maps.Y = data.Y;
% Maps.Ux = data.Ux;     Maps.Uy = data.Uy;
% Maps.Uy(Maps.Uy==0)=NaN;    Maps.Ux(Maps.Ux==0)=NaN;
% Maps.Uy = inpaint_nans(Maps.Uy);
% Maps.Ux = inpaint_nans(Maps.Ux);

[J,KI,KII,KIII] = DIC2ABAQUS(Maps);
