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
% 'E' for Elastic material
% Poisson's ratio,          Young's Modulus [Pa],      		    
  Maps.nu    = 0.3;         Maps.E  = 220E9;            Maps.type  = 'E'; 
% if 'Ramberg-Osgood' type of material input            Yield Stress [Pa] 
% Maps.Exponent = 26.67;      Maps.Yield_offset = 1.24; Maps.yield = 4E9;
% Maps.type  = 'R';           
% Maps.E  = 193E9;            Maps.nu    = 0.3;
% if 'Elastic-Anisotropic' you need to define the stifness tensor 
% Maps.type  = 'A'; 
% Maps.Stiffness = [  283  121  121   0   0   0
%                     121  283  121   0   0   0
%                     121  121  283   0   0   0
%                     0    0    0     81  0   0
%                     0    0    0     0   81  0
%                     0    0    0     0   0   81]*1e9;
% Maps.nu    = 0.00;         Maps.E  = 000E9; for non-cubic materials

%% %%%%%%%%%%%%%%%%%%%%%%% END of USER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%
Data = importdata(DataDirect);
[~,RawData ] = reshapeData(Data.data);

% to remove rigid body movement for stereo DIC
% [theta,RawData.Ux,RawData.Uy,RawData.Uz,rotCentre] = ...
%     RotRemoving('true',Maps.X,Maps.Y,RawData.Ux,RawData.Uy,RawData.Uz);
% to remove rigid body movement for 2D DIC
% [theta,RawData.Ux,RawData.Uy,rotCentre] = ...
%     RotRemoving('true',Maps.X,Maps.Y,RawData.Ux,RawData.Uy);

Maps.X  = RawData.X1;       Maps.Y = RawData.Y1;   
Maps.Ux = RawData.Ux;       Maps.Uy = RawData.Uy;
Maps.Z = RawData.Z1;      Maps.Uz = RawData.Uz;

% [data] = Cropping10(Maps.X,Maps.Y,Maps.Ux, Maps.Uy);
% Maps.X = data.X;       Maps.Y = data.Y;
% Maps.Ux = data.Ux;     Maps.Uy = data.Uy;
% Maps.Uy(Maps.Uy==0)=NaN;    Maps.Ux(Maps.Ux==0)=NaN;
% Maps.Uy = inpaint_nans(Maps.Uy);
% Maps.Ux = inpaint_nans(Maps.Ux);

[J,KI,KII,KIII] = DIC2ABAQUS(Maps);