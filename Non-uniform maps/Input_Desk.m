clc;clear;close all
addpath(genpath([pwd '\functions']));
pwdfor = 'C:\Users\ak13\Documents\Fo modelling';
resultsDir = fullfile(pwdfor,'B0045_Data.dat'); % file location
MatP.results = erase(resultsDir,'.dat');
% Domain size (square, crack tip at centre).
MatP.input_unit     = 'mm';            % meter (m) or milmeter (mm) or micrometer(um);
MatP.pixel_size   = 1;              % if DIC values are in pixel, 1 if in physical units;
MatP.Operation    = 'DIC';          % Strain, xED = xEBSD, DIC = Displacement
MatP.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
[~,MatP.unique ]     = fileparts(MatP.results);

%% INPUT MATERIAL PROPERTIES AND DATA
% 'E' for Elastic material
% Poisson's ratio,          Young's Modulus [Pa],      		Material Name     
  MatP.nu    = 0.321;         MatP.E  = 70E9;                MatP.Mat = 'Al_5052';
  MatP.type  = 'E'; 
% if 'Ramberg-Osgood' type of material input                Yield Stress [Pa] 
% MatP.Exponent = 0.2;      MatP.Yield_offset = 0.2/100;      MatP.yield = 193E6;
% MatP.type  = 'R';         MatP.Mat = 'Al_5052';
% MatP.E  = 70E9;           MatP.nu  = 0.321; % Poisson's Ratio
% if 'Elastic-Anisotropic' you need to define the stifness tensor 
% MatP.type  = 'A'; 
% MatP.Stiffness = [  283  121  121   0   0   0
%                     121  283  121   0   0   0
%                     121  121  283   0   0   0
%                     0    0    0     81  0   0
%                     0    0    0     0   81  0
%                     0    0    0     0   0   81]*1e9;
% MatP.Mat = 'Ferrite';

%% Crack an direction

[BCf, UnitOffset] = DIC2CAE_wNAN(MatP, [], resultsDir);