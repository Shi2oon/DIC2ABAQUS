% Validation from sythetic data
clc;clear;close all
addpath(genpath([pwd '\functions']));
% put the KI, KII, KII values for stereo-DIC data as in
% Calibration_2DKIII(KI,KII,KIII);for 2D, put KIII as zerp
[~,~,DIC_Data] = Calibration_2DKIII(5,1,3); % creating the synthatic field

Maps.Mat          = 'Calibration';
Maps.type         = 'E';% 'A' if u want to use anistropic matrix or 'E' for linear elastic
Maps.input_unit   = 'um';        % meter (m) or milmeter (mm) or micrometer(um); 
Maps.pixel_size   = 1;           % if DIC values are in pixel, 1 if in physical units;
Maps.Operation    = 'DIC';       % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = 'Calibration';

[J,KI,KII,KIII] = DIC2CAE(Maps);