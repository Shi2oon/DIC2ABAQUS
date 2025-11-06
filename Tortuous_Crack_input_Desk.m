% This the input desk for cruved or non striaght cracks where you need to
% exclude the crack geomtery first and save the displacement field with the
% excluded crack geomtry in a format you can use with this code, preferably
% as a .dat file,

clc;clear;close all
addpath(genpath([pwd '\Non-unifrom maps']));
addpath(genpath([pwd '\Miscellaneous']));
resultsDir = fullfile(pwd, 'Miscellaneous','Tortuous_Crack_Data.dat'); % file location
Maps.results = erase(resultsDir,'.dat');
Maps.input_unit     = 'um';         % meter (m) or milmeter (mm) or micrometer(um);
Maps.pixel_size   = 1e-3;           % if DIC values are in pixel, 1 if in physical units;
Maps.Operation    = 'DIC';          % Strain, xED = xEBSD, DIC = Displacement
Maps.stressstat   = 'plane_stress'; % 'plane_stress' OR 'plane_strain'
Maps.unique       = 'Crack_in_Al_5052';

%% INPUT MATERIAL PROPERTIES AND DATA
% 'E' for Elastic material
Maps.Mat = 'Al_5052';       %Material Name, do not use space
% Poisson's ratio,          Young's Modulus [Pa],     
Maps.nu    = 0.321;         Maps.E  = 70E9;     
Maps.type  = 'E';
%{
% if 'Ramberg-Osgood' type of material input
% The hardening exponet for Ramberg-Osgood is bigger than 1
% see https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/stm/default.htm?startat=ch04s03ath111.html
%                                                   Yield Stress [Pa]
Maps.Exponent = 26.67;  Maps.Yield_offset = 1.24;   Maps.yield = 4E9;
Maps.type  = 'R';       Maps.E  = 210E9;            Maps.nu    = 0.3;


% if 'Elastic-Anisotropic' you need to define the stifness tensor
Maps.type  = 'A';
Maps.Stiffness = [  283  121  121   0   0   0
                    121  283  121   0   0   0
                    121  121  283   0   0   0
                    0    0    0     81  0   0
                    0    0    0     0   81  0
                    0    0    0     0   0   81]*1e9;
Maps.nu    = 0.30;         
Maps.E  = 210E9; % for non-cubic materials
% Place EBSD orientation data in the same directory as this script
% Data must be of the EBSDsquare type, in a .mat file
% This script is compatibly with MTEX-6.0.0
Maps.EBSDfilename = 'EBSDsquareExample.mat'; 


% If using UMAT set the material properties in usermaterials.f
% See OXFORD-UMAT documentation (https://github.com/TarletonGroup/CrystalPlasticity)
Maps.type  = 'U';
Maps.depvar = 50;
Maps.materialID = 3; % (1) bcc, (2) fcc, (3) hcp
Maps.PROPS = 0;
Maps.UMATfilepath = fullfile(pwd, 'OXFORD-UMAT v3.1');
% Place EBSD orientation data in the same directory as this script
% Data must be of the EBSDsquare type, in a .mat file
% This script is compatible with MTEX-6.0.0
Maps.EBSDfilename = 'EBSDsquareExample.mat'; 
Maps.modelThickness = 3e-3; % in microns
Maps.zElems = 3; % number of elements through thickness
%}



tic
% DIC2CAE_wNAN(MatP, Crack, resultsDir,angle_deg)
% here "Crack" is the crack location starting with x corrdinates and then y
% you can leave it empty and then will be asked to provide it
% you can also enter the crack angle or where you think the crack will go
% using "angle_deg" variable which importnat when calculating the
% J-integral as it depended on the virtual crack extension method, see
% https://doc.comsol.com/6.2/doc/com.comsol.help.sme/sme_ug_theory.06.093.html
[BCf, Maps.UnitOffset,Maps.stepsize] = DIC2ABAQUS_wNAN(Maps, ...
                    [1375 955]*Maps.pixel_size,resultsDir,180);
ABAQUS = PrintRunCode(Maps,Maps.results);
fprintf('Execution time: %.1f minutes.\n', toc/60);
[J,K,KI,KII,Direction] = PlotKorJ(ABAQUS,Maps.E,Maps.UnitOffset);
plotJKIII(KI,KII,[],J,Maps.stepsize,Maps.input_unit)
saveas(gcf, [Maps.results '\' Maps.unique '_J_KI_II.fig']);
saveas(gcf, [Maps.results '\' Maps.unique '_J_KI_II.tif']);    close all
save([Maps.results '\' Maps.unique '_DIC2CAE.mat'],'Maps','J','KI','KII',...
    'K','Direction');

%% changing the crack direction
load([Maps.results '_DIC2CAE.mat'],'MatP','Direction');
fprintf('\nRecommended virtual crack extension for max. J-integral direction is %.1f ± %.1f\t\n',...
    round(Direction.true,1), round(Direction.div,1))
Ans = questdlg_timer(60,['Virtual crack extension is ' ...
    num2str(round(Direction.true,1)) ' ± ' num2str(round(Direction.div,1)) ...
    ', Do you want to adjust?'], ...
    [ num2str(round(Direction.true,1)) ' ± ' num2str(round(Direction.div,1))],...
    'Y','N','N');
if Ans == 'Y'
    ABAQUS = fullfile(Maps.results,['Job-' Maps.unique]);
    [ABAQUS] = Adjust4Direction(ABAQUS,round(Direction.true,1));
end
[J,K,KI,KII,Direction] = PlotKorJ(ABAQUS,Maps.E,UnitOffset);
plotJKIII(KI,KII,[],J,Maps.stepsize,Maps.input_unit)
saveas(gcf, [Maps.results '\' Maps.unique '_corrected_J_KI_II.fig']);
saveas(gcf, [Maps.results '\' Maps.unique '_corrected_J_KI_II.tif']); close all
save([Maps.results '\' Maps.unique '_DIC2CAE_corrected.mat'],'Maps',...
    'J','KI','KII','K','Direction');
