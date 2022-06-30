function [files,Operation,unit,stif] = ...
    Westergaard_Validation(StressIntensityFactor,state,Mode,stPs,file2)
	
	% a function to create a synthtic displacement data around a crack tip 
	% in steel for a specific mode I, II or both stress intesity factors 
	% which can be used with the DIC2ABAQUS code
	
	% Varialbles
		% StressIntensityFactor: in MPa
		% state: either 'plane_strain' or 'plane_stress'
		% Mode: 'I' or 'II' or 'fun' for mixed mode entered as matrix [I, II]
   
if ~exist('Mode','var');        Mode = 'I';         end
if isempty(Mode);               Mode = 'I';         end

% close all; clear; clc
if exist('file2','var')
    newfile = fullfile(file2, [num2str(StressIntensityFactor) '_' Mode ' WS']);
elseif strcmpi(Mode, 'I') || strcmpi(Mode, 'II')
    newfile = [ pwd '\2D_' num2str(StressIntensityFactor) '_' Mode ' WS']; 
elseif strcmpi(Mode, 'fun')
    newfile = [ pwd '\2D_' num2str(StressIntensityFactor(1)) 'I_' num2str(StressIntensityFactor(2)) 'II WS']; 
end
mkdir(newfile);

minGrid = -4E-3; % [m]
if exist('stPs','var')
    gridStep = 8e-3/stPs;
else
    gridStep = 0.2E-3; % [m]
end

maxGrid = 4E-3; % [m]

xvec = minGrid : gridStep : maxGrid; % [m]
yvec = minGrid : gridStep : maxGrid; % [m]

[x,y] = meshgrid(xvec,yvec); % [m]
% StressIntensityFactor = 30; % [MPa m^0.5]
fprintf('preparing synthetic Westergaard Solution Data .. ');
E = 210E9; % Young's Modulus [Pa]
nu = 0.3; % poisson ratio
mu = E/(2.*(1+nu)); % Shear Modulus [Pa]

    switch state
        case 'plane_strain'
            kappa = 3 - (4 .* nu); % [/]
        case 'plane_stress'
            kappa = (3 - nu)./(1 + nu); % [/]
    end
[th,r] = cart2pol(x,y);

switch Mode
    case 'I'
        KI = StressIntensityFactor * 1E6; KII=0;

    case 'II'
        KI = 0; KII = StressIntensityFactor * 1E6;
    case 'fun'
        KI = StressIntensityFactor(1) * 1E6;
        KII = StressIntensityFactor(2) * 1E6;
end

% calulate the displacement % Anderson p99 A2.44a
ux = ( 0.5*KI/mu*sqrt(r/(2*pi)).*(+cos(th/2).*(kappa-cos(th)))+...
              0.5*KII/mu*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa+2+cos(th))));
uy = ( 0.5*KI/mu*sqrt(r/(2*pi)).*(+sin(th/2).*(kappa-cos(th)))+...
              0.5*KII/mu*sqrt(r/(2*pi)).*(-cos(th/2).*(kappa-2+cos(th))));


subplot(1,3,1); contourf(x,y,ux); 
axis image; box off; c=colorbar;    c.Label.String= ['U_{x} [m]'];
subplot(1,3,2);contourf(x,y,uy); 
axis image; box off; c=colorbar;    c.Label.String= ['U_{y} [m]'];
subplot(1,3,3);contourf(x,y,sqrt(uy.^2+ux.^2)); 
axis image; box off; c=colorbar;    c.Label.String= ['U_{mag} [m]'];
set(gcf,'Position',[400 84 1209 1026]);
saveas(gcf, [newfile '\' Mode '_Disp_fields.tiff']);   
saveas(gcf, [newfile '\' Mode '_Disp_fields.fig']); close  
[dux_dx,dux_dy] = gradient(ux,gridStep);
[duy_dx,duy_dy] = gradient(uy,gridStep);
exx = dux_dx;
eyy = duy_dy;
exy = 0.5*(dux_dy + duy_dx);

    %% save displacement fields
    Operation{1,1} = 'DIC';               unit{1,1}  = 'm';
    alldata = [x(:) y(:) ux(:) uy(:)]; % [m]
    if exist('stPs','var')
        files{1,1} = [newfile '\m_DISP_' num2str(stPs) '.mat'];
    else
        files{1,1} = [newfile '\m_DISP.mat'];
    end
    save(files{1,1},'alldata')
    
    Operation{1,2} = 'DIC';               unit{1,2}  = 'mm';
    alldata = [x(:) y(:) ux(:) uy(:)].*1e3; % [m]
    if exist('stPs','var')
        files{1,2} = [newfile '\mm_DISP_' num2str(stPs) '.mat'];
    else
        files{1,2} = [newfile '\mm_DISP.mat'];
    end
    save(files{1,2},'alldata')
    
    Operation{1,3} = 'DIC';               unit{1,3}  = 'um';
    alldata = [x(:) y(:) ux(:) uy(:)].*1e6; % [m]
    if exist('stPs','var')
        files{1,3} = [newfile '\um_DISP_' num2str(stPs) '.mat'];
    else
        files{1,3} = [newfile '\um_DISP.mat'];
    end
    save(files{1,3},'alldata')

    %% save strain fields
    Operation{2,1} = 'Str';               unit{2,1}  = 'm';
    alldata = [x(:) y(:) exx(:) eyy(:) exy(:)]; % [m]
    if exist('stPs','var')
        files{2,1} = [newfile '\m_Strain_' num2str(stPs) '.mat'];
    else
        files{2,1} = [newfile '\m_Strain.mat'];
    end
    save(files{2,1},'alldata')
    
    Operation{2,2} = 'Str';               unit{2,2}  = 'mm';
    alldata = [x(:).*1e3 y(:).*1e3 exx(:) eyy(:) exy(:)]; % [m]
    if exist('stPs','var')
        files{2,2} = [newfile '\mm_Strain_' num2str(stPs) '.mat'];
    else
        files{2,2} = [newfile '\mm_Strain.mat'];
    end
    save(files{2,2},'alldata')
    
    Operation{2,3} = 'Str';               unit{2,3}  = 'um';
    alldata = [x(:).*1e6 y(:).*1e6 exx(:) eyy(:) exy(:)]; % [m]
    if exist('stPs','var')
        files{2,3} = [newfile '\um_Strain_' num2str(stPs) '.mat'];
    else
        files{2,3} = [newfile '\um_Strain.mat'];
    end
    save(files{2,3},'alldata')
    stif{1}='E';
    stif{2}='E';

fprintf ('DONE\n\n');