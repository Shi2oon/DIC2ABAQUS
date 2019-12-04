function [mesh,dataIn] =  loadUd(input_unit,scan_dir, scan_ID,Operation,dataIn)
mesh.Operation     = Operation;
mesh.operation     = 'Norm'; % for origiinal data and ''West' for correction data#

if mesh.Operation == 'xED'
    %% HR-EBSD data% E strain  S stress  A Deformation 
    [dataIn.Maps,alldata]   = GetGrainData(fullfile(scan_dir,[scan_ID '.mat']),dataIn.unique); 
    dataIn.results     = dataIn.Maps.SavingD;
    input_unit         = dataIn.Maps.units.xy ;
    mesh.OpS           = 'xED';
    mesh.Operation     = 'Str';
    dataIn.type        = 'A'; 
    dataIn.nu = dataIn.Maps.nu;     dataIn.E = dataIn.Maps.E;
    
else
    try
        load(fullfile(scan_dir,[scan_ID '.mat']));
    end
    if exist ('alldata')==0
        tmp = importdata(fullfile(scan_dir,[scan_ID '.dat']));
        alldata = [tmp.data(:,1) tmp.data(:,2) tmp.data(:,3) tmp.data(:,4)];        
    end
 dataIn.results  = [dataIn.fullpath '\Int2ABAQUS'];    mkdir(dataIn.results)
end
xx              = unique(alldata(:,1),'first'); 
mesh.x.length   = length(xx);
mesh.x.xi       = (max(xx))-min(xx)/2;
alldata(isnan(alldata))=0;

%% strain data
    % Convert [x y] to metres
switch input_unit
    case 'mm'
        % convert mm --> m
        alldata(:,1:2) = alldata(:,1:2).* 10^-3; 
    case 'um'
        % convert um --> m
        alldata(:,1:2) = alldata(:,1:2).* 10^-6; % convert pixels --> m
    case 'm'
        %nothing
    otherwise
        error('Invalid input_unit')
end

%%
mesh.selct = input('\nSlow with percision (S) or fast with less percision (F)?   ','s');
if mesh.selct == 'S' || mesh.selct =='s'
    if  any(ismember(fields(mesh),'OpS'))~=1
        [dataIn.Maps] = reshapeStrain (alldata);
        dataIn.Maps.SavingD = dataIn.results;
    end
    dataIn.Maps = xEBSD2DIS(dataIn.Maps);     

elseif mean(mesh.Operation=='Str') && (mesh.selct == 'F' || mesh.selct == 'f')
dataIn.stressstat  = 'plane_strain';
dataIn.input_unit  = 'm';
alldata(isnan(alldata))=0;
dispsU_1   = alldata(:,1:2);

%% Interpolate strain data on a regular grid.
straindata = [alldata(:,1:2) alldata(:,3:5)]; %exx eyy exy
% straindata = alldata(:,3:5); %exx eyy exy
str_limits = [min(dispsU_1);max(dispsU_1)];
nbpts      = [length(unique(dispsU_1(:,1))) length(unique(dispsU_1(:,2)))];
pas        = (str_limits(2,:)-str_limits(1,:))./nbpts;
xrdXq      = str_limits(1,1):pas(1):str_limits(2,1);
xrdYq      = str_limits(1,2):pas(2):str_limits(2,2);
% [Xq, Yq]   = meshgrid(xrdXq,xrdYq);

% Fxx = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), straindata(:,1));
% Fyy = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), straindata(:,2));
% Fxy = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), straindata(:,3));
% intpExx = Fxx(Xq,Yq); 
% intpEyy = Fyy(Xq,Yq); 
% intpExy = Fxy(Xq,Yq);
% straindata = [Xq(:) Yq(:) intpExx(:) intpEyy(:) -1.*intpExy(:)];

%% for HR-EBSD Data to injetc stress data 
    if  any(ismember(fields(mesh),'OpS'))==1
        if sum(strfind(mesh.OpS,'xED'))~=0
        stressdata = [alldata(:,1:2) alldata(:,6:end)];
        if      strcmp(dataIn.Maps.units.S,'KPa');         sf = 1e-3;
        elseif  strcmp(dataIn.Maps.units.S,'MPa');         sf = 1e-6;
        elseif strcmp(dataIn.Maps.units.S,'GPa');          sf = 1e-9;   
        elseif strcmp(dataIn.Maps.units.S,'Pa');           sf = 1;  
        end
%         Fxx = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), stressdata(:,1));
%         Fyy = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), stressdata(:,2));
%         Fxy = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), stressdata(:,3));
%         intpSxx = Fxx(Xq,Yq); 
%         intpSyy = Fyy(Xq,Yq); 
%         intpSxy = Fxy(Xq,Yq);
%         stressdata = [Xq(:) Yq(:) intpSxx(:) intpSyy(:) -1.*intpSxy(:)];
        dataIn.stressstat    = 'plane_stress';
        dataIn.Maps.units.xy = dataIn.input_unit;
        dataIn.nu            = dataIn.Maps.nu;          dataIn.E  = dataIn.Maps.E*sf; 
        mesh.stressdata      = stressdata.*sf;
        end
    end
    
%% Build the grid where strain points are gauss points
baseXq = [xrdXq(1)-(pas(1)/2) xrdXq+(pas(1)/2)];
baseYq = [xrdYq(1)-(pas(2)/2) xrdYq+(pas(2)/2)];
[Xbase, Ybase] = meshgrid(baseXq,baseYq);
dispsU_1 = [Xbase(:) Ybase(:)];
dispsd_1 = nan(size(dispsU_1));

%Guess the dimensions of the image
tmp(1:2,:)  = dispsU_1'; 
tmp(3:4,:)  = dispsd_1';
mesh.winDIC = GetSize(tmp(2,:));

%sort the dataset
tmp       = sortrows(tmp',[1 2])';
mesh.UDIC = tmp(1:2,:);
mesh.dDIC = tmp(3:4,:);
mesh.straindata = straindata;
mesh.winFE = [2 2; mesh.winDIC(1)-1 mesh.winDIC(2)-1];
end
dataIn.alldata = alldata;
