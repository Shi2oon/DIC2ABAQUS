function [BCf, UnitOffset,stepsize] = DIC2ABAQUS_wNAN(MatP, Crack, resultsDir,angle_deg)
% create an abaqus model
% you need >>
% M4Nodes from the created mesh, thsi is nromally an output
% M4Elements which is one lement set
% dicdata (Displacement data arranged as vectors, example
% [X(:) Y(:) Ux(:) Uy(:)]
% load([resultsDir '\Intergated Uxy.mat'])
%
% Crack which contain the crack x and y coordianted with Crack(1,:) is for
% the crack tip and Crack(2,:) for the crack end >>
% MatP includes
% 	E (Young's Modulus [Pa]) or
% 	Stiffness (Stiffness tensor), if 'Elastic-Anisotropic'
% 	nu (Poisson's ratio)
%   Mat (Material Name)
%   type 'E' for Elastic or 'R' for Ramberg-Osgood or 'A' for Elastic-Anisotropic
%   if 'Ramberg-Osgood' type of material input
%   Exponent
%   Yield_offset
%   yield (Yield Stress [Pa] )
%   input_unit %'m', 'mm', 'um'
%

DICdata = importdata(resultsDir);
% if the data is not regulary grides use the code below
%{
X = linspace(min(DICdata.data(:,1)),max(DICdata.data(:,1)),ceil(sqrt(length(DICdata.data(:,1)))));
Y = linspace(min(DICdata.data(:,2)),max(DICdata.data(:,2)),ceil(sqrt(length(DICdata.data(:,2)))));
[x,y] = meshgrid(X,Y);
F = scatteredInterpolant(DICdata.data(:,1), DICdata.data(:,2), DICdata.data(:,3),'natural');
Ux = F(x, y);
F = scatteredInterpolant(DICdata.data(:,1), DICdata.data(:,2), DICdata.data(:,4),'natural');
Uy = F(x, y);
DICdata.data = [x(:),y(:),Ux(:),Uy(:)];
%}

[alldata,dataDot] = reshapeData(DICdata.data*MatP.pixel_size);
MatP.stepsize = mean(unique(diff(dataDot.X1(1,:))));
stepsize = MatP.stepsize;
[resultsDir,Kl,~] = fileparts(resultsDir);
resultsDir = [resultsDir '\' Kl];     mkdir(resultsDir);
MatP.results = resultsDir;
alldata = [alldata alldata(:,3)];


%% Create the mesh object:
% function M4 = FE_OOM(alldata,ShapeFunOrder,resultsDir)
% last update 30/11/2024.
fprintf ('Started Meshing ... ');
if size(alldata,2) == 9
    % alldata = [X(:) Y(:) Z(:) e11(:) e22(:) e33(:) e12(:) e13(:) e23(:)];

    
    NDIM  = 3;              % 3D
    NNODE = 8;
    NGP   = 8;
    [M4.Nodes,M4.Elements,M4.Xall,M4.Yall,M4.Zall,M4.E11,M4.E22,M4.E33,M4.E12,...
        M4.E13,M4.E23,M4.ScaleYN] = HexMeshAbaqus(alldata); % non uniform data
    saveas(gcf,[resultsDir '\Meshed ' num2str(NGP) '.fig']);
    saveas(gcf,[resultsDir '\Meshed ' num2str(NGP) '.png']); close
else
    % alldata = [X(:) Y(:) e11(:) e22(:) e12(:)];
    
    alldata2Dsze = size(alldata, 1);
    if MatP.type == 'U'
        NDIM  = 3;              % 3D
        z_coords = [0:MatP.zElems]'*MatP.modelThickness/MatP.zElems;
        z_coords_all = kron(z_coords, ones(size(alldata, 1), 1));
        alldata = repmat(alldata, MatP.zElems+1, 1);
        alldata = [alldata(:, [1,2]), z_coords_all, alldata(:, [3,4,3])];        
        
        [M4.Xall,M4.Yall,M4.E11,M4.E22,M4.E12,M4.Nodes,M4.Elements] = ...
            MeshPlaneStress3D(alldata); % non uniform data  
    else    
        NDIM = 2;      % 2D
        alldata = sortrows(alldata,[1,2]);
        %     try
        %     catch
        %     fprintf('FE_Mesh_Generator failed ! .. Meshing ');

        [M4.Xall,M4.Yall,M4.E11,M4.E22,M4.E12,M4.Nodes,M4.Elements] = ...
            Meshing(alldata);  % support both uniform and non uniform data
    end
    [M4.Elements,M4.E11,M4.E22,M4.E12,M4.X1,M4.X2,M4.ScaleYN] = ... % remove outliers
        WhatNaN(M4.Elements,M4.E11,M4.E22,M4.E12,M4.Xall,M4.Yall);
    
    % Mesh plot
    %{
fill(M4.X1,M4.X2,'w');
set(gca,'Ydir','normal');
s1.XDir='reverse';   s1.YDir='reverse'; axis image;axis xy; colormap jet;%axis off;
xlabel('X[\mum]');  	ylabel('Y[\mum]');
title('2D Mesh');    set(gca, 'YAxisLocation', 'right')

saveas(gcf,[resultsDir '\Meshed ' num2str(NGP) '.fig']);
saveas(gcf,[resultsDir '\Meshed ' num2str(NGP) '.png']);
    %}
end

%% missing variable
fclose all;
disp('1.  Collecting missing variables and adjust units')

Um = (dataDot.Ux.^2+dataDot.Uy.^2).^0.5;
close all;
contourf(dataDot.X1,dataDot.Y1,Um,'LineStyle','none','HandleVisibility','off');
axis image; axis off; colormap jet
set(gcf,'position',[30 50 1300 950]); c=colorbar;
if strcmpi(MatP.input_unit, 'um'); UnI = '\mum';
else  UnI = MatP.input_unit; end
c.Label.String=['U_{Mag.} [' UnI ']'];
if isempty(Crack)
    try
        xo=MatP.xo;  yo=MatP.yo;
    catch
        [xo,yo] = ginput(2);
    end

    %     title('U_Y :: Select the crack tip start from crack tip');[xo,yo] = ginput(2);
    xLin       = unique(dataDot.X1);
    [~, index] = min(abs(xLin-xo(1)));      xo(1) = xLin(index);
    yLin       = unique(dataDot.Y1);
    [~, index] = min(abs(yLin-yo(1)));      yo(1) = yLin(index);
    Crack = [xo(1) yo(1); xo(2) yo(2)];
else
    xo = Crack(1);      yo = Crack(2);
        xLin       = unique(dataDot.X1);
    [~, index] = min(abs(xLin-xo(1)));      xo(1) = xLin(index);
    yLin       = unique(dataDot.Y1);
    [~, index] = min(abs(yLin-yo(1)));      yo(1) = yLin(index);
    Crack = [xo(1) yo(1); 0 0];
end
hold on; plot(xo(1),yo(1),'pk','LineStyle','-.','MarkerSize',14,...
    'MarkerEdgeColor','w','MarkerFaceColor','k'); hold off;
legend('Crack Tip','location','best')
saveas(gcf,[resultsDir '\' MatP.unique '_DIC2ABAQUS Coodrinate.fig']);
saveas(gcf,[resultsDir '\' MatP.unique '_DIC2ABAQUS Coodrinate.tif'],'tiffn');  close all
dataDot.X1 = dataDot.X1(:);         dataDot.Y1=dataDot.Y1(:);
dataDot.Ux=dataDot.Ux(:);        dataDot.Uy=dataDot.Uy(:);
dataDot.X1(isnan(dataDot.Uy))=[];   dataDot.Y1(isnan(dataDot.Uy))=[];
dataDot.Ux(isnan(dataDot.Uy))=[];   dataDot.Uy(isnan(dataDot.Uy))=[];
DICdata = [dataDot.X1(:) dataDot.Y1(:) dataDot.Ux(:) dataDot.Uy(:)];
M4Nodes    = M4.Nodes;
M4Elements = M4.Elements;

% unit set
switch MatP.input_unit
    case 'm'   % conver to m
        UnitOffset = 1;
    case 'mm'   % conver to m
        UnitOffset = 1e-3;
    case 'um'	% conver to m
        UnitOffset = 1e-6;
    case 'nm'	% conver to m
        UnitOffset = 1e-9;
end

if MatP.type == 'U'
    UnitOffset = 1e-6; % um required by UMATs
end


%% Nodes
disp('2.  Writing Nodes');
Nodes     = cell(length(M4Nodes)+9,1);
Nodes(1)  = cellstr('*Heading');
Nodes(2)  = cellstr(['** Job name: ' Kl ' Model name: ' Kl]);
Nodes(3)  = cellstr('** Generated by: Abaqus/CAE 2016');
Nodes(4)  = cellstr('*Preprint, echo=NO, model=NO, history=NO, contact=NO');
Nodes(5)  = cellstr('**');
Nodes(6)  = cellstr('** PARTS');
Nodes(7)  = cellstr('**');
Nodes(8)  = cellstr('*Part, name=sample');
Nodes(9)  = cellstr('*Node');   %Generate Nodes in Input File

[NNode, ND] = size(M4Nodes);
if ND == 3  %2D
    for i=1:1:NNode
        Nodes(9+i) = cellstr(['      ' num2str(M4Nodes(i,1)) ',   ' num2str(M4Nodes(i,2)) ...
            ',   ' num2str(M4Nodes(i,3)) ]);
    end
elseif ND==4  %3D
    for i=1:1:NNode
        Nodes(9+i) = cellstr(['      ' num2str(M4Nodes(i,1)) ',   ' num2str(M4Nodes(i,2)) ...
            ',   ' num2str(M4Nodes(i,3)) ',   ' num2str(M4Nodes(i,4)) ]);
    end
    M4NodesBack = M4Nodes(M4Nodes(:,4)==0, :);
end

%% Generate Elements in Input File
NElem = size(M4Elements, 1);
if ND==3
    NElemInPlane = NElem; 
else
    NElemInPlane = NElem/MatP.zElems;
end


% Interpolate EBSD orientations onto element centroids
if MatP.type == 'U' || MatP.type == 'A'
    % check if MTEX is installed
    if exist('mtexdata','file') ~= 2
    %     run([pwd, '\mtex-6.0.0\startup_mtex.m'])
        error('MTEX installation not found. Please install MTEX 6.0.0 to read EBSD data.')
        return
    end
    
    % check if EBSD data is present
    filePath = fullfile(pwd, MatP.EBSDfilename);
    if exist(filePath, 'file')
        load(filePath);
    else
        error('File does not exist: %s', filePath);
    end

    if ND == 3
        M4NodeCoordX = (M4.X1)';
        M4NodeCoordY = (M4.X2)';
    else
        M4ElementsBackX1 = M4.X1(:, 1:NElemInPlane);
        M4ElementsBackX2 = M4.X2(:, 1:NElemInPlane);
        M4NodeCoordX = (M4ElementsBackX1)';
        M4NodeCoordY = (M4ElementsBackX2)';
    end

    M4.centroidsX = mean(M4NodeCoordX, 2);
    M4.centroidsY = mean(M4NodeCoordY, 2);

    EBSDinterp = interp(ebsd, M4.centroidsX, M4.centroidsY);
end


disp('3.  Writing elements');
Ele = cell(NElem + 1,1);
if strcmpi(MatP.stressstat, 'plane_strain')
    Ele(1)  =  cellstr('*ELEMENT, ELSET=Set_All, TYPE=CPE4');
elseif MatP.type == 'U'
    Ele(1)  =  cellstr('*ELEMENT, ELSET=Set_All, TYPE=C3D8');
else
    Ele(1)  =  cellstr('*ELEMENT, ELSET=Set_All, TYPE=CPS4');
end

for j=1:length(M4Elements)      % Loop for the elements in the elements set
    clearvars NNN
    NNN = num2str(M4Elements(j,1));
    for k=2:length(M4Elements(j,:))
        NNN = [NNN ', ' num2str(M4Elements(j,k)) ];
    end
    Ele(j+1)  =  cellstr(NNN);
end

%%
disp('4.  Writing element sets');
% StP = sort(M4Elements(:,2));
% SelectSet    = cell(ceil(length(StP)/16)+17,1);
% SelectSet(1) =  cellstr('*Elset, elset=Set_All, internal');
% for k=1:16:length(StP)-1
%     clearvars NNN
%     NNN = [num2str(StP(k+1)) ', '];
%     if k < length(StP)-16
%         for ik = 2:15
%             NNN = [NNN num2str(StP(k+ik)) ', '];
%         end
%         SelectSet(ceil(k/16)+1) = cellstr(strcat([NNN, '' num2str(StP(k+16))]));
%     else
%         for ik = 2:length(StP)-k-1
%             NNN = [NNN num2str(StP(k+ik)) ', '];
%         end
%         SelectSet(ceil(k/16)+1) = cellstr(strcat([NNN, '' num2str(StP(length(StP)))]));
%     end
% end

SelectSet  = cell(NElemInPlane*2,1);
if ND == 3
    for i = 1:NElemInPlane
        SelectSet(2*i-1) = cellstr(['*Elset, elset=Set_', num2str(M4Elements(i, 1))]);
        SelectSet(2*i) = cellstr([' ', num2str(M4Elements(i, 1))]);
    end
else
    for i = 1:NElemInPlane
        SelectSet(2*i-1) = cellstr(['*Elset, elset=Set_', num2str(M4Elements(i, 1))]);
        vals = strings(MatP.zElems,1);
        for j = 1:MatP.zElems
            vals(j) = num2str(M4Elements(i+NElemInPlane*(j-1), 1));
        end
        SelectSet(2*i) = cellstr(strjoin(vals, ", "));
    end
end



%%

if MatP.type == 'U' || MatP.type == 'A'
    disp('5a.  Writing Section Assignments');
    
    OrientationDef = cell(3, 1);
    OrientationDef(1)  =  cellstr('*Orientation, name=Ori-All');
    OrientationDef(2)  =  cellstr('1., 0., 0., 0., 1., 0.');
    OrientationDef(3)  =  cellstr('3, 0.');
    
    SectAss = cell(NElemInPlane, 1);
    for i = 1:NElemInPlane
        SectAss(i) = cellstr(['*Solid Section, elset=Set_', num2str(M4Elements(i, 1)), ', orientation=Ori-All, material=Grain_', num2str(M4Elements(i, 1)), '_MAT']);
    end
       
    
    disp('5b.  Writing Assembly');
    Assmp     =  cell(13,1);

    Assmp(1)  =  cellstr(',');
    Assmp(2)  =  cellstr('*End Part');
    Assmp(3)  =  cellstr('**');
    Assmp(4)  =  cellstr('**');
    Assmp(5) =  cellstr('** ASSEMBLY');
    Assmp(6) =  cellstr('**');
    Assmp(7) =  cellstr('*Assembly, name=Assembly');
    Assmp(8) =  cellstr('**');
    Assmp(9) =  cellstr('*Instance, name=sample-1, part=sample');
    Assmp(10) =  cellstr('*End Instance');
    Assmp(11) =  cellstr('**');
    Assmp(12) =  cellstr('**');
    Assmp(13) =  cellstr('**');

    Assmp = [OrientationDef; SectAss; Assmp];
else
    disp('5.  Writing Materials Orientation and Assembly');
    Assmp     =  cell(16,1);
    Assmp(1)  =  cellstr('*Orientation, name=Ori-All');
    Assmp(2)  =  cellstr('1., 0., 0., 0., 1., 0.');
    Assmp(3)  =  cellstr('3, 0.');
    Assmp(4)  =  cellstr('** Section: Section-1');
    Assmp(5)  =  cellstr(['*Solid Section, elset=Set_All, orientation=Ori-All, material=' MatP.Mat]);
    Assmp(6)  =  cellstr(',');
    Assmp(7)  =  cellstr('*End Part');
    Assmp(8)  =  cellstr('**');
    Assmp(9)  =  cellstr('**');
    Assmp(10) =  cellstr('** ASSEMBLY');
    Assmp(11) =  cellstr('**');
    Assmp(12) =  cellstr('*Assembly, name=Assembly');
    Assmp(13) =  cellstr('**');
    Assmp(14) =  cellstr('*Instance, name=sample-1, part=sample');
    Assmp(15) =  cellstr('*End Instance');
    Assmp(16) =  cellstr('**');
end


%%
disp('6.  Writing instances and node lists');
% M4Nodes = M4.Nodes;       M4Elements = M4.Elements;

if ND == 4
    z_coords_DIC = kron(z_coords, ones(size(DICdata, 1), 1));
    DICdata = repmat(DICdata, MatP.zElems+1, 1);
    DICdata = [DICdata(:, 1:2), z_coords_DIC, DICdata(:, 3:4)];
end
    
[IX,~] = ismember(M4Nodes(:,2:3),DICdata(:,1:2));
ix     = sum(IX,2);
iC = ix;    iC(iC~=2)=[];

if length(iC)==length(DICdata)
    M4Nodes(ix~=2,:) = [];
    Datum  = [M4Nodes(:,1) DICdata];
else
    for io=2-1:length(DICdata)
        if ~isnan(DICdata(io,3)) && ND ==3
            id = M4Nodes(:,2:3) == DICdata(io,1:2);
            id = sum(id,2);
            Datum(io,1) = find(id==2);
        elseif ~isnan(DICdata(io,3))            
            id = M4Nodes(:,2:4) == DICdata(io,1:3);
            id = sum(id,2);
            Datum(io,1) = find(id==3);
        end
    end
    Datum = [Datum(:,1) DICdata];
end
%{
Datum  = M4Nodes(M4.Elements(:,1),1);
Fx = scatteredInterpolant(dicdata(:,1),dicdata(:,2),dicdata(:,3),'natural');
Fy = scatteredInterpolant(dicdata(:,1),dicdata(:,2),dicdata(:,4),'natural');
Datum(:,4) = Fx(Datum(:,2),Datum(:,3));
Datum(:,5) = Fy(Datum(:,2),Datum(:,3));
Datum(Datum(:,5)==0,:)=NaN;
Datum(isnan(Datum(:,5)),:)=[];
%}

if ND == 4
    Datum3D = Datum;
    Datum = Datum(1:size(Datum,1)/(MatP.zElems+1), :);
    Datum = Datum(:, [1:3,5,6]);
end

patchASSEM = cell(2*length(Datum)+4,1);

for k=1:size(Datum,1)
    if ND == 3
        if ~isnan(Datum(k,4))
            patchASSEM(k*2-1) =  cellstr(strcat(['*Nset, nset=_PickedSet',...
                num2str(Datum(k,1)), ', internal, instance=sample-1']));
            patchASSEM(k*2)   =  cellstr([num2str(Datum(k,1)),',']);
        end
    else
        if ~isnan(Datum(k,4))
            patchASSEM(k*2-1) =  cellstr(strcat(['*Nset, nset=_PickedSet',...
                num2str(Datum(k,1)), ', internal, instance=sample-1']));

            vals = strings(MatP.zElems + 1, 1);
            for j = 1:MatP.zElems + 1
                vals(j) = num2str(Datum(k,1) + alldata2Dsze*(j-1));
            end
            patchASSEM(k*2) = cellstr(strjoin(vals, ", "));
        end
    end
end

if ND == 3
    kk = Datum(k,1);
    % crack tip
    id = M4Nodes(:,2:3) == Crack(1,:);	id = sum(id,2);  Crack1(1) = find(id==2);
    patchASSEM(k*2+1) =  cellstr(strcat('*Nset, nset=_PickedSet',num2str(kk+1), ...
        ', internal, instance=sample-1'));
    patchASSEM(k*2+2) =  cellstr([num2str(Crack1(1)),',']);
    patchASSEM(k*2+3) =  cellstr(strcat('*Nset, nset=_PickedSet',num2str(kk+2), ...
        ', internal, instance=sample-1'));
    patchASSEM(k*2+4) =  cellstr([num2str(Crack1(1)),',']);
elseif ND == 4
    kk = Datum(k,1);
    % crack tip
    id = M4Nodes(:,2:3) == Crack(1,:);	id = sum(id,2);  Crack1 = find(id==2);
    for i = 1:length(Crack1)
        patchASSEM(k*2+1+(i-1)*4) = cellstr(strcat('*Nset, nset=crackSet',num2str(kk+1), ...
        '-', num2str(i), ' , instance=sample-1'));
        patchASSEM(k*2+2+(i-1)*4) = cellstr([num2str(Crack1(i)),',']);
        patchASSEM(k*2+3+(i-1)*4) = cellstr(strcat('*Nset, nset=crackSet',num2str(kk+2), ...
        '-', num2str(i), ' , instance=sample-1'));
        patchASSEM(k*2+4+(i-1)*4) =  cellstr([num2str(Crack1(i)),',']);        
    end
    
end


if ND == 4
    % node set for u3 constraint
    patchASSEMu3 = cell(2,1);
    patchASSEMu3(1) =  cellstr(strcat(['*Nset, nset=_PickedSet',...
        num2str(Datum(k,1) + 1), ', internal, instance=sample-1']));
    patchASSEMu3(2)   =  cellstr([num2str(Datum(k,1)),',']);
    patchASSEM = [patchASSEM; patchASSEMu3];
end

%% Materials Prop

disp('7.  Writing Materials propreties');
Mate     =  cell(17,1);
Mate(1)  =  cellstr('*End Assembly');
Mate(2)  =  cellstr('** ');
Mate(3)  =  cellstr('** MATERIALS');
Mate(4)  =  cellstr('** ');
Mate(5)  =  cellstr(['*Material, name=' MatP.Mat]);
    iNum = -1;
if  MatP.type == 'A'
    clear('Mate')
    MateIntro = cell(4,1);
    MateIntro(1)  =  cellstr('*End Assembly');
    MateIntro(2)  =  cellstr('** ');
    MateIntro(3)  =  cellstr('** MATERIALS');
    MateIntro(4)  =  cellstr('** ');
    
    Mate = cell(NElemInPlane*5, 1);    
    for i = 1:NElemInPlane
        Mate(5*i - 4)  =  cellstr(['*Material, name=Grain_', num2str(M4Elements(i, 1)), '_MAT']);
        Mate(5*i - 3)  =  cellstr('*Elastic, type=ANISOTROPIC');
        
        C = MatP.Stiffness;
        if C(1,1) < 1e6; 	C = C.*1e9;       	disp('Check Modulus Units .. Corrected');    end
        C = C.*UnitOffset^2;
        ori = EBSDinterp.orientations(i);
        phi1 = rad2deg(ori.phi1);
        Phi = rad2deg(ori.Phi);
        phi2 = rad2deg(ori.phi2);

        R = buildRotationMatrix(deg2rad(phi1), deg2rad(Phi), deg2rad(phi2));
        C_Voigt = C([1 2 3 6 5 4], [1 2 3 6 5 4]); % convert abaqus ordering to Voigt
        C4 = voigt2Tensor(C_Voigt); % write C as 4th rank tensor
        C4_rotated = rotate_stiffness(C4, R);
        C_Voigt_rotated = tensor2Voigt(C4_rotated);
        C_rotated = C_Voigt_rotated([1 2 3 6 5 4], [1 2 3 6 5 4]);
        C = C_rotated;      
        
        Mate(5*i - 2)  =  cellstr([' ',num2str(C(1,1)),', ',num2str(C(1,2)),', ',num2str(C(2,2)),', ',num2str(C(1,3)),...
            ', ',num2str(C(2,3)),', ',num2str(C(3,3)),', ',num2str(C(1,4)),', ',num2str(C(2,4))]);
        Mate(5*i - 1)  =  cellstr([' ',num2str(C(3,4)),', ',num2str(C(4,4)),', ',num2str(C(1,5)),', ',num2str(C(2,5)),...
            ', ',num2str(C(3,5)),', ',num2str(C(4,5)),', ',num2str(C(5,5)),', ',num2str(C(1,6))]);
        Mate(5*i - 0)  =  cellstr([' ',num2str(C(2,6)),', ',num2str(C(3,6)),', ',num2str(C(4,6)),', ',num2str(C(5,6)),...
            ', ',num2str(C(6,6))]);
    end
    Mate = [MateIntro; Mate];  

elseif MatP.type == 'R'
    Mate(6)  =  cellstr('*Deformation Plasticity');
    if MatP.E < 1e6;    MatP.E=MatP.E*1e9;	disp('Check Modulus Units .. Corrected');    end
    MatP.E = MatP.E*UnitOffset^2;
    MatP.yield=MatP.yield*UnitOffset^2;
    Mate(7)  =  cellstr([' ',num2str(MatP.E),', ',num2str(MatP.nu),...
        ', ',num2str(MatP.yield),', ' ...
        num2str(MatP.Exponent),', ',num2str(MatP.Yield_offset)]);

elseif MatP.type == 'E'
    Mate(6)  =  cellstr('*Elastic');
    if MatP.E < 1e6;    MatP.E=MatP.E*1e9;	disp('Check Modulus Units .. Corrected');    end
    MatP.E = MatP.E*UnitOffset^2;
    Mate(7)  =  cellstr([' ',num2str(MatP.E),', ',num2str(MatP.nu)]);
    Mate(8)  =  cellstr('** ');
    Mate(9)  =  cellstr('** ');

elseif MatP.type == 'U'
    clear('Mate')
    MateIntro = cell(4,1);
    MateIntro(1)  =  cellstr('*End Assembly');
    MateIntro(2)  =  cellstr('** ');
    MateIntro(3)  =  cellstr('** MATERIALS');
    MateIntro(4)  =  cellstr('** ');

    Mate = cell(NElemInPlane*5, 1);
    for i = 1:NElemInPlane
        Mate(5*i - 4)  =  cellstr(['*Material, name=Grain_', num2str(M4Elements(i, 1)), '_MAT']);
        Mate(5*i - 3)  =  cellstr(['*Depvar']);
        Mate(5*i - 2)  =  cellstr(['     ', num2str(MatP.depvar), ',']);
        Mate(5*i - 1)  =  cellstr(['*User Material, constants=6']);

        ori = EBSDinterp.orientations(i);
        phi1 = rad2deg(ori.phi1);
        Phi = rad2deg(ori.Phi);
        phi2 = rad2deg(ori.phi2);

        Mate(5*i - 0)  =  cellstr([num2str(phi1), ', ',num2str(Phi), ', ', num2str(phi2), ', ', num2str(M4Elements(i, 1)), ', ', num2str(MatP.materialID), ', ', num2str(MatP.PROPS)]);
    end

    Mate = [MateIntro; Mate];
end


Steps = cell(8, 1);
Steps(1) =  cellstr('** ----------------------------------------------------------------');
Steps(2) =  cellstr('** ');
Steps(3) =  cellstr('** STEP: Step-1');
Steps(4) =  cellstr('** ');
Steps(5) =  cellstr('*Step, name=Step-1, nlgeom=YES, inc=10000000');
Steps(6) =  cellstr('*Static');
Steps(7) =  cellstr('0.01, 1., 1e-09, 1.');
Steps(8) =  cellstr('** ');

Mate = [Mate; Steps];

%%
disp('8.  Writing Boundary Conditions');
patchBC    =  cell(4*length(Datum)+2,1);
patchBC(1) =  cellstr('** BOUNDARY CONDITIONS');
patchBC(2) =  cellstr('**');
for k=1:size(Datum,1)
    if ~isnan(Datum(k,4))
        patchBC(2+k*4-3) = cellstr(strcat(['** Name: BC-', num2str(Datum(k,1)),...
            ' Type: Displacement/Rotation']));
        patchBC(2+k*4-2) = cellstr('*Boundary');
        patchBC(2+k*4-1) = cellstr(strcat(['_PickedSet',num2str(Datum(k,1)),...
            ', 1, 1, ',num2str(Datum(k,4))]));
        patchBC(2+k*4)   = cellstr(strcat(['_PickedSet',num2str(Datum(k,1)),...
            ', 2, 2, ',num2str(Datum(k,5))]));
    end
end

if ND == 4
    % u3 constraint for 3D model
    BC_Z = cell(3,1);
    BC_Z(1) = cellstr(strcat('** Name: BC-Z Type: Displacement/Rotation'));
    BC_Z(2) = cellstr('*Boundary');
    BC_Z(3) = cellstr(strcat(['_PickedSet',num2str(Datum(k,1) + 1),', 3, 3, 0']));
    patchBC = [patchBC; BC_Z];
end

%% calculate angle
% Define the coordinates of the two points
if ~exist('angle_deg','var')
    if isempty(angle_deg)
    Crack = [xo(1) yo(1); xo(2) yo(2)];
    point1 = [Crack(1,1), Crack(1,2)]; % Coordinates of the first point
    point2 = [Crack(2,1), Crack(2,2)]; % Coordinates of the second point

    % Calculate the differences in coordinates
    vector = point2 - point1;
    angle_rad = atan2(vector(2), vector(1));

    % Convert the angle to degrees
    angle_deg = rad2deg(angle_rad);
%     if angle_deg < 1
%         angle_deg = 90;
%     end
    end
end
%%
disp('9.  Writing Output request');

OutputIntro = cell(5, 1);
OutputIntro(1)  =  cellstr('** ');
OutputIntro(2)  =  cellstr('** OUTPUT REQUESTS');
OutputIntro(3)  =  cellstr('** ');
OutputIntro(4)  =  cellstr('*Restart, write, frequency=0');
OutputIntro(5)  =  cellstr('** ');

OutputField = cell(10 ,1);
if MatP.type == 'U'
    OutputField(1)  =  cellstr('** FIELD OUTPUT: F-Output-1');
    OutputField(2)  =  cellstr('** ');
    OutputField(3)  =  cellstr('*Output, field');
    OutputField(4)  =  cellstr('*Node Output');
    OutputField(5)  =  cellstr('CF, RF, U');
    OutputField(6)  =  cellstr('*Element Output, directions=YES');
    OutputField(7)  =  cellstr('LE, PE, PEEQ, PEMAG, S, SDV');
    OutputField(8)  =  cellstr('*Contact Output');
    OutputField(9)  =  cellstr('CDISP, CSTRESS');
    OutputField(10)  =  cellstr('**');
else
    OutputField(1)  =  cellstr('** FIELD OUTPUT: F-Output-1');
    OutputField(2)  =  cellstr('** ');
    OutputField(3)  =  cellstr('*Output, field, variable=PRESELECT');
    OutputField(4)  =  cellstr('** ');
end

OutputHisJ = cell(8 + length(Crack) ,1);
OutputHisJ(1) =  cellstr('** HISTORY OUTPUT: H-Output-1');
OutputHisJ(2) =  cellstr('** ');
OutputHisJ(3) =  cellstr('*Output, history, variable=PRESELECT');
OutputHisJ(4) =  cellstr('** ');
OutputHisJ(5) =  cellstr('** HISTORY OUTPUT: OutpJint');
OutputHisJ(6) =  cellstr('** ');
Coun = round(min([length(unique(DICdata(:,1))) length(unique(DICdata(:,2)))])*0.5-1,0);
OutputHisJ(7) =  cellstr(['*Contour Integral, crack name=OutpJint_Crack-1, ' ...
    'contours=' num2str(Coun) ', crack tip nodes']);
if (ND == 3)
    OutputHisJ(8) = cellstr(['_PickedSet' num2str(kk+1) ', _PickedSet' num2str(kk+2)...
        ', ' num2str(cosd(angle_deg)) ', ' num2str(sind(angle_deg))  ', 0.']);
    OutputHisJ(9) =  cellstr('*End Step');
else
    for i = 1:length(Crack1)
        OutputHisJ(8 + (i-1)) = cellstr(['crackSet',num2str(kk+1),'-', num2str(i),...
            ', crackSet' num2str(kk+2),'-', num2str(i),', ' num2str(cosd(angle_deg)) ...
            ', ' num2str(sind(angle_deg))  ', 0.']);
    end
    OutputHisJ(8 + length(Crack1)) =  cellstr('*End Step');
end
OutputHis = OutputHisJ;

if (MatP.type ~= 'R') && (MatP.type ~= 'U')
    OutputHisJ(9) =  cellstr('** ');
    OutputHisK = cell(5, 1);
    OutputHisK(1) =  cellstr('** HISTORY OUTPUT: OutpKval');
    OutputHisK(2) =  cellstr('** ');
    OutputHisK(3) =  cellstr(['*Contour Integral, crack name=OutpKval_Crack-1, '...
        'contours=' num2str(Coun) ', crack tip nodes, type=K FACTORS']);
    OutputHisK(4) = cellstr(['_PickedSet' num2str(kk+1) ', _PickedSet' num2str(kk+2)...
        ', ' num2str(cosd(angle_deg)) ', ' num2str(sind(angle_deg))  ', 0.']);
    OutputHisK(5) =  cellstr('*End Step');
    OutputHis = [OutputHisJ; OutputHisK];
end
% end

Outputd = [OutputIntro; OutputField; OutputHis];
%% Find where to inject the BC patch
finalform = [Nodes;Ele;SelectSet;Assmp;patchASSEM;Mate;patchBC;Outputd];
% finalform = [Nodes;Ele;Assmp;patchASSEM;Mate;patchBC;Outputd];

%% Write outputfile
fprintf('10. Writing .inp file .. ');
BCf = [resultsDir '\' MatP.unique '.inp'];
fileID = fopen(BCf,'w');
for i=1:size(finalform,1)
    stri = finalform(i);
    if ~cellfun('isempty',stri)
        fprintf(fileID,'%s\n',char(stri));
    end
end
fclose(fileID);
fprintf('Done\nCheck %s for the abaqus .inp model\n',BCf);
BCf = erase(BCf,'.inp');

%% Copy UMAT files to results directory

if MatP.type == 'U'
    UMATfolderpath = fileparts(MatP.UMATfilepath);
    copyfile(UMATfolderpath, resultsDir)
end


end


