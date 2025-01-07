function [ alldata,dataum ] = reshapeStrainData( raw_data )
%PROCESS_DATA Summary of this function goes here
%   Detailed explanation goes here

x  = raw_data(:,1);
y  = raw_data(:,2);
Exx = raw_data(:,3);
Eyy = raw_data(:,4);
Exy = raw_data(:,5);
    
xVec = unique(x);
yVec = unique(y);

% nDataPoints = length(x);
% if length(xVec) <	length(raw_data(:,1))*0.5

%Define grid
[xMap,yMap] = meshgrid(xVec,yVec);
[nRows, nCols] = size(xMap);

% nGridPoints = length(xMap(:));

ExxMap = NaN(nRows, nCols); %Initialise
EyyMap = NaN(nRows, nCols); %Initialise
ExyMap = NaN(nRows, nCols); %Initialise

if size(raw_data,2) == 9
    z  = raw_data(:,3);
    zVec = unique(z);
    Exx = raw_data(:,4);
    Eyy = raw_data(:,5);
    Ezz = raw_data(:,6);
    Exy = raw_data(:,7);
    Exz = raw_data(:,8);
    Eyz = raw_data(:,9);
    [xMap,yMap,zMap] = meshgrid(xVec,yVec,zVec);
    [nRows, nCols , nDep] = size(xMap);
    ExxMap = NaN(nRows, nCols, nDep); %Initialise
    EyyMap = NaN(nRows, nCols, nDep); %Initialise
    EzzMap = NaN(nRows, nCols, nDep); %Initialise
    ExyMap = NaN(nRows, nCols, nDep); %Initialise
    ExzMap = NaN(nRows, nCols, nDep); %Initialise
    EyzMap = NaN(nRows, nCols, nDep); %Initialise
end

for iRow = 1:nRows % loop rows
    for iCol = 1:nCols % loop cols
        
        if size(raw_data,2) == 9 %% 3D
            for iDep = 1:nDep
                xt = xMap(iRow,iCol,iDep);
                yt = yMap(iRow,iCol,iDep);
                zt = zMap(iRow,iCol,iDep);
                idx = find(x==xt & y==yt & z==zt); 
                if ~isempty(idx)
                    Exxt = Exx(idx(1));
                    Eyyt = Eyy(idx(1));
                    Ezzt = Ezz(idx(1));
                    Exyt = Exy(idx(1));
                    Exzt = Exz(idx(1));
                    Eyzt = Eyz(idx(1));
                    ExxMap(iRow,iCol,iDep) = Exxt;
                    EyyMap(iRow,iCol,iDep) = Eyyt;
                    EzzMap(iRow,iCol,iDep) = Ezzt;
                    ExyMap(iRow,iCol,iDep) = Exyt;
                    ExzMap(iRow,iCol,iDep) = Exzt;
                    EyzMap(iRow,iCol,iDep) = Eyzt;
                end
            end 
            
        else %% 2D
            xt = xMap(iRow,iCol);
            yt = yMap(iRow,iCol);
            idx = find(and(x==xt,y==yt)); %find linear index of point corresponding to xt,yt;
            if ~isempty(idx)
                Exxt = Exx(idx(1));
                Eyyt = Eyy(idx(1));
                Exyt = Exy(idx(1));
                ExxMap(iRow,iCol) = Exxt;
                EyyMap(iRow,iCol) = Eyyt;
                ExyMap(iRow,iCol) = Exyt;
            end
        end
    end
end

dataum.X1 = xMap;
dataum.Y1 = yMap;
% dataum.Uy = uyMap;
% dataum.Ux = uxMap;
% threshold = 0.95;
% [ ExxMap ] = dispFieldSmoothing( ExxMap, threshold );
dataum.Exx = ExxMap;
% [ EyyMap ] = dispFieldSmoothing( EyyMap, threshold );
dataum.Eyy = EyyMap;
% [ ExyMap ] = dispFieldSmoothing( ExyMap, threshold );
dataum.Exy = ExyMap;
alldata = [dataum.X1(:)  dataum.Y1(:)  dataum.Exx(:) dataum.Eyy(:) dataum.Exy(:)];
if size(raw_data,2) == 9
    dataum.Z1 = zMap;
    dataum.Ezz = EzzMap;
    dataum.Exy = ExyMap;
    dataum.Exz = ExzMap;
    dataum.Eyz = EyzMap;
    alldata = [dataum.X1(:)  dataum.Y1(:)  dataum.Z1(:)  dataum.Exx(:) dataum.Eyy(:) ...
               dataum.Ezz(:) dataum.Exy(:) dataum.Exz(:) dataum.Eyz(:)];
end
    
%{
else
	disp('the data is highly non-uniform, I will now re-arrange the data');
%     scatter3(alldata(:,1), alldata(:,2), alldata(:,3),[],alldata(:,3)); view([0 90])
	Fx = scatteredInterpolant(raw_data(:,1), raw_data(:,2), raw_data(:,3),'natural','nearest');
	Fy = scatteredInterpolant(raw_data(:,1), raw_data(:,2), raw_data(:,4),'natural','nearest');
	X = linspace(min(raw_data(:,1)),max(raw_data(:,1)),300);
	Y = min(raw_data(:,2)):abs(X(2)-X(1)):max(raw_data(:,2));
	[dataum.X1,dataum.Y1] = meshgrid(X,Y);
	Ux = Fx(dataum.X1(:),dataum.Y1(:));
	Uy = Fy(dataum.X1(:),dataum.Y1(:));
    dataum.Ux = reshape(Ux,length(Y),length(X));
    dataum.Uy = reshape(Uy,length(Y),length(X));
	% scatter3(x(:),y(:),Ux,[],Ux); view([0 90])
end
%}


