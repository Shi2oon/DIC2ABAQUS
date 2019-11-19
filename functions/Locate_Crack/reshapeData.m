function [ dataum ] = reshapeData( raw_data )
%PROCESS_DATA Summary of this function goes here
%   Detailed explanation goes here
x  = raw_data(:,1);
y  = raw_data(:,2);
ux = raw_data(:,3);
uy = raw_data(:,4);

xVec = unique(x);
yVec = unique(y);

% nDataPoints = length(x);

%Define grid
[xMap,yMap] = meshgrid(xVec,yVec);
[nRows, nCols] = size(xMap);

% nGridPoints = length(xMap(:));

uxMap = NaN(nRows, nCols); %Initialise
uyMap = NaN(nRows, nCols); %Initialise

for iRow = 1:nRows % loop rows
    for iCol = 1:nCols % loop cols
        xt = xMap(iRow,iCol);
        yt = yMap(iRow,iCol);
        idx = find(and(x==xt,y==yt)); %find linear index of point corresponding to xt,yt;
        if ~isempty(idx)
            uxt = ux(idx(1));
            uyt = uy(idx(1));
            uxMap(iRow,iCol) = uxt;
            uyMap(iRow,iCol) = uyt;
        end
    end
end

dataum.X = xMap;
dataum.Y = yMap;
% dataum.Uy = uyMap;
% dataum.Ux = uxMap;
threshold = 0.95;
[ uxMap ] = dispFieldSmoothing( uxMap, threshold );
dataum.Ux = uxMap;
[ uyMap ] = dispFieldSmoothing( uyMap, threshold );
dataum.Uy = uyMap;
end

