function [ dataum ] = reshapeAbqusData( raw_data )
%PROCESS_DATA Summary of this function goes here
%   Detailed explanation goes here
x  = raw_data(:,1);
y  = raw_data(:,2);
S11 = raw_data(:,3);
S22 = raw_data(:,4);
S12 = raw_data(:,5);
Mises = raw_data(:,6);

xVec = unique(x);
yVec = unique(y);

% nDataPoints = length(x);

%Define grid
[xMap,yMap] = meshgrid(xVec,yVec);
[nRows, nCols] = size(xMap);

% nGridPoints = length(xMap(:));

S11Map = NaN(nRows, nCols); %Initialise
S22Map = NaN(nRows, nCols); %Initialise
S12Map = NaN(nRows, nCols); %Initialise
MisesMap = NaN(nRows, nCols); %Initialise

for iRow = 1:nRows % loop rows
    for iCol = 1:nCols % loop cols
        xt = xMap(iRow,iCol);
        yt = yMap(iRow,iCol);
        idx = find(and(x==xt,y==yt)); %find linear index of point corresponding to xt,yt;
        if ~isempty(idx)
            S11Map(iRow,iCol)  = S11(idx(1));
            S22Map(iRow,iCol)  = S22(idx(1));
            S12Map(iRow,iCol)  = S12(idx(1));
            MisesMap(iRow,iCol) = Mises(idx(1));            
        end
    end
end

dataum.X   = xMap;
dataum.Y   = yMap;
threshold  = 0.95;
[ S11 ] = dispFieldSmoothing( S11Map, threshold );
if mean(mean((S11)))==0
else
    [ S11Map ]    = dispFieldSmoothing( S11Map, threshold );
    [ S22Map ]    = dispFieldSmoothing( S22Map, threshold );
    [ S12Map ]    = dispFieldSmoothing( S12Map, threshold );
    [ MisesMap ]    = dispFieldSmoothing( MisesMap, threshold );
end

    dataum.S11    = S11Map;
    dataum.S22    = S22Map;
    dataum.S12    = S12Map;
    dataum.Mises  = MisesMap;
    
end

