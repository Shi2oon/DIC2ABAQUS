function [ smoothedData ] = dispFieldSmoothing( rawData, threshold )
%DISPFIELDSMOOTHING Summary of this function goes here
%   Detailed explanation goes here

tmp1 = rawData;

% fill in NaNs
tmp2 = inpaint_nans(tmp1,1);

% outlier deletion
medfiltSpacing = 10;

[tmp3] = OutDel(tmp2, medfiltSpacing, threshold);
smoothedData = inpaint_nans(tmp3,1);

end

