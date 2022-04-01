function [Output] = Crop(Input,Xcrop,Ycrop)
% Takes map 'Input' and crops it using parameters Xcrop,Ycrop which are the
% row and column indices of the cropped data in the original map


Output = Input(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));