function [a] = DirectionCosine(theta)
% Creates matrix of direction cosines due to a transformation of 3D axes
% through an angle 'theta' anticlockwise about the 3 axis.

% 'theta' must be specified in radians

a = zeros(3,3); % initialise a
% Rz = @(theta)[cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];

% this is actually the transformation of Rx
% Rx=@(theta) [1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
% Rz = transpose(Rx(theta));
a(1,1) = cos(theta);
a(1,2) = sin(theta);
a(2,1) = -sin(theta);
a(2,2) = cos(theta);
a(3,3) = 1;