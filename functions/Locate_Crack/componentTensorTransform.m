function [NewTensorMap] = componentTensorTransform(OldTensorMap,a)

% Input: 
% series of 9 maps representing the spatial variation of the 9
% components of a second rank tensor, with respect to an old co-ordinate
% system.  The transformation from the old to the new co-ordinate system is
% defined by the direction cosine matrix a.

% Output: 
% series of 9 maps representing the spatial variation of the
% components of the same system, this time with respect to the new
% co-ordinate system

% R = @(theta)[cosd(theta) sind(theta) 0;-sind(theta) cosd(theta) 0;0 0 1];
% A0 = R(theta)*A0*R(theta)';
[fi,fj] = size(OldTensorMap);

NewTensorMap = cell(fi,fj);

parfor i=1:fi
    for j=1:fj
        NewTensorMap{i,j} = zeros(size(OldTensorMap{i,j}));
    end
end

parfor i=1:fi
    for j=1:fj
        % Perform einstein summation for component {i,j} 
        for k = 1:fi
            for l = 1:fj
                NewTensorMap{i,j} = NewTensorMap{i,j} + a(i,k).*a(j,l).*OldTensorMap{k,l};
            end
        end
    end
end