function Crot = rotate_stiffness(C, R)
%ROTATE_STIFFNESS Rotate a 4th rank stiffness tensor
%   Crot = rotate_stiffness(C, R)
%
%   Inputs:
%       C : 3x3x3x3 tensor
%       R : 3x3 rotation matrix
%
%   Output:
%       Crot : rotated 3x3x3x3 tensor

Crot = zeros(3,3,3,3);

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                val = 0;
                for m = 1:3
                    for n = 1:3
                        for p = 1:3
                            for q = 1:3
                                val = val + R(i,m) * R(j,n) * C(m,n,p,q) * R(k,p) * R(l,q);
                            end
                        end
                    end
                end
                Crot(i,j,k,l) = val;
            end
        end
    end
end
end
