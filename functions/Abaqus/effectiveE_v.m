function [E,v,G,Co] = effectiveE_v(C) 
% this function caclulates effective Youn modulus and Possion ratio for a
% ansitropic material based on this paper
% Reference: https://doi.org/10.3390/cryst8080307
BV = (C(1,1)+2*C(1,2))/3;               % Voigt bulk modulus
GV = (C(1,1)-C(1,2)+3*C(4,4))/5;        % Voigt shear modulus

S = C^-1;
BR = 1/(3*S(1,1)+6*S(1,2));             % Reuss bulk modulus
GR = 5/(4*S(1,1)-4*S(1,2)+3*S(4,4));   % Reuss shear modulus

B = (BR+BV)/2;                          % Hill’s average bulk modulus
G = (GR+GV)/2;                          % Hill’s average shear modulus
E = 9*B*G/(3*B+G);                      % Young’s modulus (E)
v = (3*B-E)/(6*B);                      % Poisson’s ratio
Co = [];

% K = (C(1,1)+C(2,2)+C(3,3)+2*(C(1,2)+C(2,3)+C(1,2)))/9; % istropic shear Modulus
% Gv = (C(1,1)+C(2,2)+C(3,3)-(C(1,2)+C(2,3)+C(1,2))+2*(C(4,4)+C(5,5)+C(6,6)))/15; % Bulk Modulus

%% Paper: What is the Young’s Modulus of Silicon?
Cc =C^-1;
Co.Ex = 1/Cc(1,1);
Co.Ey = 1/Cc(2,2);
Co.Ez = 1/Cc(3,3);
Co.Gxy = 1/Cc(4,4);
Co.Gxz = 1/Cc(5,5);
Co.Gyz = 1/Cc(6,6);
Co.vxy = -Co.Ey*Cc(1,2);
Co.vxz = -Co.Ez*Cc(1,3);
Co.vyz = -Co.Ez*Cc(2,3);

Co.C = [ 1/Co.Ex       -Co.vxy/Co.Ey   -Co.vxz/Co.Ez 0   0   0
    -Co.vxy/Co.Ex    1/Co.Ey      -Co.vyz/Co.Ez 0   0   0
    -Co.vxz/Co.Ex   -Co.vyz/Co.Ey  1/Co.Ez     0   0   0
    0          0           0 	1/Co.Gyz 0   0
    0          0           0   0   1/Co.Gxz 0
    0          0           0   0   0   1/Co.Gxy];
Co.C = Co.C^-1;

%% a different approach as sometime the first approach sometimes 
% delivers minus results!
if G<0 || E<0 || v<0 || v > 0.5
    v = 1+(Co.vxz+Co.vyz)/2;
    if v > 0.5 % for metals
        v = abs((3*BV-E)/(6*BV));
    end
end