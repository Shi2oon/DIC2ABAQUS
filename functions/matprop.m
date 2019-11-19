function [mat] = matprop(E,nu,yield,stressstat,input_unit)
mat.E           = E;
mat.nu          = nu;
mat.yield       = yield;
mat.stressstate = stressstat;

if ( strcmp(mat.stressstate,'plane_strain') ) 
    % Plane Strain case
      mat.D = (mat.E/((1+mat.nu)*(1-2*mat.nu)))*[ 1-mat.nu mat.nu 0; ...
          mat.nu 1-mat.nu 0; 0 0 (1-2*mat.nu)/2];
else  
    % Plane Stress case
      mat.D = (mat.E/(1-mat.nu^2))*[1 mat.nu 0; mat.nu 1 0; 0  0 (1-mat.nu)/2];
end

%Element type (deprec)
mat.eltype  = 'Q4';

if input_unit =='m'
    mat.input_unit = 1;
elseif input_unit == 'mm'
    mat.input_unit = 1e-3;
elseif input_unit == 'um'
    mat.input_unit = 1e-6;
end



