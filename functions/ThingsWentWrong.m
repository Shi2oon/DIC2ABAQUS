function [alldata,Dir] = ThingsWentWrong(Dir,Maps)
% this function works if things failed after intergation and you don't want
% to repeat things again. first make sure you uploaded the Data Uxy.mat
% file before running this function
    Dir.results     = Maps.SavingD;
%     Dir.input_unit  = Maps.units.xy ;
    Dir.Stiffness   = Maps.Stiffness;
    Dir.type        = 'A'; 
	Dir.nu          = Maps.nu;     
	Dir.E           = Maps.E;
    alldata         = [Maps.X1(:), Maps.Y1(:),Maps.Ux(:),Maps.Uy(:)];
    Dir.Maps        = Maps;
end