function MaterialsAndMapDefination(fileID,SaveD,maskdim,len,folder,MatP,offset)
fprintf(fileID,'# An error at ALL ELLEMENTS indicate you need to re-define \n');
fprintf(fileID,'# the crak vector Change and submit the job again\n');
fprintf(fileID,'# for reports check report, field output (produce values vs. elements \n');
fprintf(fileID,'# you can find x and y values in .inp files under Nodes\n');
dicInputPath = pythonFileName(SaveD);
fprintf(fileID,'dicInputPath = "%s"; \n',dicInputPath);
% ATTENTION : crackpoints should be defined with the crack tip in first position.....
% If a crackpoint is on a border, it's coordinate has to be EXACTLY the 
% coordinate of the border defined in the DICinput file
fprintf(fileID,'crackpoints = ((%f,%f),(%f,%f));\n',....
                maskdim.xo(1),maskdim.yo(1),maskdim.xo(2),maskdim.yo(2));
% Area without boundary conditions;masksVal are Xmin Ymin Xmax Ymax
fprintf(fileID,'masksVal=((%f,%f,%f,%f),);\n',....
                maskdim.xm(2),maskdim.ym(2),maskdim.xm(1),maskdim.ym(1));
%Nb of subsets masked around the crack
fprintf(fileID,'dangerZone = %d;\n', abs(maskdim.ds1-maskdim.ds2)); 
% NB of J-integral contours
fprintf(fileID,'nbCtrJint  = %d;\n', len);
% Where to write the results
outputPth = pythonFileName(fullfile(folder, 'KJ_Output.txt'));
fprintf(fileID,'outputPth  = "%s";\n', outputPth);

% material's parameters Si (mm) units, Elastic modulus Pa, Poisson's ratio
% Yield stress Pa, Exponent, Yield offset
fprintf(fileID,'MaterialName = "%s";\n',MatP.Mat);

if      MatP.type == 'E'
        fprintf(fileID,'matLaw = "Elastic";\n');
% Extract K or not using the interaction integral method / Only valid for 'Elastic Model'
        fprintf(fileID,'extractK   = 1;\n');   
        if MatP.E < 1e6;    MatP.E = MatP.E*1e9;	
            disp('Check Modulus Unit .. Corretced');    end
        MatP.E = MatP.E*offset^2;
        fprintf(fileID,'matParams=((%d,%.3f),);\n',MatP.E,MatP.nu);
elseif  MatP.type == 'R'
        fprintf(fileID,'matLaw = "Ramberg-Osgood";\n');
        if MatP.E < 1e6;    MatP.E = MatP.E*1e9;	MatP.yield = MatP.yield*1e9;
            disp('Check Modulus and Yield stress Units .. Corretced');    end
        MatP.E = MatP.E*offset^2;       MatP.yield = MatP.yield*offset^2;
        fprintf(fileID,'matParams=((%d,%.3f,%d,%.3f,%.3f),);\n',...
                MatP.E,MatP.nu,MatP.yield,MatP.Exponent, MatP.Yield_offset);
        fprintf(fileID,'extractK   = 0;\n'); 
elseif  MatP.type == 'A'
        fprintf(fileID,'matLaw = "Elastic-Anisotropic";\n');
        fprintf(fileID,'extractK   = 1;\n'); 
        
        fprintf(fileID,'# You will need to define local coordinate to do so\n');
        fprintf(fileID,'# Select the Proprety Module, then in the toolbar select Assign\n');
        fprintf(fileID,'# Material Orientation, select all the Module, use Default CSYS, OK \n');
        %stiffness tensors [Pa]
    if MatP.Stiffness(1,1) < 1e6; 	MatP.Stiffness = MatP.Stiffness.*1e9;       	
        disp('Check stiffness tensors Units .. Corretced');    end
        MatP.Stiffness = MatP.Stiffness.*offset^2;
        fprintf(fileID,'C11 = %d;\tC12 = %d;\tC13 = %d;\tC14 = %d;\tC15 = %d;\tC16 = %d;\n',...
        MatP.Stiffness(1,1),MatP.Stiffness(1,2),MatP.Stiffness(1,3),...
        MatP.Stiffness(1,4),MatP.Stiffness(1,5),MatP.Stiffness(1,6));		
                   fprintf(fileID,'C22 = %d;\tC23 = %d;\tC24 = %d;\tC25 = %d;\tC26 = %d;\n',...
        MatP.Stiffness(2,2),MatP.Stiffness(2,3),MatP.Stiffness(2,4),...
                            MatP.Stiffness(2,5),MatP.Stiffness(2,6));
                              fprintf(fileID,'C33 = %d;\tC34 = %d;\tC35 = %d;\tC36 = %d;\n',...
        MatP.Stiffness(3,3),MatP.Stiffness(3,4),MatP.Stiffness(3,5),MatP.Stiffness(3,6));
                                         fprintf(fileID,'C44 = %d;\tC45 = %d;\tC46 = %d;\n',...
        MatP.Stiffness(4,4),MatP.Stiffness(4,5),MatP.Stiffness(4,6));
                                                    fprintf(fileID,'C55 = %d;\tC56 = %d;\n',...
        MatP.Stiffness(5,5),MatP.Stiffness(5,6));
                                                               fprintf(fileID,'C66 = %d;\n',...
        MatP.Stiffness(6,6));
end
end