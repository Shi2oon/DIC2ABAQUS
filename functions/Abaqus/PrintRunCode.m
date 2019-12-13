function [folder ] = PrintRunCode(MatP, maskdim,SaveD,len)
[folder,~,~] = fileparts(SaveD);
folder = [folder '\Abaqus Output'];     mkdir(folder);
fileid = fullfile(folder, 'Abaqus_script.py');
fileID = fopen(fileid,'w');

ImportModules(fileID,folder)            %% prepare used modules
MaterialsAndMapDefination(fileID,SaveD,maskdim,len,folder,MatP) %% Materials and Map definiation 
DefineFunctions(fileID)                 %% define python functions
ReadingInputFile(fileID)                %% Reading input file
CreatPartsPrelimMeshing(fileID)       	%% ABAQUS PART CREATION & PRELIM MESHING
SaveModel = DeleteCrackedRegion(fileID,folder,MatP.unique); %% DELETE CRACKED REGION
RemeshCrackedRegion(fileID)             %% REMESH CRACKED REGION
OutputRemeshedRegion(fileID)            %% OUTPUT REMESHED REGION
MergeModelsMeshes(fileID,SaveModel)     %% MERGE MODELS MESHES
JobMaterialHistory(fileID,MatP.unique,MatP.type,len)     %% JOB, MATERIAL, HISTORYOUT
prepareBC(fileID,SaveModel)             %% PREPARE FOR BOUNDARY CONDITIONS
ApplyBCs(fileID,SaveModel)              %% APPLY BOUNDARY CONDITIONS
JobSyubmitResults(fileID,folder,MatP.unique)   %% JOB SUBMIT AND RESULT PARSING
SaveReport(fileID,folder,MatP.unique)          %% save report 
printTiffs(fileID,folder,'N');                 %% save tiff images for U and S
printTiffs(fileID,folder,'D');                 %% save tiff images for deformaed U and S

fclose(fileID);

%% backup to run
if  MatP.type == 'A'
    fileid = fullfile(folder, 'Run_After_Fix.py');
    fileID = fopen(fileid,'w');
    re_submit(fileID,folder,MatP.unique)           % resubmit job
    SaveReport(fileID,folder,MatP.unique)          % save report 
    printTiffs(fileID,folder,'N');                 % save tiff images for U and S
    printTiffs(fileID,folder,'D');                 % save tiff images for deformaed U and S
    fclose(fileID);
end

% fprintf('\n Printing is Complete .. Open using Notepad ++\n');
% fprintf('Done, For  results and python code Check\n%s \n',folder); close all

PWD =pwd;           cd(folder)
system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
% unix(['abaqus cae ','noGUI','=Abaqus_script.py']);   % Unix system
cd(PWD);