function [folder ] = Print4OUROMA(MatP, maskdim,SaveD,len)
[folder,~,~] = fileparts(SaveD);
folder = [folder '\Abaqus Output'];     mkdir(folder);
fileid = fullfile(folder, 'Abaqus_script.py');
fileID = fopen(fileid,'w');

ImportModules(fileID,folder)            %% prepare used modules
MaterialsAndMapDefination(fileID,SaveD,maskdim,len,folder,MatP) %% Materials and Map definiation 
DefineFunctions(fileID)                 %% define python functions
ReadingInputFile(fileID)                %% Reading input file
CreatPartsPrelimMeshing(fileID)       	%% ABAQUS PART CREATION & PRELIM MESHING
SaveModel = DeleteCrackedRegion(fileID,folder,MatP); %% DELETE CRACKED REGION
RemeshCrackedRegion(fileID)             %% REMESH CRACKED REGION
OutputRemeshedRegion(fileID)            %% OUTPUT REMESHED REGION
MergeModelsMeshes(fileID,SaveModel)     %% MERGE MODELS MESHES
JobMaterialHistory(fileID,MatP,len)     %% JOB, MATERIAL, HISTORYOUT
prepareBC(fileID,SaveModel)             %% PREPARE FOR BOUNDARY CONDITIONS
ApplyBCs(fileID,SaveModel)              %% APPLY BOUNDARY CONDITIONS
JobSyubmitResults(fileID,folder,MatP)   %% JOB SUBMIT AND RESULT PARSING
SaveReport(fileID,folder,MatP)          %% save report 

fclose(fileID);

%% backup to run
if  MatP.type == 'A'
fileid = fullfile(folder, 'Run_After_Fix.py');
fileID = fopen(fileid,'w');

ApplyBCs(fileID,SaveModel)              %% APPLY BOUNDARY CONDITIONS
JobSyubmitResults(fileID,folder,MatP)   %% JOB SUBMIT AND RESULT PARSING
SaveReport(fileID,folder,MatP)          %% save report 

fclose(fileID);
end

fprintf('\n Printing is Complete .. Open using Notepad ++\n');
fprintf('Done, For  results and python code Check %s \n',folder); close all