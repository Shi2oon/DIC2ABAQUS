function [folder,CAE ] = PrintRunCode(MatP, maskdim,SaveD,len,offset)
[folder,~,~] = fileparts(SaveD);
folder = [folder '\Abaqus Output'];     mkdir(folder);
fileid = fullfile(folder, 'Abaqus_script.py');
fileID = fopen(fileid,'w');

ImportModules(fileID,folder)            %% prepare used modules
% Materials and Map definiation 
MaterialsAndMapDefination(fileID,SaveD,maskdim,len,folder,MatP,offset) 
DefineFunctions(fileID)                 %% define python functions
ReadingInputFile(fileID)                %% Reading input file
CreatPartsPrelimMeshing(fileID,MatP.stressstat)       	%% ABAQUS PART CREATION & PRELIM MESHING
SaveModel = DeleteCrackedRegion(fileID,folder,MatP.unique); %% DELETE CRACKED REGION
RemeshCrackedRegion(fileID)             %% REMESH CRACKED REGION
OutputRemeshedRegion(fileID)            %% OUTPUT REMESHED REGION
MergeModelsMeshes(fileID,SaveModel)     %% MERGE MODELS MESHES
JobMaterialHistory(fileID,MatP.unique,MatP.type,len)     %% JOB, MATERIAL, HISTORYOUT
prepareBC(fileID,SaveModel)             %% PREPARE FOR BOUNDARY CONDITIONS
ApplyBCs(fileID,SaveModel)              %% APPLY BOUNDARY CONDITIONS
JobSyubmitResults(fileID,folder,MatP.unique)   %% JOB SUBMIT AND RESULT PARSING
SaveReport(fileID,folder,MatP.unique)          %% save report 
printTiffs(fileID,folder,'N',MatP.unique,2);   %% save tiff images for U and S
printTiffs(fileID,folder,'D',MatP.unique,2);   %% save tiff images for deformaed U and S
fclose(fileID);

%% backup to run
% if  MatP.type == 'A'
%     fileid = fullfile(folder, 'Run_After_Fix.py');
%     fileID = fopen(fileid,'w');
%     re_submit(fileID,folder,MatP.unique)           % resubmit job
%     SaveReport(fileID,folder,MatP.unique)          % save report 
%     printTiffs(fileID,folder,'N');                 % save tiff images for U and S
%     printTiffs(fileID,folder,'D');                 % save tiff images for deformaed U and S
%     fclose(fileID);
% end

% fprintf('\n Printing is Complete .. Open using Notepad ++\n');
% fprintf('Done, For  results and python code Check\n%s \n',folder); close all

%%
PWD =pwd;           cd(folder)
system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
% unix(['abaqus cae ','noGUI','=Abaqus_script.py']);   % Unix system

delete([folder '\' MatP.unique '.log']);
delete([folder '\' MatP.unique '.msg']);
delete([folder '\' MatP.unique '.prt']);
delete([folder '\' MatP.unique '.sim']);
delete([folder '\' MatP.unique '.sta']);
delete([folder '\' MatP.unique '.com']);
delete([folder '\' MatP.unique 'Done.tmp']);

CAE = exist([MatP.unique '.rpt'],'file');
if ~exist([MatP.unique '.rpt'],'file')
    fprintf('\nYou can find a python script to run in Abaqus in\n%s\n',folder);
end

cd(PWD);

folder = fullfile(folder,MatP.unique);