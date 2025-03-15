function [folder,CAE ] = PrintRunCode(MatP, maskdim,SaveD,len,offset)
[folder,~,~] = fileparts(SaveD);
folder = fullfile(folder, 'Abaqus Output');     mkdir(folder);
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

% % unix(['abaqus cae ','noGUI','=Abaqus_script.py']);   % Unix system
% Pause to ensure the script has time to execute
pause(10); % Adjust the pause duration as needed

% Close the Abaqus CAE GUI
[status, cmdout] = system('tasklist /FI "IMAGENAME eq abq*.*" /FO CSV');

if status == 0
    % Parse the output to find the process ID (PID) of Abaqus CAE
    lines = strsplit(cmdout, '\n');
    for i = 2:length(lines)
        if contains(lines{i}, 'abq')
            tokens = strsplit(lines{i}, ',');
            pid = str2double(tokens{2});
            % Terminate the Abaqus CAE process
            system(['taskkill /PID ' num2str(pid) ' /F']);
        end
    end
end

delete(fullfile(folder , [MatP.unique '.log']));
delete(fullfile(folder , [MatP.unique '.msg']));
delete(fullfile(folder , [MatP.unique '.prt']));
delete(fullfile(folder , [MatP.unique '.sim']));
delete(fullfile(folder , [MatP.unique '.sta']));
delete(fullfile(folder , [MatP.unique '.com']));
delete(fullfile(folder , [MatP.unique 'Done.tmp']));

CAE = exist([MatP.unique '.rpt'],'file');
if ~exist([MatP.unique '.rpt'],'file')
    fprintf('\nYou can find a python script to run in Abaqus in\n%s\n',folder);
end

cd(PWD);

folder = fullfile(folder,MatP.unique);