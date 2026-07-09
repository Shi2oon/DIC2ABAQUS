function CAE = PrintRunCode(MatP,folder)
fileid = fullfile(folder, 'Abaqus_script.py');
fileID = fopen(fileid,'w');
submitDUC2CAE_wNAN(fileID,folder,MatP.unique)
% printTiffs(fileID,folder,'N',MatP.unique,2);   %% save tiff images for U and S
% printTiffs(fileID,folder,'D',MatP.unique,2);   %% save tiff images for deformaed U and S
fclose(fileID);

%%
PWD =pwd;           cd(folder)
system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
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
delete([folder '\Job-' MatP.unique '.log']);
delete([folder '\Job-' MatP.unique '.msg']);
delete([folder '\Job-' MatP.unique '.prt']);
delete([folder '\Job-' MatP.unique '.sim']);
delete([folder '\Job-' MatP.unique '.sta']);
delete([folder '\Job-' MatP.unique '.com']);
delete([folder '\Job-' MatP.unique 'Done.tmp']);
if ~exist([MatP.unique '.rpt'],'file')
    fprintf('\nYou can find a python script to run in Abaqus in\n%s\n',folder);
end
cd(PWD);
CAE = fullfile(folder,['Job-' MatP.unique]);
end