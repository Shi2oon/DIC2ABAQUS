function CAE = PrintRunCode(MatP,folder)
fileid = fullfile(folder, 'Abaqus_script.py');
fileID = fopen(fileid,'w');
submitDUC2CAE_wNAN(fileID,folder,MatP.unique,MatP.type)
% printTiffs(fileID,folder,'N',MatP.unique,2);   %% save tiff images for U and S
% printTiffs(fileID,folder,'D',MatP.unique,2);   %% save tiff images for deformaed U and S
fclose(fileID);

%%
PWD =pwd;           cd(folder)

system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
% % unix(['abaqus cae ','noGUI','=Abaqus_script.py']);   % Unix system
% Pause to ensure the script has time to execute
pause(10); % Adjust the pause duration as needed


if MatP.type == 'U'
    % PowerShell script to find ifort target and argument from abaqus shortcut
    psCmd = sprintf(['$wsh = New-Object -ComObject WScript.Shell; ' ...
                     '$sc = $wsh.CreateShortcut(''%s''); ' ...
                     'Write-Output ($sc.TargetPath); ' ...
                     'Write-Output ($sc.Arguments);'], MatP.abqCmdShortcutPath);
    
    % Run PowerShell from MATLAB
    [status, result] = system(['powershell -NoProfile -ExecutionPolicy Bypass -Command "', psCmd, '"']);
    
    % format command for running in cmd
    parts = splitlines(strtrim(result));    % split output into lines
    target = strtrim(parts{1});             % first line = filepath
    args = '';
    if numel(parts) > 1
        args = strtrim(strjoin(parts(2:end),' '));  % remaining lines -> args
    end
    
    % remove anything from the first & onward (we don’t want cmd.exe /k part)
    args = regexprep(args, '\s*&.*$', '');
    
    % build final command string
    cmd = ['"' target '"'];                 % quote the filepath
    if ~isempty(args)
        cmd = [cmd ' ' args];               % append args if present
    end
%     fprintf('Command for cmd.exe:\n%s\n', cmd);

    % Create command to run job
    [~, UMATname, ext] = fileparts(MatP.UMATfilepath);
    jobCmd = ['abaqus job=Job-',MatP.unique,' user=',UMATname,ext];

    % run job
    cmd = [cmd, ' && ' jobCmd];
    system(cmd);
end

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
