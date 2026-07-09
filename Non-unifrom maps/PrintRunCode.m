function CAE = PrintRunCode(MatP,folder)
fileid = fullfile(folder, 'Abaqus_script.py');
fileID = fopen(fileid,'w');
submitDUC2CAE_wNAN(fileID,folder,MatP.unique,MatP.type, MatP)
% printTiffs(fileID,folder,'N',MatP.unique,2);   %% save tiff images for U and S
% printTiffs(fileID,folder,'D',MatP.unique,2);   %% save tiff images for deformaed U and S
fclose(fileID);
CAE = fullfile(folder,['Job-' MatP.unique]);
%%
PWD =pwd;           cd(folder)

if MatP.type ~= 'U'
    system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
    % % unix(['abaqus cae ','noGUI','=Abaqus_script.py']);   % Unix system
    % Pause to ensure the script has time to execute
    pause(10); % Adjust the pause duration as needed
    % Close the Abaqus CAE GUI
[status, cmdout] = system('tasklist /FI "IMAGENAME eq abq*.*" /FO CSV');
end

if MatP.type == 'U'
%     status = runAbaqus(MatP.runMode, MatP.unique,pwd,MatP.abqCmdShortcutPath);
    %
    % PowerShell script to find ifort target and argument from abaqus shortcut
    psCmd = sprintf(['$wsh = New-Object -ComObject WScript.Shell; ' ...
                     '$sc = $wsh.CreateShortcut(''%s''); ' ...
                     'Write-Output ($sc.TargetPath); ' ...
                     'Write-Output ($sc.Arguments);'], MatP.abqCmdShortcutPath);
    
    % Run PowerShell from MATLAB
    [status, result] = system(['powershell -NoProfile -ExecutionPolicy Bypass -Command "', psCmd, '"']);
    
%     if status == 0
%         disp('Code is not running, just open abaqus and run it yourself!')
%         return
%     end
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

    system([cmd,' && ','abaqus cae ','noGUI','=Abaqus_script.py'])
    pause(10)
% Close the Abaqus CAE GUI
[status, cmdout] = system('tasklist /FI "IMAGENAME eq abq*.*" /FO CSV');
    %}
end

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

    % Create a timer object
%     t = timer('StartDelay', 10, ...
%               'TimerFcn', @(~,~)closeDialog('YesNoCancel'));   start(t);
choice = 'No';%questdlg('Delete job files except the ODB file? (Auto-closes in 10s)', ...
%                   'Timeout Option', ...
%                   'Yes','No','No');

if strcmp(choice,'Yes')
    delete([folder '\Job-' MatP.unique '.log']);
    delete([folder '\Job-' MatP.unique '.msg']);
    delete([folder '\Job-' MatP.unique '.prt']);
    delete([folder '\Job-' MatP.unique '.sim']);
    delete([folder '\Job-' MatP.unique '.sta']);
    delete([folder '\Job-' MatP.unique '.com']);
    delete([folder '\Job-' MatP.unique 'Done.tmp']);
    delete(fullfile(folder,'*.for')); % delete any umats
end

if ~exist([MatP.unique '.rpt'],'file')
    fprintf('\nYou can find a python script to run in Abaqus in\n%s\n',folder);
end
cd(PWD);

end

% Callback function for timer to close the dialog
function closeDialog(dlgName)
    % Find all questdlg figures
    fig = findall(0, 'Type', 'figure', 'Name', dlgName);
    if ~isempty(fig)
        delete(fig); % Closes the figure
    end
end

function status = runAbaqus(runMode, trialDir,job,abqCmdShortcutPath)
cmdInner = sprintf('cd /d "%s" && abaqus job="%s" input="%s" user="%s" interactive',trialDir,job,[job '.inp'],"umat.for");
[status,out] = runCmd(runMode,cmdInner,trialDir,abqCmdShortcutPath);
% writeText(fullfile(trialDir,[job '_abaqus_command.txt']),sprintf('COMMAND:\n%s\n\nOUTPUT:\n%s\n',cmdInner,out));
% ok = (status==0);
end

function [status,out] = runCmd(runMode,cmdInner,workDir,abqCmdShortcutPath)
if strcmpi(runMode,'plain')
    [status,out] = system(cmdInner); return;
elseif strcmpi(runMode,'shortcut')
    bat = fullfile(workDir,'run_abaqus_from_shortcut.bat');
%     writeText(bat, sprintf('@echo on\n%s\nexit /b %%ERRORLEVEL%%\n',cmdInner));
    ps = sprintf('powershell -NoProfile -ExecutionPolicy Bypass -Command "Start-Process -FilePath ''%s'' -ArgumentList ''/c \"\"%s\"\"'' -Wait -PassThru | %% { exit $_.ExitCode }"',abqCmdShortcutPath,bat);
    [status,out] = system(ps); return;
else
    error('Unknown runMode: %s',runMode);
end
end

% function writeText(f,s), fid=fopen(f,'w'); if fid<0, error('Cannot write %s',f); end; fwrite(fid,s); fclose(fid); end
