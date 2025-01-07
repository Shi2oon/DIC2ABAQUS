function [folder]=Adjust4Direction(saveTo,Angle)
% a code to calculate J, K-I, K-II for a specific gemoetry using declerared
% angles (Ang) as vector for virtual crack extension (q).
% The code takes:
% SaveTo: The .inp file which was used at least once for 1 q direction to
%           get the desired results.
% Stifness: anistropic stifness matrix or young modulus
% OfBy: for the units conversion to meter, 1 for meter, 1e-3 for mm, 1e-6 for um
% Ang: Range of desired angles for q, the range is + from 0 to 180 and
%       negative for 181 (-179) to 359 (-1)
% ToS: used if decleared angeles are reversed in Y direction.

% Written by Abdalrhaman Koko (abdo.koko@materials.ox.ac.uk) as part of 
% DPhil project
% last edit: 17 Dec. 2020 (Abdo)

    fid = fopen([erase(saveTo, '.inp') '.inp'], 'r');
    C = textscan(fid, '%s', 'Delimiter', '\n');     fclose(fid);
    S = regexp(C{1}{end-6}, ',', 'split');
    setN = erase(S{1},'_PickedSet');    setN = str2double(setN);
    if isnan(setN)
        setN = erase(S{1},'_PICKEDSET');    setN = str2double(setN);
    end
    C{1}{end-6} = ['_PickedSet' num2str(setN) ', _PickedSet' num2str(setN+1)...
        ', ' num2str(cosd(Angle)) ', ' num2str(sind(Angle))  ', 0.'];
    C{1}{end-1} = ['_PickedSet' num2str(setN) ', _PickedSet' num2str(setN+1)...
        ', ' num2str(cosd(Angle)) ', ' num2str(sind(Angle))  ', 0.'];
%     C{1}(end-5:end-1)=[];
    finalform = C{1}(1:end);
    iv=1;
    fileID = fopen([fileparts(saveTo) '\Job-' num2str(iv) '.inp'],'w');
    for i=1:size(finalform,1);          stri = finalform(i);
        if ~cellfun('isempty',stri);    fprintf(fileID,'%s\n',char(stri)); end
    end
    fclose(fileID); 
	folder = Js_Code(iv,fileparts(saveTo));

end

%%
function [OtF] = Js_Code(NumB,folder)
fileid = fullfile(folder,'Abaqus_script.py');
fileID = fopen(fileid,'w');
fprintf(fileID,'# -*- coding: mbcs -*- \n');
fprintf(fileID,'from part import * \n');
fprintf(fileID,'from material import * \n');
fprintf(fileID,'from section import * \n');
fprintf(fileID,'from assembly import * \n');
fprintf(fileID,'from step import * \n');
fprintf(fileID,'from interaction import * \n');
fprintf(fileID,'from load import * \n');
fprintf(fileID,'from mesh import * \n');
fprintf(fileID,'from optimization import * \n');
fprintf(fileID,'from job import * \n');
fprintf(fileID,'from sketch import * \n');
fprintf(fileID,'from visualization import * \n');
fprintf(fileID,'from connectorBehavior import * \n');
fprintf(fileID,'mdb.ModelFromInputFile(inputFileName= \n');
SaveTmp = pythonFileName([folder '/Job-' num2str(NumB) '.inp']);
fprintf(fileID,'    "%s",  \n',SaveTmp);
fprintf(fileID,'    name="Job-%d") \n',NumB);
fprintf(fileID,'from part import *\n');
fprintf(fileID,'from material import *\n');
fprintf(fileID,'from section import *\n');
fprintf(fileID,'from assembly import *\n');
fprintf(fileID,'from step import *\n');
fprintf(fileID,'from interaction import *\n');
fprintf(fileID,'from load import *\n');
fprintf(fileID,'from mesh import *\n');
fprintf(fileID,'from optimization import *\n');
fprintf(fileID,'from job import *\n');
fprintf(fileID,'from sketch import *\n');
fprintf(fileID,'from visualization import *\n');
fprintf(fileID,'from connectorBehavior import *\n');
fprintf(fileID,'mdb.Job(atTime=None, contactPrint=OFF, description="", echoPrint=OFF, \n');
fprintf(fileID,'    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, \n');
fprintf(fileID,'    memory=90, memoryUnits=PERCENTAGE, model="Job-%d", modelPrint=OFF, \n',NumB);
fprintf(fileID,'    multiprocessingMode=DEFAULT, name="Number_%d", nodalOutputPrecision=SINGLE, \n',NumB);
fprintf(fileID,'    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch="", type=\n');
fprintf(fileID,'    ANALYSIS, userSubroutine="", waitHours=0, waitMinutes=0)\n');
fprintf(fileID,'mdb.jobs["Number_%d"].submit(consistencyChecking=OFF)\n',NumB);
%{
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STARTED, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "clientHost": "OUMS-TJM10", "handle": 0, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(WARNING, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "message": "THE REQUEST FOR MISES OUTPUT WILL BE REPLACED BY A REQUEST FOR S OUTPUT", \n');
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(ODB_FILE, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "file": "A:\\OneDrive - Nexus365\\Work\\ABAQUS\\Python and ABAQUS\\Dual\\Number_%d.odb", \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(COMPLETED, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "message": "Analysis phase complete", "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STARTED, {"phase": STANDARD_PHASE, \n',NumB);
fprintf(fileID,'    "clientHost": "OUMS-TJM10", "handle": 14604, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STEP, {"phase": STANDARD_PHASE, "stepId": 1, \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(ODB_FRAME, {"phase": STANDARD_PHASE, "step": 0, \n',NumB);
fprintf(fileID,'    "frame": 0, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(MEMORY_ESTIMATE, {"phase": STANDARD_PHASE, \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d", "memory": 35.0})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(ODB_FRAME, {"phase": STANDARD_PHASE, "step": 0, \n',NumB);
fprintf(fileID,'    "frame": 1, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STATUS, {"totalTime": 1.0, "attempts": 1, \n',NumB);
fprintf(fileID,'    "timeIncrement": 1.0, "increment": 1, "stepTime": 1.0, "step": 1, \n');
fprintf(fileID,'    "jobName": "Number_%d", "severe": 0, "iterations": 1, \n',NumB);
fprintf(fileID,'    "phase": STANDARD_PHASE, "equilibrium": 1})\n');
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(END_STEP, {"phase": STANDARD_PHASE, "stepId": 1, \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(COMPLETED, {"phase": STANDARD_PHASE, \n',NumB);
fprintf(fileID,'    "message": "Analysis phase complete", "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'# Save by scro3511 on 2020_05_10-17.16.29; build 2016 2015_09_24-21.31.09 126547\n');
%}
fclose(fileID);

%% Excute
PWD =pwd;           cd(folder)
system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
cd(PWD);
OtF = [folder '\Number_' num2str(NumB)];
end