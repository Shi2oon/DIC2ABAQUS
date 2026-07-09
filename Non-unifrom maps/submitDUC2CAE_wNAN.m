function submitDUC2CAE_wNAN(fileID,folder,Unique)
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
FoldrOut = pythonFileName(folder);
fprintf(fileID,'os.chdir(r"%s");\n',FoldrOut);

fprintf(fileID,'mdb.ModelFromInputFile(inputFileName= \n');
SaveF = pythonFileName(fullfile(folder, [Unique '.inp']));
fprintf(fileID,'    "%s", name= \n', SaveF);
fprintf(fileID,'    "%s") \n', Unique);
fprintf(fileID,'del mdb.models["Model-1"] \n');
fprintf(fileID,'mdb.Job(atTime=None, contactPrint=OFF, description="", echoPrint=OFF,  \n');
fprintf(fileID,'    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,  \n');
fprintf(fileID,'    memory=90, memoryUnits=PERCENTAGE, model="%s", modelPrint=OFF,  \n', Unique);
fprintf(fileID,'    multiprocessingMode=DEFAULT, name="Job-%s", nodalOutputPrecision=SINGLE,  \n', Unique);
fprintf(fileID,'    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch="", type= \n');
fprintf(fileID,'    ANALYSIS, userSubroutine="", waitHours=0, waitMinutes=0) \n');
fprintf(fileID,'mdb.jobs["Job-%s"].submit(consistencyChecking=OFF) \n', Unique);
% wait for the job to be complete
% fprintf(fileID,'mdb.jobs["%s"].waitForCompletion(); \n',Unique);
% % read output database to retrieve the J-integral values
% fprintf(fileID,'time.sleep(3); \n');
% SaveTmp = pythonFileName(fullfile(folder, [Unique 'Done.tmp']));
% fprintf(fileID,'finish_file = open("Job-%s","a"); \n',SaveTmp);
% fprintf(fileID,'finish_file.write("DONE"); \n');
% fprintf(fileID,'finish_file.close(); \n');
end