 function JobSyubmitResults(fileID,folder,MatPunique)
    % submit the job 
    fprintf(fileID,'	mdb.jobs["%s"].submit(consistencyChecking=OFF) \n',MatPunique);
    % wait for the job to be complete 
    fprintf(fileID,'	mdb.jobs["%s"].waitForCompletion(); \n',MatPunique);
    % read output database to retrieve the J-integral values 
    fprintf(fileID,'	time.sleep(3); \n');
%     fprintf(fileID,'	odb = session.openOdb("%s.odb"); \n',MatPunique);
%     fprintf(fileID,'	timestep = odb.steps["Step-1"]; \n');
%     fprintf(fileID,'	alloutputs = timestep.historyRegions["ElementSet . ALL ELEMENTS"].historyOutputs; \n');
%     fprintf(fileID,'	Jval = []; K1val = []; K2val = []; \n');
%     fprintf(fileID,'	for i in alloutputs.keys(): \n');
%     fprintf(fileID,'		if i.split()[0] == "J": \n');
%     fprintf(fileID,'			Jval.append(alloutputs[i].data[-1][1]); \n');
%     fprintf(fileID,'		elif i.split()[0] == "K1": \n');
%     fprintf(fileID,'			K1val.append(alloutputs[i].data[-1][1]); \n');
%     fprintf(fileID,'		elif i.split()[0] == "K2": \n');
%     fprintf(fileID,'			K2val.append(alloutputs[i].data[-1][1]); \n');
%     % write values to the file 
%     fprintf(fileID,'	outfile = open(outputPth,"a"); \n');	
%     fprintf(fileID,'	for i in Jval: \n');
%     fprintf(fileID,'		outfile.write(str(i)); \n');
%     fprintf(fileID,'		outfile.write("\\t"); \n');
%     fprintf(fileID,'	outfile.write("\\n"); \n');	
%     fprintf(fileID,'	if(extractK):		 \n');
%     fprintf(fileID,'		for i in K1val: \n');
%     fprintf(fileID,'			outfile.write(str(i)); \n');
%     fprintf(fileID,'			outfile.write("\\t"); \n');
%     fprintf(fileID,'		outfile.write("\\n"); \n');		
%     fprintf(fileID,'		for i in K2val: \n');
%     fprintf(fileID,'			outfile.write(str(i)); \n');
%     fprintf(fileID,'			outfile.write("\\t"); \n');
%     fprintf(fileID,'		outfile.write("\\n"); \n');	
%     fprintf(fileID,'	outfile.write("\\n"); \n');
%     fprintf(fileID,'	outfile.close(); \n');	
%     fprintf(fileID,'	if counting!=len(masksVal): \n');
%     fprintf(fileID,'		odb.close(); \n');
%     fprintf(fileID,'		mdb.close() \n');
    % Creates a dummy file for Matlab to know that the analysis is over 
    SaveTmp = pythonFileName(fullfile(folder, [MatPunique 'Done.tmp']));
    fprintf(fileID,'finish_file = open("%s","a"); \n',SaveTmp);
    fprintf(fileID,'finish_file.write("DONE"); \n');
    fprintf(fileID,'finish_file.close(); \n');
 end  