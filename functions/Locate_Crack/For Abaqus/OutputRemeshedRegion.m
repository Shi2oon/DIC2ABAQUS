function OutputRemeshedRegion(fileID)
    % extract nodes and elements for the next step and write them to a file
    fprintf(fileID,'file_out1 = "%s"; \n','C:\Temp\meshoutnod.OMA');
    fprintf(fileID,'file_out2 = "%s"; \n','C:\Temp\meshoutelem.OMA');
    fprintf(fileID,'allNodes = mdb.models["Model-1"].rootAssembly.instances["temp-1"].nodes \n');
    fprintf(fileID,'allElems = mdb.models["Model-1"].rootAssembly.instances["temp-1"].elements \n');
    fprintf(fileID,'coor_spyx_sec=[]; coor_spyy_sec=[]; \n');
    fprintf(fileID,'pointer=open(file_out1,"w"); \n');
    fprintf(fileID,'for i in allNodes : \n');
    fprintf(fileID,'	pointer.write(str(i.coordinates[0])+" "+str(i.coordinates[1])+"\\n"); \n');
    fprintf(fileID,'	coor_spyx_sec.append(i.coordinates[0]); coor_spyy_sec.append(i.coordinates[1]); \n');
    fprintf(fileID,'pointer.close(); \n');
    fprintf(fileID,'coor_spyx_sec = list(set(coor_spyx_sec)); coor_spyy_sec = list(set(coor_spyy_sec)); \n');
    fprintf(fileID,'sec_oma_lim = [(min(coor_spyx_sec),min(coor_spyy_sec)),(max(coor_spyx_sec),max(coor_spyy_sec))]; \n');	
    fprintf(fileID,'pointer=open(file_out2,"w"); \n');
    fprintf(fileID,'for i in allElems : \n');
    fprintf(fileID,'	pointer.write(str(i.connectivity[0])+" "+str(i.connectivity[1])+" "+str(i.connectivity[2])+" "+str(i.connectivity[3])+"\\n"); \n');
    fprintf(fileID,'pointer.close(); \n');
    fprintf(fileID,'mdb.close(); \n');
end