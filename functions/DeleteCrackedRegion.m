function SaveModel = DeleteCrackedRegion(fileID,folder,MatPunique)
    % get a list of the existing nodes \n');
    fprintf(fileID,'allNodes = mdb.models["Model-1"].rootAssembly.instances["sample-1"].nodes \n');
    % get the limit coordinates of the nodes affected by OMA
    fprintf(fileID,'oma_lim = [];	 \n');
    fprintf(fileID,'oma_tru_lim = getOMAlimit(crackpoints,dangerZone*unitsizex); \n');
    fprintf(fileID,'oma_lim.append( (min((i for i in x), key=lambda var1:abs(var1-oma_tru_lim[0][0])),min((i for i in y), key=lambda var1:abs(var1-oma_tru_lim[0][1])))); \n');
    fprintf(fileID,'oma_lim.append( (min((i for i in x), key=lambda var1:abs(var1-oma_tru_lim[1][0])),min((i for i in y), key=lambda var1:abs(var1-oma_tru_lim[1][1])))); \n');
    % get a list of the mesh nodes and elements concerned by OMA
    fprintf(fileID,'oma_elem ={}; \n');
    fprintf(fileID,'coor_spyx = []; coor_spyy = []; \n');
    fprintf(fileID,'for i in allNodes: \n');
    fprintf(fileID,'	if round(i.coordinates[0],rndparam)>=round(oma_lim[0][0],rndparam) and round(i.coordinates[0],rndparam)<=round(oma_lim[1][0],rndparam) and round(i.coordinates[1],rndparam)>=round(oma_lim[0][1],rndparam) and round(i.coordinates[1],rndparam)<=round(oma_lim[1][1],rndparam): \n');
    fprintf(fileID,'		for j in i.getElements(): \n');
    fprintf(fileID,'			if j.label in oma_elem: \n');
    fprintf(fileID,'				oma_elem[j.label] += 1; \n');
    fprintf(fileID,'			else: \n');
    fprintf(fileID,'				oma_elem[j.label] = 1; \n');
    % spy the coordinates to get the exact sides bounding box (i.e. with Abaq shity rounding)
    fprintf(fileID,'		coor_spyx.append(i.coordinates[0]); coor_spyy.append(i.coordinates[1]); \n');
    % tidy up the exact BBox data 
    fprintf(fileID,'coor_spyx = list(set(coor_spyx)); coor_spyy = list(set(coor_spyy)); \n');
    fprintf(fileID,'prim_oma_lim = [(min(coor_spyx),min(coor_spyy)),(max(coor_spyx),max(coor_spyy))]; \n');
    % get a list of the elements to delete 
    fprintf(fileID,'elemtodel = []; \n');
    fprintf(fileID,'for i in oma_elem.keys(): \n');
    fprintf(fileID,'	if(oma_elem[i]==4): \n');
    fprintf(fileID,'		elemtodel.append(i); \n');
    % and delete elements 
    fprintf(fileID,'mdb.models["Model-1"].parts["sample"].deleteElement(elements=mdb.models["Model-1"].parts["sample"].elements.sequenceFromLabels(elemtodel),deleteUnreferencedNodes=ON); \n');
    % save the model 
    SaveModel = pythonFileName(fullfile(folder, [MatPunique '.cae']));
    fprintf(fileID,'mdb.saveAs("%s"); \n',SaveModel);
    fprintf(fileID,'mdb.close() \n');
end