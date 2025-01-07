function ApplyBCs(fileID,SaveModel)
%     fprintf(fileID,'outfile = open(outputPth,"w"); \n');
%     fprintf(fileID,'outfile.write("J-integral\\nK-I value\\nK-II value\\n"); \n');
%     fprintf(fileID,'outfile.close(); \n');
    fprintf(fileID,'counting = 1; \n');
    fprintf(fileID,'for diff_masks in masksVal: \n');
    fprintf(fileID,'	openMdb("%s"); \n',SaveModel);
    fprintf(fileID,'	allNodes = mdb.models["Model-1"].rootAssembly.instances["sample-1"].nodes \n');
    % Define crack in the FE model: 
    fprintf(fileID,'	crktip = allNodes.getByBoundingSphere((crackpoints[0][0],crackpoints[0][1],0), pow(10,-rndparam)); \n');
    fprintf(fileID,'	crkvec = allNodes.getByBoundingSphere((crackpoints[1][0],crackpoints[1][1],0), pow(10,-rndparam)); \n');
    fprintf(fileID,'	mdb.models["Model-1"].rootAssembly.engineeringFeatures.ContourIntegral( \n');
    fprintf(fileID,'		crackFront=Region(nodes=crktip), \n');
    fprintf(fileID,'		crackTip=Region(nodes=crktip), \n');
    fprintf(fileID,'		extensionDirectionMethod=Q_VECTORS	, name="Crack-1", qVectors=((crkvec[0],crktip[0]), )) \n');
    % get a correspondence node label/node coordinates 
    fprintf(fileID,'	aba_lab = []; aba_x = []; aba_y = []; \n');
    fprintf(fileID,'	for i in allNodes: \n');
    fprintf(fileID,'		aba_lab.append(i.label); \n');
    fprintf(fileID,'		aba_x.append(i.coordinates[0]); \n');
    fprintf(fileID,'		aba_y.append(i.coordinates[1]); \n');	
    fprintf(fileID,'	aba_x = [ round(elem, 8) for elem in aba_x ]; \n');
    fprintf(fileID,'	aba_y = [ round(elem, 8) for elem in aba_y ]; \n');	
    fprintf(fileID,'	aba_all = zip(aba_lab,aba_x,aba_y); \n');
    fprintf(fileID,'	del aba_lab; del aba_x; del aba_y; \n');
    fprintf(fileID,'	aba_sort = sorted(aba_all, key=itemgetter(1,2)); \n');
    fprintf(fileID,'	del aba_all; \n');	
    fprintf(fileID,'	bcarray = []; \n');
    fprintf(fileID,'	for i in oriclean: \n');
    fprintf(fileID,'		for j in aba_sort:				 \n');
    fprintf(fileID,'			if round(i[0], rndparam)==round(j[1], rndparam) and round(i[1], rndparam)==round(j[2], rndparam): \n');
    % test if the BC is in the free mask zone or not
    fprintf(fileID,'				if (j[1]>diff_masks[0] and j[1]<diff_masks[2] and j[2]>diff_masks[1] and j[2]<diff_masks[3]): \n');
    fprintf(fileID,'					continue; \n');
    % test if the BC is in the dangerZone
    fprintf(fileID,'				elif(j[1]>prim_oma_lim[0][0]+(unitsizex/2) and j[1]<prim_oma_lim[1][0]+(unitsizey/2) and j[2]>prim_oma_lim[0][1]-(unitsizex/2) and j[2]<prim_oma_lim[1][1]-(unitsizey/2)): \n');
    fprintf(fileID,'					continue; \n');
    fprintf(fileID,'				else: \n');
    fprintf(fileID,'					bcarray.append([j[0], i[2], i[3]]); \n');
    fprintf(fileID,'					aba_sort.remove(j); \n');
    fprintf(fileID,'					break; \n');
    % Find masked nodes labels and the elements including them and add them to be deleted
    fprintf(fileID,'	mskelelbl = []; \n');
    fprintf(fileID,'	for i in nodemsk: \n');
    fprintf(fileID,'		for j in aba_sort: \n');
    fprintf(fileID,'			if i[0]==j[1] and i[1]==j[2]: \n');
    % test if the node is in the dangerZone
    fprintf(fileID,'				if(j[1]>=prim_oma_lim[0][0] and j[1]<=prim_oma_lim[1][0] and j[2]>=prim_oma_lim[0][1] and j[2]<=prim_oma_lim[1][1]): \n');
    fprintf(fileID,'					continue; \n');
    fprintf(fileID,'				else: \n');
    % find constitutive elements 
    fprintf(fileID,'					curnodmsk = allNodes.getFromLabel(j[0]); \n');
    fprintf(fileID,'					curelemsk = curnodmsk.getElements(); \n');
    fprintf(fileID,'					for k in curelemsk: \n');
    fprintf(fileID,'						mskelelbl.append(k.label); \n');
    fprintf(fileID,'	mskelelbl = list(set(mskelelbl)); \n');				
    fprintf(fileID,'	for i in bcarray: \n');
    % find co-ordinates of the node 
    fprintf(fileID,'		myNodes = allNodes.sequenceFromLabels([i[0],]); \n');
    % create a BC at each node 
    fprintf(fileID,'		mdb.models["Model-1"].DisplacementBC(createStepName="Step-1", name="BC-"+str(i[0]), region=Region(nodes=myNodes),u1=i[1], u2=i[2], ur3=UNSET); \n');
end 