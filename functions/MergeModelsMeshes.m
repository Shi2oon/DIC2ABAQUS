function MergeModelsMeshes(fileID,SaveModel)
    % Open saved model 
    fprintf(fileID,'openMdb("%s"); \n',SaveModel);
    % Read meshing data 
    fprintf(fileID,'pointer=open(file_out1,"r"); \n');
    fprintf(fileID,'mesh_nodes = []; \n');
    fprintf(fileID,'for line in iter(pointer): \n');
    fprintf(fileID,'	temp = line.split(); \n');
    fprintf(fileID,'	mesh_nodes.append((float(temp[0]),float(temp[1]))); \n');
    fprintf(fileID,'pointer.close(); \n');
    fprintf(fileID,'pointer=open(file_out2,"r"); \n');
    fprintf(fileID,'mesh_elem= []; \n');
    fprintf(fileID,'for line in iter(pointer): \n');
    fprintf(fileID,'	temp = line.split(); \n');
    fprintf(fileID,'	mesh_elem.append((int(temp[0]),int(temp[1]),int(temp[2]),int(temp[3]))); \n');
    fprintf(fileID,'pointer.close() \n');
    % merge the nodes, does not create a new node if it is on the border
    fprintf(fileID,'allNodes = mdb.models["Model-1"].rootAssembly.instances["sample-1"].nodes; \n');
    fprintf(fileID,'new_labels = {};lablist=[]; \n');
    fprintf(fileID,'for cur_nod in mesh_nodes: \n');
    % if the node belong to border, then it is already created; get its label 
    fprintf(fileID,'	if abs(cur_nod[0]-prim_oma_lim[0][0])<pow(10,-rndparam) or abs(cur_nod[0]-prim_oma_lim[1][0])<pow(10,-rndparam) or abs(cur_nod[1]-prim_oma_lim[1][1])<pow(10,-rndparam) or abs(cur_nod[1]-prim_oma_lim[0][1])<pow(10,-rndparam): \n');	
    fprintf(fileID,'		nodecur = allNodes.getByBoundingSphere(cur_nod+(0,), pow(10,-rndparam)); \n');
    fprintf(fileID,'		if(len(nodecur) == 0):	 \n');		
    fprintf(fileID,'			a = mdb.models["Model-1"].parts["sample"].Node((cur_nod+(0,)),None); \n');
    fprintf(fileID,'			new_labels[a.label] = [cur_nod[0],cur_nod[1]]; \n');
    fprintf(fileID,'			lablist.append(a.label); \n');
    fprintf(fileID,'		else: \n');
    fprintf(fileID,'			new_labels[nodecur[0].label] = [cur_nod[0],cur_nod[1]]; \n');
    fprintf(fileID,'			lablist.append(nodecur[0].label); \n');	
    % if no, create the node with a new label 
    fprintf(fileID,'	else:		 \n');
    fprintf(fileID,'		a = mdb.models["Model-1"].parts["sample"].Node((cur_nod+(0,)),None); \n');
    fprintf(fileID,'		new_labels[a.label] = [cur_nod[0],cur_nod[1]]; \n');
    fprintf(fileID,'		lablist.append(a.label); \n');
    % create the elements as specified in the file 
    fprintf(fileID,'allNodes = mdb.models["Model-1"].parts["sample"].nodes; \n');
    fprintf(fileID,'for cur_elem in mesh_elem:	 \n');
    fprintf(fileID,'	nodelem = (allNodes.getFromLabel(lablist[cur_elem[0]]),allNodes.getFromLabel(lablist[cur_elem[1]]),allNodes.getFromLabel(lablist[cur_elem[2]]), \n');
    fprintf(fileID,'		allNodes.getFromLabel(lablist[cur_elem[3]])); \n');
    fprintf(fileID,'	mdb.models["Model-1"].parts["sample"].Element(nodelem,QUAD4); \n');
end 