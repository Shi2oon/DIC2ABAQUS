function RemeshCrackedRegion(fileID)
    % Create deleted part as a rectangle 
    fprintf(fileID,'mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=200.0); \n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].rectangle(point1=prim_oma_lim[0], point2=prim_oma_lim[1]); \n');
    fprintf(fileID,'mdb.models["Model-1"].Part(dimensionality=TWO_D_PLANAR, name="temp", type=DEFORMABLE_BODY); \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].BaseShell(sketch=mdb.models["Model-1"].sketches["__profile__"]); \n');
    fprintf(fileID,'del mdb.models["Model-1"].sketches["__profile__"]; \n');
    % create crack points as Datum points
    fprintf(fileID,'dt_pts = findDatumCrack(crackpoints,prim_oma_lim,unitsizex,rndparam); \n');
    fprintf(fileID,'plandef_dat = []; \n');
    fprintf(fileID,'for i in dt_pts: \n');
    fprintf(fileID,'	plandef_dat.append(mdb.models["Model-1"].parts["temp"].DatumPointByCoordinate((i[0],i[1],0))); \n');
    % partition part using the datum points 
    fprintf(fileID,'for segnum in range(len(plandef_dat)-1): \n');
    fprintf(fileID,'	mdb.models["Model-1"].parts["temp"].PartitionFaceByShortestPath(faces=mdb.models["Model-1"].parts["temp"].faces[0],  \n');
    fprintf(fileID,'	point1=mdb.models["Model-1"].parts["temp"].datums[plandef_dat[segnum].id], \n');
    fprintf(fileID,'	point2=mdb.models["Model-1"].parts["temp"].datums[plandef_dat[segnum+1].id]); \n');
    % create independent instance 
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.Instance(dependent=OFF, name="temp-1", part= mdb.models["Model-1"].parts["temp"]); \n');
    % select crack edge
    fprintf(fileID,'crkedg=[]; \n');
    fprintf(fileID,'vertlist = mdb.models["Model-1"].rootAssembly.instances["temp-1"].vertices; \n');
    fprintf(fileID,'for i in range(1,len(dt_pts)-1): \n');
    fprintf(fileID,'	curpt = vertlist.getByBoundingSphere([dt_pts[i][0],dt_pts[i][1],0],0.1); \n');
    fprintf(fileID,'	for kk in curpt[0].getEdges(): \n');
    fprintf(fileID,'		crkedg.append(kk); \n');
    fprintf(fileID,'crkedg = list(set(crkedg)) \n');
    % all other edges have to be seeded 
    fprintf(fileID,'alledges = mdb.models["Model-1"].rootAssembly.instances["temp-1"].edges; \n');
    fprintf(fileID,'for i in range(len(alledges)): \n');
    fprintf(fileID,'	if i not in crkedg: \n');
    fprintf(fileID,'		mdb.models["Model-1"].rootAssembly.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1,edges=(alledges[i],), size=unitsizex); \n');
    % apply crack seam 
    fprintf(fileID,'crkdec = {} \n');
    fprintf(fileID,'for i in crackpoints: \n');
    fprintf(fileID,'	curpt = vertlist.getByBoundingSphere([i[0],i[1],0],0.1); \n');
    fprintf(fileID,'	for kk in curpt[0].getEdges(): \n');
    fprintf(fileID,'		if kk in crkdec: \n');
    fprintf(fileID,'			crkdec[kk]+=1; \n');
    fprintf(fileID,'		else: \n');
    fprintf(fileID,'			crkdec[kk] = 1; \n');
    fprintf(fileID,'for i in crkdec.items(): \n');
    fprintf(fileID,'	if(i[1]>1): \n');
    fprintf(fileID,'		mdb.models["Model-1"].rootAssembly.Set(edges=alledges[i[0]:i[0]+1], name="SetCrk"+str(i[0])) \n');
    fprintf(fileID,'		mdb.models["Model-1"].rootAssembly.engineeringFeatures.assignSeam(regions=mdb.models["Model-1"].rootAssembly.sets["SetCrk"+str(i[0])]) \n');
    % apply the mesh type and mesh
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.setMeshControls(elemShape=QUAD, regions= \n');
    fprintf(fileID,'	mdb.models["Model-1"].rootAssembly.instances["temp-1"].faces) \n');
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.generateMesh(regions=( \n');
    fprintf(fileID,'	mdb.models["Model-1"].rootAssembly.instances["temp-1"], )) \n');
 end