function Print4OUROMA(MatP, maskdim,SaveD,len)
[folder,~,~] = fileparts(SaveD);
fileid = fullfile(folder, 'Ur_Code.py');
fileID = fopen(fileid,'w');

%% prepare used modules
    fprintf(fileID,'from __future__ import division  \n');
    % ABAQUS
    fprintf(fileID,'from abaqusConstants import * \n');
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
    fprintf(fileID,'from odbAccess import * \n');
    fprintf(fileID,'\n');
    % NOT ABAQUS BUT EXISTING LIBS 
    fprintf(fileID,'from operator import itemgetter \n');
    fprintf(fileID,'from math import ceil \n');
    fprintf(fileID,'import numpy  as np \n');
    fprintf(fileID,'import os \n');
    fprintf(fileID,'\n');

%% Materials and Map definiation 
fprintf(fileID,'#Set Abaqus Work Directory to the data file  \n');
fprintf(fileID,'# An error at ALL ELLEMENTS indicate you need to re-define \n');
fprintf(fileID,'# the crak vector Change and submit the job again\n');
dicInputPath = pythonFileName(SaveD);
fprintf(fileID,'dicInputPath = "%s"; \n',dicInputPath);
% ATTENTION : crackpoints should be defined with the crack tip in first position.....
% If a crackpoint is on a border, it's coordinate has to be EXACTLY the 
% coordinate of the border defined in the DICinput file
fprintf(fileID,'crackpoints = ((%f,%f),(%f,%f));\n',....
                maskdim.xo(1),maskdim.yo(1),maskdim.xo(2),maskdim.yo(2));
% Area without boundary conditions;masksVal are Xmin Ymin Xmax Ymax
fprintf(fileID,'masksVal=((%f,%f,%f,%f),);\n',....
                maskdim.xm(2),maskdim.ym(2),maskdim.xm(1),maskdim.ym(1));
%Nb of subsets masked around the crack
fprintf(fileID,'dangerZone = %d;\n', abs(maskdim.ds1-maskdim.ds2)); 
% NB of J-integral contours
fprintf(fileID,'nbCtrJint  = %d;\n', round(len*0.47-1,0));
% Where to write the results
outputPth = pythonFileName(fullfile(folder, 'KJ_Output.txt'));
fprintf(fileID,'outputPth  = "%s";\n', outputPth);

% material's parameters Si (mm) units, Elastic modulus Pa, Poisson's ratio
% Yield stress Pa, Exponent, Yield offset
fprintf(fileID,'MaterialName = "%s";\n',MatP.Mat);

if      MatP.type == 'E'
        fprintf(fileID,'matLaw = "Elastic";\n');
% Extract K or not using the interaction integral method / Only valid for 'Elastic Model'
        fprintf(fileID,'extractK   = 1;\n');   
        fprintf(fileID,'matParams=((%d,%.3f),);\n',MatP.E,MatP.nu);
elseif  MatP.type == 'R'
        fprintf(fileID,'matLaw = "Ramberg-Osgood";\n');
        fprintf(fileID,'matParams=((%d,%.3f,%d,%.3f,%.3f),);\n',...
                MatP.E,MatP.nu,MatP.yield,MatP.Exponent, MatP.Yield_offset);
        fprintf(fileID,'extractK   = 0;\n'); 
elseif  MatP.type == 'A'
        fprintf(fileID,'matLaw = "Elastic-Anisotropic";\n');
        fprintf(fileID,'matParams=((%d,%.3f),);\n',MatP.E*1e9,MatP.nu);
        fprintf(fileID,'extractK   = 1;\n'); 
        C = MatP.Maps.Stiffness.*1e9;
%stiffness tensors [Pa]        
fprintf(fileID,'# You will need to define local coordinate to do so\n');
fprintf(fileID,'# Select the Proprety Module, then in the toolbar select Assign\n');
fprintf(fileID,'# Material Orientation, select all the Module, use Default CSYS, OK \n');
fprintf(fileID,'C11 = %d;\tC12 = %d;\tC13 = %d;\n',C(1,1),C(1,2),C(1,2));		
fprintf(fileID,'C22 = %d;\tC23 = %d;\tC33 = %d;\n',C(2,2),C(2,3),C(3,3));
fprintf(fileID,'C44 = %d;\tC55 = %d;\tC66 = %d;\n',C(4,4),C(5,5),C(6,6));
end

%% define python functions
    % DETERMINE RECTANGULAR CRACK BOUNDING BOX
    fprintf(fileID,'def getOMAlimit(crackp,offset): \n');
    fprintf(fileID,'	xmin = min(crackp,key=itemgetter(0))[0] - offset; \n');
    fprintf(fileID,'	ymin = min(crackp,key=itemgetter(1))[1] - offset; \n');
    fprintf(fileID,'	xmax = max(crackp,key=itemgetter(0))[0] + offset; \n');
    fprintf(fileID,'	ymax = max(crackp,key=itemgetter(1))[1] + offset; \n');
    fprintf(fileID,'	return ((xmin,ymin),(xmax,ymax)); \n');
    % DEFINE DATUM POINTS COORDINATES TO DEFINE THE CRACK
    fprintf(fileID,'def findDatumCrack(crackpts,oma_lim,lc,rndparam): \n');
    % bounding box edge labels \n');
    fprintf(fileID,'	ptidx = np.array([[1,4],[3,2]]); \n');
    % Find closest edge of the first and last crack seam \n');
    % points on the rectangle \n');
    fprintf(fileID,'	oma_lim = np.array(oma_lim); \n');
    fprintf(fileID,'	firS = crackpts[0]; \n');
    fprintf(fileID,'	lasS = crackpts[len(crackpts)-1]; \n');
    fprintf(fileID,'	disFi = abs(oma_lim-firS); \n');
    fprintf(fileID,'	disLa = abs(oma_lim-lasS); \n');
    fprintf(fileID,'	firsPProj = ptidx[np.where(disFi==np.amin(disFi))][0]; \n');
    fprintf(fileID,'	lastPProj = ptidx[np.where(disLa==np.amin(disLa))]; \n');
    fprintf(fileID,'	frealProj=[]; \n');
    fprintf(fileID,'	if (np.amin(disFi)<pow(10,-rndparam)): \n');
    fprintf(fileID,'		frealProj = firS; \n');
    fprintf(fileID,'	else: \n');
    fprintf(fileID,'		if firsPProj == 1: \n');
    fprintf(fileID,'			ycoo = firS[1]+(oma_lim[0][1]%%lc-firS[1]%%lc); \n');
    fprintf(fileID,'			frealProj = [oma_lim[0][0],ycoo]; \n');
    fprintf(fileID,'		elif firsPProj == 2: \n');
    fprintf(fileID,'			ycoo = firS[0]+(oma_lim[0][0]%%lc-firS[0]%%lc); \n');
    fprintf(fileID,'			frealProj = [ycoo,oma_lim[1][1]]; \n');
    fprintf(fileID,'		elif firsPProj == 3: \n');
    fprintf(fileID,'			ycoo = firS[1]+(oma_lim[0][1]%%lc-firS[1]%%lc); \n');
    fprintf(fileID,'			frealProj = [oma_lim[1][0],ycoo]; \n');
    fprintf(fileID,'		else: \n');
    fprintf(fileID,'			ycoo = firS[0]+(oma_lim[0][0]%%lc-firS[0]%%lc); \n');
    fprintf(fileID,'			frealProj = [ycoo,oma_lim[0][1]]; \n');
    fprintf(fileID,'	lrealProj=[]; \n');
    fprintf(fileID,'	if (np.amin(disLa)<pow(10,-rndparam)): \n');
    fprintf(fileID,'		lrealProj=lasS; \n');
    fprintf(fileID,'	else: \n');
    fprintf(fileID,'		if lastPProj == 1: \n');
    fprintf(fileID,'			ycoo = lasS[1]+(oma_lim[0][1]%%lc-lasS[1]%%lc); \n');
    fprintf(fileID,'			lrealProj = [oma_lim[0][0],ycoo]; \n');
    fprintf(fileID,'		elif lastPProj == 2: \n');
    fprintf(fileID,'			ycoo = lasS[0]+(oma_lim[0][0]%%lc-lasS[0]%%lc); \n');
    fprintf(fileID,'			lrealProj = [ycoo,oma_lim[1][1]]; \n');
    fprintf(fileID,'		elif lastPProj == 3: \n');
    fprintf(fileID,'			ycoo = lasS[1]+(oma_lim[0][1]%%lc-lasS[1]%%lc); \n');
    fprintf(fileID,'			lrealProj = [oma_lim[1][0],ycoo]; \n');
    fprintf(fileID,'		else: \n');
    fprintf(fileID,'			ycoo = lasS[0]+(oma_lim[0][0]%%lc-lasS[0]%%lc); \n');
    fprintf(fileID,'			lrealProj = [ycoo,oma_lim[0][1]]; \n');
    fprintf(fileID,'	datumPoints = []; \n');
    fprintf(fileID,'	datumPoints.append(tuple(frealProj)); \n');
    fprintf(fileID,'	for pc in crackpts: \n');
    fprintf(fileID,'		datumPoints.append(pc); \n');
    fprintf(fileID,'	datumPoints.append(tuple(lrealProj)); \n');
    fprintf(fileID,'	return f7(datumPoints); \n');
    fprintf(fileID,'def f7(seq): \n');
    fprintf(fileID,'    seen = set(); \n');
    fprintf(fileID,'    seen_add = seen.add; \n');
    fprintf(fileID,'    return [ x for x in seq if not (x in seen or seen_add(x))]; \n');
    
%% Reading input file
    fprintf(fileID,'x=[]; y=[]; \n');
    fprintf(fileID,'vx=[]; vy=[]; \n');
    fprintf(fileID,'x0=[]; y0=[]; \n');
    % Open input file and read nodes positions and displacements
    fprintf(fileID,'pointer=open(dicInputPath,"r") \n');
    fprintf(fileID,'for line in iter(pointer): \n');
    fprintf(fileID,'	temp = line.split() \n');
    fprintf(fileID,'	try: \n');
    fprintf(fileID,'		x.append(float(temp[0])) \n');
    fprintf(fileID,'	except: \n');
    fprintf(fileID,'		continue \n');
    fprintf(fileID,'	y.append(float(temp[1])) \n');
    fprintf(fileID,'	vx.append(float(temp[2])) \n');
    fprintf(fileID,'	vy.append(float(temp[3])) \n');
    fprintf(fileID,'pointer.close() \n');
    % Determine number of nodes in each direction, mesh is supposed quad regular 
    fprintf(fileID,'nonodesx=len(set(x));				nonodesy=len(set(y)); \n');
    % Round nodes position to get rid off machine error 
    fprintf(fileID,'x = [ round(elem, 8) for elem in x ] \n');
    fprintf(fileID,'y = [ round(elem, 8) for elem in y ] \n');
    % Determine mesh size in each direction 
    fprintf(fileID,'unitsizex=abs((max(x)-min(x))/(nonodesx-1)) \n');
    fprintf(fileID,'unitsizey=abs((max(y)-min(y))/(nonodesy-1)) \n');
    % Determine the rounding order to apply later; it is define as 5%% of the unitsize 
    fprintf(fileID,'rndparam = -int(floor((log10(0.05*unitsizex)))); \n');
    
%% ABAQUS PART CREATION & PRELIM MESHING
    % Sketching and creating the part as a rectangle
    fprintf(fileID,'mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=200.0); \n');
    fprintf(fileID,'mdb.models["Model-1"].sketches["__profile__"].rectangle(point1=(min(x), max(y)), point2=(max(x), min(y))); \n');
    fprintf(fileID,'mdb.models["Model-1"].Part(dimensionality=TWO_D_PLANAR, name="temp", type=DEFORMABLE_BODY); \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].BaseShell(sketch=mdb.models["Model-1"].sketches["__profile__"]); \n');
    fprintf(fileID,'del mdb.models["Model-1"].sketches["__profile__"]; \n');
    % Selecting the 4 edges and meshing them 
    fprintf(fileID,'edge1=mdb.models["Model-1"].parts["temp"].edges.findAt([(max(x)+min(x))/2,max(y),0]) \n');
    fprintf(fileID,'edge2=mdb.models["Model-1"].parts["temp"].edges.findAt([max(x),(max(y)+min(y))/2,0]) \n');
    fprintf(fileID,'edge3=mdb.models["Model-1"].parts["temp"].edges.findAt([(max(x)+min(x))/2,min(y),0]) \n');
    fprintf(fileID,'edge4=mdb.models["Model-1"].parts["temp"].edges.findAt([min(x),(max(y)+min(y))/2,0]) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge1],number=(nonodesx-1),constraint=FIXED) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge3],number=(nonodesx-1),constraint=FIXED) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge2],number=(nonodesy-1),constraint=FIXED) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge4],number=(nonodesy-1),constraint=FIXED) \n');
    % Defining mesh type, elements types and mesh the full part 
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].setMeshControls(elemShape=QUAD, regions= mdb.models["Model-1"].parts["temp"].faces, technique=STRUCTURED) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].setElementType(elemTypes=(ElemType(elemCode=CPS4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, hourglassControl=DEFAULT, distortionControl=DEFAULT), ), regions=(mdb.models["Model-1"].parts["temp"].faces, )) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].generateMesh(); \n');
    % Creating orphan mesh
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].PartFromMesh(name="sample") \n');
    fprintf(fileID,'del mdb.models["Model-1"].parts["temp"]; \n');
    % Creating calculation step; axis system and assembly instance 
    fprintf(fileID,'mdb.models["Model-1"].StaticStep(name="Step-1", previous="Initial") \n');
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.DatumCsysByDefault(CARTESIAN) \n');
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.Instance(dependent=ON, name="sample-1", part=mdb.models["Model-1"].parts["sample"]) \n');

%% DELETE CRACKED REGION
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
    fprintf(fileID,'mdb.saveAs("C:\\Temp\\Step1OMA.cae"); \n');
    fprintf(fileID,'mdb.close() \n');

%% REMESH CRACKED REGION
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
    
%% OUTPUT REMESHED REGION
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
    
%% MERGE MODELS MESHES
    % Open saved model 
    fprintf(fileID,'openMdb("C:\\Temp\\Step1OMA.cae"); \n');
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
    
%% JOB, MATERIAL, HISTORYOUT
    % creating job
    fprintf(fileID,'mdb.Job(atTime=None, contactPrint=OFF, description="", echoPrint=OFF,  \n');
    fprintf(fileID,'	explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,  \n');
    fprintf(fileID,'	memory=90, memoryUnits=PERCENTAGE, model="Model-1", modelPrint=OFF, \n'); 
    fprintf(fileID,'	multiprocessingMode=DEFAULT, name="final", nodalOutputPrecision=SINGLE,  \n');
    fprintf(fileID,'	numCpus=1, queue=None, scratch="", type=ANALYSIS, userSubroutine="",  \n');
    fprintf(fileID,'	waitHours=0, waitMinutes=0) \n');
    % Creating the material 
    fprintf(fileID,'mdb.models["Model-1"].Material(name=MaterialName) \n');
    fprintf(fileID,'if matLaw=="Elastic": \n');
    fprintf(fileID,'	mdb.models["Model-1"].materials[MaterialName].Elastic(table=((matParams[0][0], matParams[0][1]), )) \n');
    fprintf(fileID,'elif matLaw=="Ramberg-Osgood": \n');
    fprintf(fileID,'	mdb.models["Model-1"].materials[MaterialName].DeformationPlasticity(table=((matParams[0][0],  \n');
    fprintf(fileID,'		matParams[0][1], matParams[0][2], matParams[0][3], matParams[0][4]), )) \n');
    fprintf(fileID,'elif matLaw == "Elastic-Anisotropic": \n');
    fprintf(fileID,'	mdb.models["Model-1"].materials[MaterialName].Elastic(table=((C11, C12, C22, C13, C23, C33,  \n');
    fprintf(fileID,'		0.0, 0.0, 0.0, C44, 0.0, 0.0, 0.0, 0.0, C55, 0.0, 0.0,  \n');
    fprintf(fileID,'		0.0, 0.0, 0.0, C66), ), type=ANISOTROPIC) \n');
    fprintf(fileID,'mdb.models["Model-1"].HomogeneousSolidSection(material=MaterialName, name="Section-1", thickness=None) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["sample"].SectionAssignment(offset=0.0, offsetField="", offsetType=MIDDLE_SURFACE, region=Region( \n');
    fprintf(fileID,'	elements=mdb.models["Model-1"].parts["sample"].elements), sectionName="Section-1",  \n');
    fprintf(fileID,'	thicknessAssignment=FROM_SECTION) \n');
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.regenerate() \n');
    fprintf(fileID,'mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1",  \n');
    fprintf(fileID,'	createStepName="Step-1", name="OutpJint", numberOfContours=nbCtrJint, rebar= \n');
    fprintf(fileID,'	EXCLUDE, sectionPoints=DEFAULT) \n');
    fprintf(fileID,'if(extractK): \n');
    fprintf(fileID,'	mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1",  \n');
    fprintf(fileID,'	createStepName="Step-1", name="OutpKval", numberOfContours=nbCtrJint, rebar= \n');
    fprintf(fileID,'	EXCLUDE, contourType=K_FACTORS, sectionPoints=DEFAULT) \n');
    
%% PREPARE FOR BOUNDARY CONDITIONS
    % sort the original bc the same way
    fprintf(fileID,'ori_all = zip(x,y,vx,vy); \n');
    fprintf(fileID,'ori_sort = sorted(ori_all, key=itemgetter(0,1)); \n');
    fprintf(fileID,'del ori_all; \n');
    fprintf(fileID,'oriclean=[]; nodemsk = []; \n');
    fprintf(fileID,'for i in ori_sort: \n');
    % if the BC is not null then it is not masked in DIC so we can apply it 
    fprintf(fileID,'	if i[2]!=0 and i[3]!=0: \n');
    fprintf(fileID,'		oriclean.append(i); \n');
    fprintf(fileID,'	else: \n');
    fprintf(fileID,'		nodemsk.append(i); \n');
    fprintf(fileID,'del ori_sort; \n');
    fprintf(fileID,'mdb.saveAs("C:\\Temp\\Step1OMA.cae"); \n');
    fprintf(fileID,'mdb.close() \n');
    
%% APPLY BOUNDARY CONDITIONS
    fprintf(fileID,'outfile = open(outputPth,"w"); \n');
    fprintf(fileID,'outfile.write("J-integral\\nK-I value\\nK-II value\\n"); \n');
    fprintf(fileID,'outfile.close(); \n');
    fprintf(fileID,'counting = 1; \n');
    fprintf(fileID,'for diff_masks in masksVal: \n');
    fprintf(fileID,'	openMdb("%s"); \n','C:\Temp\Step1OMA.cae');
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
    
%% JOB SUBMIT AND RESULT PARSING
    % submit the job 
    fprintf(fileID,'	mdb.jobs["final"].submit(consistencyChecking=OFF) \n');
    % wait for the job to be complete 
    fprintf(fileID,'	mdb.jobs["final"].waitForCompletion(); \n');
    % read output database to retrieve the J-integral values 
    fprintf(fileID,'	time.sleep(3); \n');
    fprintf(fileID,'	odb = session.openOdb("final.odb"); \n');
    fprintf(fileID,'	timestep = odb.steps["Step-1"]; \n');
    fprintf(fileID,'	alloutputs = timestep.historyRegions["ElementSet . ALL ELEMENTS"].historyOutputs; \n');
    fprintf(fileID,'	Jval = []; K1val = []; K2val = []; \n');
    fprintf(fileID,'	for i in alloutputs.keys(): \n');
    fprintf(fileID,'		if i.split()[0] == "J": \n');
    fprintf(fileID,'			Jval.append(alloutputs[i].data[-1][1]); \n');
    fprintf(fileID,'		elif i.split()[0] == "K1": \n');
    fprintf(fileID,'			K1val.append(alloutputs[i].data[-1][1]); \n');
    fprintf(fileID,'		elif i.split()[0] == "K2": \n');
    fprintf(fileID,'			K2val.append(alloutputs[i].data[-1][1]); \n');
    % write values to the file 
    fprintf(fileID,'	outfile = open(outputPth,"a"); \n');	
    fprintf(fileID,'	for i in Jval: \n');
    fprintf(fileID,'		outfile.write(str(i)); \n');
    fprintf(fileID,'		outfile.write("\\t"); \n');
    fprintf(fileID,'	outfile.write("\\n"); \n');	
    fprintf(fileID,'	if(extractK):		 \n');
    fprintf(fileID,'		for i in K1val: \n');
    fprintf(fileID,'			outfile.write(str(i)); \n');
    fprintf(fileID,'			outfile.write("\\t"); \n');
    fprintf(fileID,'		outfile.write("\\n"); \n');		
    fprintf(fileID,'		for i in K2val: \n');
    fprintf(fileID,'			outfile.write(str(i)); \n');
    fprintf(fileID,'			outfile.write("\\t"); \n');
    fprintf(fileID,'		outfile.write("\\n"); \n');	
    fprintf(fileID,'	outfile.write("\\n"); \n');
    fprintf(fileID,'	outfile.close(); \n');	
    fprintf(fileID,'	if counting!=len(masksVal): \n');
    fprintf(fileID,'		odb.close(); \n');
    fprintf(fileID,'		mdb.close() \n');
    % Creates a dummy file for Matlab to know that the analysis is over 
    fprintf(fileID,'finish_file = open("C:\\Temp\\OUROMAdone.tmp","a"); \n');
    fprintf(fileID,'finish_file.write("DONE"); \n');
    fprintf(fileID,'finish_file.close(); \n');

fclose(fileID);
fprintf('\n Printing is Complete .. Open using Notepad ++\n');

% copyfile([pwd '\OUROMA_v3.1.py'],folder)