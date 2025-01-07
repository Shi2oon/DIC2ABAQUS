function CreatPartsPrelimMeshing(fileID,stressstat)
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
    if strcmpi(stressstat, 'plane_strain')
        fprintf(fileID,'mdb.models["Model-1"].parts["temp"].setElementType(elemTypes=(ElemType(elemCode=CPE4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, hourglassControl=DEFAULT, distortionControl=DEFAULT), ), regions=(mdb.models["Model-1"].parts["temp"].faces, )) \n');
    else
        fprintf(fileID,'mdb.models["Model-1"].parts["temp"].setElementType(elemTypes=(ElemType(elemCode=CPS4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, hourglassControl=DEFAULT, distortionControl=DEFAULT), ), regions=(mdb.models["Model-1"].parts["temp"].faces, )) \n');
    end
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].generateMesh(); \n');
    % Creating orphan mesh
    fprintf(fileID,'mdb.models["Model-1"].parts["temp"].PartFromMesh(name="sample") \n');
    fprintf(fileID,'del mdb.models["Model-1"].parts["temp"]; \n');
    % Creating calculation step; axis system and assembly instance 
    fprintf(fileID,'mdb.models["Model-1"].StaticStep(name="Step-1", previous="Initial") \n');
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.DatumCsysByDefault(CARTESIAN) \n');
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.Instance(dependent=ON, name="sample-1", part=mdb.models["Model-1"].parts["sample"]) \n');
end