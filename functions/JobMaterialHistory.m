function JobMaterialHistory(fileID,MatPunique,MatPtype,len)
    % creating job
    fprintf(fileID,'mdb.Job(atTime=None, contactPrint=OFF, description="", echoPrint=OFF,  \n');
    fprintf(fileID,'	explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,  \n');
    fprintf(fileID,'	memory=90, memoryUnits=PERCENTAGE, model="Model-1", modelPrint=OFF, \n'); 
    fprintf(fileID,'	multiprocessingMode=DEFAULT, name="%s", nodalOutputPrecision=SINGLE,  \n',MatPunique);
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
    fprintf(fileID,'		C14, C24, C34, C44, C15, C25, C35, C45, C55, C16, C26,  \n');
    fprintf(fileID,'		C36, C46, C56, C66), ), type=ANISOTROPIC) \n');
    fprintf(fileID,'mdb.models["Model-1"].HomogeneousSolidSection(material=MaterialName, name="Section-1", thickness=None) \n');
    fprintf(fileID,'mdb.models["Model-1"].parts["sample"].SectionAssignment(offset=0.0, offsetField="", offsetType=MIDDLE_SURFACE, region=Region( \n');
    fprintf(fileID,'	elements=mdb.models["Model-1"].parts["sample"].elements), sectionName="Section-1",  \n');
  	fprintf(fileID,'	thicknessAssignment=FROM_SECTION) \n');
    fprintf(fileID,'mdb.models["Model-1"].rootAssembly.regenerate() \n');
    
  	if  MatPtype == 'A'
        fprintf(fileID,'p = mdb.models["Model-1"].parts["sample"];\n');
        fprintf(fileID,'e = p.elements;\n');
        fprintf(fileID,'elements = e.getSequenceFromMask(mask=("[#ffffffff:%d #1fffff ]", ), );\n',round(len*len/5,0));
        fprintf(fileID,'region = regionToolset.Region(elements=elements);\n');
        fprintf(fileID,'orientation=None;\n');
        fprintf(fileID,'mdb.models["Model-1"].parts["sample"].MaterialOrientation(region=region,\n'); 
        fprintf(fileID,'    orientationType=GLOBAL, axis=AXIS_3, \n');
        fprintf(fileID,'    additionalRotationType=ROTATION_NONE, localCsys=None, fieldName="", \n');
        fprintf(fileID,'    stackDirection=STACK_3)\n');
        
        fprintf(fileID,'mdb.models["Model-1"].parts["sample"].MaterialOrientation( \n');
        fprintf(fileID,'	additionalRotationType=ROTATION_NONE, axis=AXIS_3, fieldName="", localCsys= \n');
        fprintf(fileID,'	None, orientationType=GLOBAL, region=Region( \n');
        fprintf(fileID,'	elements=mdb.models["Model-1"].parts["sample"].elements.getSequenceFromMask( \n');
        fprintf(fileID,'	mask=("[#ffffffff:302 #1fff ]", ), )), stackDirection=STACK_3) \n');
        fprintf(fileID,'mdb.models["Model-1"].rootAssembly.regenerate() \n');
    end
    
    fprintf(fileID,'mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1",  \n');
    fprintf(fileID,'	createStepName="Step-1", name="OutpJint", numberOfContours=nbCtrJint, rebar= \n');
    fprintf(fileID,'	EXCLUDE, sectionPoints=DEFAULT) \n');
    fprintf(fileID,'if(extractK): \n');
    fprintf(fileID,'	mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1",  \n');
    fprintf(fileID,'	createStepName="Step-1", name="OutpKval", numberOfContours=nbCtrJint, rebar= \n');
    fprintf(fileID,'	EXCLUDE, contourType=K_FACTORS, sectionPoints=DEFAULT) \n');
 end