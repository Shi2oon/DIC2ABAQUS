function printTiffs(fileID,folder,Stat)

fprintf(fileID,'session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues(title=OFF, \n');
    fprintf(fileID,'    state=OFF, annotations=OFF, compass=OFF); \n');
    fprintf(fileID,'session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues( \n');
    fprintf(fileID,'    legendFont="-*-verdana-medium-r-normal-*-*-120-*-*-p-*-*-*"); \n');
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.commonOptions.setValues(  \n');
    fprintf(fileID,'    visibleEdges=NONE);  \n');
    % U1
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  \n');
    fprintf(fileID,'    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U1"), ); \n');
    if Stat == 'N' % for normal view
        fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.display.setValues(plotState=(  \n');
        fprintf(fileID,'    CONTOURS_ON_UNDEF, ));  \n'); % not deform
    elseif Stat == 'D' % for deformed view
        fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.display.setValues(plotState=( \n');
        fprintf(fileID,'    CONTOURS_ON_DEF, ));\n'); % deformed view
    end
    SaveF = pythonFileName(fullfile(folder, [ 'U1' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');    
    % U2
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  \n');
    fprintf(fileID,'    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U2"), ); \n');
    SaveF = pythonFileName(fullfile(folder, [ 'U2' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    
%% save stress   
    % Mises
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    INVARIANT, "Mises"), ); \n');
    SaveF = pythonFileName(fullfile(folder, ['SM' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S11
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S11"), );\n');
    SaveF = pythonFileName(fullfile(folder, [ 'S11' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S12
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S12"), );\n');
    SaveF = pythonFileName(fullfile(folder, [ 'S12' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');
    % S22
    fprintf(fileID,'session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(\n');
    fprintf(fileID,'    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(\n');
    fprintf(fileID,'    COMPONENT, "S22"), );\n');
    SaveF = pythonFileName(fullfile(folder, [ 'S22' Stat]));
    fprintf(fileID,'session.printToFile(fileName="%s", format=TIFF, canvasObjects=(\n',SaveF);
    fprintf(fileID,'    session.viewports["Viewport: 1"], ));\n');