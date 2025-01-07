function ReadingInputFile(fileID)
    fprintf(fileID,'x=[];   y=[]; \n');
    fprintf(fileID,'vx=[];  vy=[]; \n');
    fprintf(fileID,'x0=[];  y0=[]; \n');
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
%     fprintf(fileID,'# nonodesx=y.count(y[100]); nonodesy=x.count(x[1700]; \n');
    fprintf(fileID,'nonodesx=len(set(x));				nonodesy=len(set(y)); \n');
    % Round nodes position to get rid off machine error 
    fprintf(fileID,'x = [ round(elem, 8) for elem in x ] \n');
    fprintf(fileID,'y = [ round(elem, 8) for elem in y ] \n');
    % Determine mesh size in each direction 
    fprintf(fileID,'unitsizex=abs((max(x)-min(x))/(nonodesx-1)) \n');
    fprintf(fileID,'unitsizey=abs((max(y)-min(y))/(nonodesy-1)) \n');
    % Determine the rounding order to apply later; it is define as 5%% of the unitsize 
    fprintf(fileID,'rndparam = -int(floor((log10(0.05*unitsizex)))); \n');
 end 