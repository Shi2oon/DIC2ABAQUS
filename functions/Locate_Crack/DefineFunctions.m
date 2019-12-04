function DefineFunctions(fileID)
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
end