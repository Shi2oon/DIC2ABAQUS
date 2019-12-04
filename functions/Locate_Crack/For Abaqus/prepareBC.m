function prepareBC(fileID,SaveModel)
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
    fprintf(fileID,'mdb.saveAs("%s",); \n',SaveModel);
    fprintf(fileID,'mdb.close() \n');
end