function [OutJ] = JustJ(filrname,oh)
    
fid = fopen(filrname,'rt') ;
S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;

%% J 
try
    idxS = strfind(S, 'OUTPJINT_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
    idx2 = find(not(cellfun('isempty', idxS)));
    J_in = S(idx1(3)+1:idx2(1)-3) ;% pick  nodes  
catch
    idxS = strfind(S, 'CRACK-1_');
    idx1 = find(not(cellfun('isempty', idxS)));
    idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
    idx2 = find(not(cellfun('isempty', idxS)));
    J_in = S(idx1(3)+1:idx2(1)-3) ;% pick  nodes  
end

count=0;
for iv=1:size(J_in,1)
    J{iv,:} = cell2mat(cellfun(@str2num,J_in(iv),'UniformOutput',false));
    for ij = 1:length(J{iv,:})
        count = count+1;
        OutJ(count)  = abs(J{iv}(ij));
    end     
end
OutJ = OutJ(2:end)';
fclose all;
if ~exist("oh",'var')
close all; plot(OutJ)
set(gcf,'position',[98 311 1481 667])
text(1:length(OutJ),OutJ,string([1:length(OutJ)]))
pause(0.1)
oh = input('where to cut the contour? ');     
end
OutJ=OutJ(1:oh);
close