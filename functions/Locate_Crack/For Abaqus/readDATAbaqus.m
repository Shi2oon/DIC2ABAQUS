function [OutJ,KI, KII, JK] = readDATAbaqus(filrname)
    
fid = fopen(filrname,'rt') ;
S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;

%% J 
idxS = strfind(S, 'OUTPJINT_CRACK-1');
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
idx2 = find(not(cellfun('isempty', idxS)));

% pick  nodes  
J_in = S(idx1(3)+1:idx2(1)-3) ;
count=0;
for iv=1:size(J_in,1)
    J{iv,:} = cell2mat(cellfun(@str2num,J_in(iv),'UniformOutput',false));
    for ij = 1:length(J{iv,:})
        count = count+1;
        OutJ(count)  = abs(J{iv}(ij));
    end     
end
OutJ = OutJ(:);
%% KII and KII
idxS = strfind(S, 'OUTPKVAL_CRACK-1');
idx1 = find(not(cellfun('isempty', idxS)));

K_in = S(idx1(3)+1:idx2(2)-3) ;
for iv = length(K_in):-1:1
    if isempty(K_in{iv,1})
        K_in(iv)=[];
    end
end

idxS = strfind(K_in, 'MTS   DIRECTION (DEG)');
idx1 = find(not(cellfun('isempty', idxS)));
K_in(idx1)=[];
B  = convertCharsToStrings(K_in);
count=0;
for iv = 1:length(B)
    C{iv}  = textscan(B(iv),'%s');
    C{iv}  = C{iv}{1,1};
    format = 0;
    count = count+1;
    for ic = 1:length(C{iv})
        Conv = str2double(C{iv}(ic));
        if isnan(Conv)==0
            format = format+1;
            OutKJ(count,format) = str2double(C{iv}(ic));
        end
    end
end
count=1;
for iv = 1:3:size(OutKJ,1)
    KI(count:count+size(OutKJ,2)-1)  = OutKJ(iv,:);          
    KII(count:count+size(OutKJ,2)-1) = OutKJ(iv+1,:);          
    JK(count:count+size(OutKJ,2)-1)  = OutKJ(iv+2,:);    
    count = count+size(OutKJ,2);
end
KI = KI(:);         KI(KI==0)=[];
KII = KII(:);       KII(KII==0)=[];
JK = JK(:);         JK(JK==0)=[];