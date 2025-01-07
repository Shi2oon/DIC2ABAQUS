
function [OutJ,KI,KII,JK,Der_Deg] = readDATAbaqus(filrname)
fid = fopen(filrname,'rt') ;
S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;
fclose all;
%% DJ_DOWN
% try
%     idxS = strfind(S, 'OUTPJINT_CRACK-1');
%     idx1 = find(not(cellfun('isempty', idxS)));
%     idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
%     idx2 = find(not(cellfun('isempty', idxS)));
%     J_in = S(idx1(3)+1:idx2(1)-3) ;% pick  nodes  
% catch
%     idxS = strfind(S, 'CRACK-1_');
%     idx1 = find(not(cellfun('isempty', idxS)));
%     idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
%     idx2 = find(not(cellfun('isempty', idxS)));
%     J_in = S(idx1(3)+1:idx2(2)-3) ;
% end
% 
% count=0;
% for iv=1:size(J_in,1)
%     J{iv,:} = cell2mat(cellfun(@str2num,J_in(iv),'UniformOutput',false));
%     for ij = 1:length(J{iv,:})
%         count = count+1;
%         OutJ(count)  = abs(J{iv}(ij));%
%     end     
% end
% OutJ = OutJ(2:end)';

try
    idxS = strfind(S, 'OUTPJINT_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
    idx2 = find(not(cellfun('isempty', idxS)));
    J_in = S(idx1(3)+1:idx2(1)-3) ;% pick  nodes  
catch
    idxS = strfind(S, 'H-OUTPUT-2_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    idxS = strfind(S, 'LABELS REFERENCED IN THE ABOVE TABLE');
    idx2 = find(not(cellfun('isempty', idxS)));
    J_in = S(idx1(3)+1:idx2(1)-3) ;
end
J_in = convertCharsToStrings(J_in);
countX=0; j=0;
for iv=1:size(J_in,1)
%     Jvalue{iv,:} = cell2mat(cellfun(@str2num,J_in(iv),'UniformOutput',false));
    if isempty(J_in{iv})
        break;
    end
    Jvalue = textscan(J_in{iv},'%s');
    for ij = 1:length(Jvalue{1})
        if isnan(str2double(Jvalue{1}(ij)))
            j = j+1;countX=0;
        else
            countX = countX+1;
            OutJ(j,countX) = str2double(Jvalue{1}(ij));
        end
    end   
    clear Jvalue
end
%% KII and KII
try
    idxS = strfind(S, 'OUTPKVAL_CRACK-1');
    idx1 = find(not(cellfun('isempty', idxS)));
    K_in = S(idx1(3)+1:idx2(2)-3) ;
catch
    idxS = strfind(S, 'CRACK-2_');
    idx1 = find(not(cellfun('isempty', idxS)));
    K_in = S(idx1(3)+1:idx2(2)-3) ;
end

for iv = length(K_in):-1:1
    if isempty(K_in{iv,1})
        K_in(iv)=[];
    end
end

idxS = strfind(K_in, 'MTS   DIRECTION (DEG)');
idx1 = find(not(cellfun('isempty', idxS)));
if isempty(idx1)
    idxS = strfind(K_in, 'MERR  DIRECTION (DEG)');
    idx1 = find(not(cellfun('isempty', idxS)));
end
if ~isempty(idx1)
    DIRECTION = convertCharsToStrings(K_in(idx1));
    count=0;
    for iv=1:size(DIRECTION,1)
        O{iv}  = textscan(DIRECTION(iv),'%s');
        O{iv}  = O{iv}{1,1};
        format = 0;
        count = count+1;
        for ic = 1:length(O{iv})
            Conv = str2double(O{iv}(ic));
            if isnan(Conv)==0
                format = format+1;
                Der_Deg(count,format) = str2double(O{iv}(ic));
            end
        end   
    end
    Der_Deg = Der_Deg'; 
    Der_Deg = Der_Deg(:);
else
    Der_Deg = [];
end

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
fclose('all');
%{
%% acessive countors
KI(KI==0)  =[];            KII(KII==0)=[];                     JK(JK==0)  =[];  
KLM = OutJ(end);    [A,~] = ismember(OutJ,KLM);	  OutJ(A)= [];   OutJ(end+1) = KLM;
KLM = KI(end);      [A,~] = ismember(KI,KLM);	  KI(A)  = [];   KI(end+1)   = KLM;
KLM = KII(end);     [A,~] = ismember(KII,KLM);	  KII(A) = [];   KII(end+1)  = KLM;
KLM = JK(end);      [A,~] = ismember(JK,KLM);	  JK(A)  = [];   JK(end+1)   = KLM;

Kon = min([length(OutJ) length(JK) length(KI) length(KII)]);
OutJ(Kon:end)=[];   KI(Kon:end)=[]; KII(Kon:end)=[];    JK(Kon:end)=[];
%}
close all; plot(OutJ); hold on; plot(JK); legend('J','J_K')%trim acess 
set(gcf,'position',[98 311 1481 667])
text(1:length(JK),JK,string([1:length(JK)]))
pause(0.1)
oh = input('where to cut the contour? ');               
OutJ=OutJ(1:oh);KI=KI(1:oh);KII=KII(1:oh);JK=JK(1:oh);
close
if ~isempty(Der_Deg)
    Der_Deg=Der_Deg(1:oh);
end
end
