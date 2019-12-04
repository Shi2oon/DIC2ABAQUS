function Maps = readAbaqusOutput(Dirresults,Dirunique)
%% iinp file
filrname = [Dirresults '\' Dirunique '.inp'];
fid = fopen(filrname,'rt') ;
S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;
%%Get the line number of mises 
idxS = strfind(S, 'Node');
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(S, 'Element');
idx2 = find(not(cellfun('isempty', idxS)));
idxS = strfind(S, 'End');
idx3 = find(not(cellfun('isempty', idxS)));
idxS = strfind(S, 'Elset');
idx4 = find(not(cellfun('isempty', idxS)));
if isempty(idx4)==0; idx3(1) = idx4; end

% pick  nodes  
nodes = S(idx1+1:idx2-1) ;
nodes = cell2mat(cellfun(@str2num,nodes,'UniformOutput',false));
% nodes = dlmread(fname,',',[idx1,0,idx2-2,2]); % another way to do it

% pick elements 
% elements = S(idx2+1:idx3(1)-1) ;
% count = 0 ;
% ele = cell(length(elements)/2,1) ;
% for i = 1:2:length(elements)
%     count = count+1 ;
%     ele{count,1} = [elements{i} elements{i+1}];
% end
% ele = cell2mat(cellfun(@str2num,ele,'UniformOutput',false));

%% rpt file rpt fileGet stress 
filrname = [Dirresults '\' Dirunique '.rpt'];
fid = fopen(filrname,'rt') ;
S = textscan(fid,'%s','Delimiter','\n');
S = S{1} ;
% Mises (element    Node    Stress)
idxS = strfind(S, 'S.Mises');
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(S, 'Minimum');
idx2 = find(not(cellfun('isempty', idxS)));
Data = S(idx1(1)+3:idx2(1)-3) ;
Mises = cell2mat(cellfun(@str2num,Data,'UniformOutput',false));
% S11
idxS = strfind(S, 'S.S11');
idx1 = find(not(cellfun('isempty', idxS)));
Data = S(idx1(1)+3:idx2(2)-3) ;
S11  = cell2mat(cellfun(@str2num,Data,'UniformOutput',false));
% S22
idxS = strfind(S, 'S.S22');
idx1 = find(not(cellfun('isempty', idxS)));
Data = S(idx1(1)+3:idx2(3)-3);
S22  = cell2mat(cellfun(@str2num,Data,'UniformOutput',false));
% S12
idxS = strfind(S, 'S.S12');
idx1 = find(not(cellfun('isempty', idxS)));
Data = S(idx1(1)+3:idx2(4)-3);
S12  = cell2mat(cellfun(@str2num,Data,'UniformOutput',false));
fclose('all');

%% organise data 
alldata = zeros(length(nodes),6);
for iv = 1:length(nodes)
    [fk,~]=ind2sub(size(Mises),find(Mises(:,2)==nodes(iv,1)));
    alldata(iv,:) =[nodes(iv,2) nodes(iv,3)  S11(fk(1),end) S22(fk(1),end) ...
                    S12(fk(1),end)  Mises(fk(1),end)];
end

[ Maps ] = reshapeAbqusData( alldata );
