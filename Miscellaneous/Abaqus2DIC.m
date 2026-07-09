function Abaqus2DIC(filrname)
% filrname = 'C:\Users\ak13\OneDrive - National Physical Laboratory\Papers\DIC2SIF\Rplastic';
fid = fopen([filrname,'.rpt'],'rt') ;
U = textscan(fid,'%s','Delimiter','\n');
U = U{1} ;
idxS = strfind(U, '----------');
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(U, 'Minimum');
idx2 = find(not(cellfun('isempty', idxS)));
Uxy = U(idx1(1)+1:idx2(1)-3) ;% pick  nodes
Uxy = cell2mat(cellfun(@str2num,Uxy,'UniformOutput',false));

fid = fopen([filrname '.inp'],'rt') ;
N = textscan(fid,'%s','Delimiter','\n');
N = N{1} ;
idxS = strfind(N, '*Node');
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(N, '*Element');
idx2 = find(not(cellfun('isempty', idxS)));
Nodes_XY = N(idx1(1)+1:idx2(1)-1) ;% pick  nodes
Nodes_XY = cell2mat(cellfun(@str2num,Nodes_XY,'UniformOutput',false));
fclose all
alldata = [Nodes_XY(:,2:end) Uxy(:,2:end)];

X = linspace(min(alldata(:,1)),max(alldata(:,1)),ceil(sqrt(length(alldata(:,1)))));
Y = linspace(min(alldata(:,2)),max(alldata(:,2)),ceil(sqrt(length(alldata(:,2)))));
[x,y] = meshgrid(X,Y);
F = scatteredInterpolant(alldata(:,1), alldata(:,2), alldata(:,3),'natural');
Ux = F(x, y);
F = scatteredInterpolant(alldata(:,1), alldata(:,2), alldata(:,4),'natural');
Uy = F(x, y);
eyy = diff(Uy);
eyy = [eyy;eyy(1,:)];
eyy = eyy(:,1:88);
eyy(abs(eyy)>8e-5) = NaN;
uy = Uy(:,1:88);
uy(isnan(eyy)) = NaN;
Uy=[uy Uy(:,89:end)];
pcolor(Uy);axis image;set(gca,'Ydir','reverse');shading interp;colormap jet; axis off
Ux(isnan(Uy)) = NaN;    Ux = Ux(:);     Ux(isnan(Ux)) = [];
x(isnan(Uy)) = NaN;     x = x(:);      x(isnan(x)) = [];
y(isnan(Uy)) = NaN;     y = y(:);      y(isnan(y)) = [];
Uy = Uy(:);     Uy(isnan(Uy)) = [];

alldata = table(x(:), y(:), Ux(:), Uy(:),...
    'VariableNames',{'x','y', 'displacement_x', 'displacement_y'});
SaveD = fullfile(filrname,'30Kef_Data.dat');
writetable(alldata, SaveD, 'Delimiter',' ');