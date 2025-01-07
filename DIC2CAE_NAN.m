function [el,mesh] = DIC2CAE_NAN(scan_ID)
% Mesh DIC for nodes & elements
%Number of elements in the FE mesh
 kj = importdata(scan_ID);
        alldata = [kj.data(:,1) kj.data(:,2) kj.data(:,3) kj.data(:,4)];  
[ alldata ] = reshapeData( alldata );
xx              = unique(alldata(:,1),'first'); 
mesh.x.length   = length(xx);
mesh.x.xi       = (max(xx))-min(xx)/2;

dispsU_1   = alldata(:,1:2); % x and y
dispsd_1   = alldata(:,3:4); % ux an duy
tmp(1:2,:) = dispsU_1';
tmp(3:4,:) = dispsd_1';
%Guess the dimensions of the image
% mesh.winDIC = [length(unique(alldata(:,1),'first')),length(unique(alldata(:,2),'first'))]
mesh.winDIC = GetSize(tmp(2,:));
% sort the dataset
tmp = sortrows(tmp',[1 2])';
mesh.UDIC   = tmp(1:2,:);
mesh.dDIC   = tmp(3:4,:);
mesh.winFE  = [2 2; mesh.winDIC(1)-1 mesh.winDIC(2)-1];

szelnFE = (mesh.winFE(2,2)-mesh.winFE(1,2))*(mesh.winFE(2,1)-mesh.winFE(1,1));
count = 0;
%For each FE element, find constitutive DIC nodes
elnFE = zeros(szelnFE,4);
for c =  mesh.winFE(1,2):1:mesh.winFE(2,2)-1
    for r =  mesh.winFE(1,1):1:mesh.winFE(2,1)-1
        count = count + 1;
        elnFE(count,:) = ...
            [c+(r)*mesh.winDIC(2), ...
            c+1+(r)*mesh.winDIC(2), ...
            c+1+(r-1)*mesh.winDIC(2), ...
            c+(r-1)*mesh.winDIC(2)];
    end
end

%% Assemble element nodes for FE
idx=1;
for i=1:length(elnFE)
    %check if the element is not zero displacement everywhere
    curDx = mesh.dDIC(1,elnFE(i,:)); % Ux
    curDy = mesh.dDIC(2,elnFE(i,:)); % Uy
    dxy = (sqrt(curDx.^2+curDy.^2)~=0);
    if(dxy)
        curUx = mesh.UDIC(1,elnFE(i,:));
        curUy = mesh.UDIC(2,elnFE(i,:));
        
        el.n(i,:)  = elnFE(i,:);
        el.Ux(i,:) = curUx;
        el.Uy(i,:) = curUy;
        el.dx(i,:) = curDx;
        el.dy(i,:) = curDy;
               idx = idx+1;
    end
end

usf_nod  = unique(el.n(:));                     
usf_nod(usf_nod==0) = [];
mesh.UFE = mesh.UDIC(:,usf_nod);
mesh.dFE = mesh.dDIC(:,usf_nod);

% Remove 0 disp elements
msk   = ~sum(el.n,2)==0;
el.n  = el.n(msk,:);
el.Ux = el.Ux(msk,:);
el.Uy = el.Uy(msk,:);
el.dx = el.dx(msk,:);
el.dy = el.dy(msk,:);

%other info 
mesh.elFE  = size(el.n,1);
mesh.nFE   = length(unique(el.n));

end
function [Size] = GetSize(Aray)

nb=Aray(1);
count = 1;

for i = 2 : length(Aray)
    if(Aray(i) == nb)
        count = count + 1; 
%     else
%         break;
    end
end

 Size(1)=count;
 Size(2)=length(Aray)/count;
end