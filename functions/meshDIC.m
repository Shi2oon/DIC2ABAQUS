function [el,mesh] = meshDIC(mesh)
% Mesh DIC for nodes & elements
%Number of elements in the FE mesh
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


% function [el,mesh] = meshDIC(mesh)
% 
% if mesh.Operation=='Dis'
% %% Mesh DIC for nodes & elements
% %Number of elements in the FE mesh
% szelnFE = (mesh.winFE(2,2)-mesh.winFE(1,2))*(mesh.winFE(2,1)-mesh.winFE(1,1));
% count = 0;
% %For each FE element, find constitutive DIC nodes
% elnFE = zeros(szelnFE,4);
% for c =  mesh.winFE(1,2):1:mesh.winFE(2,2)-1
%     for r =  mesh.winFE(1,1):1:mesh.winFE(2,1)-1
%         count = count + 1;
%         elnFE(count,:) = ...
%             [c+(r)*mesh.winDIC(2), ...
%             c+1+(r)*mesh.winDIC(2), ...
%             c+1+(r-1)*mesh.winDIC(2), ...
%             c+(r-1)*mesh.winDIC(2)];
%     end
% end
% 
% %% Assemble element nodes for FE
% idx=1;
% for i=1:length(elnFE)
%     %check if the element is not zero displacement everywhere
%     curDx = mesh.dDIC(1,elnFE(i,:)); % Ux
%     curDy = mesh.dDIC(2,elnFE(i,:)); % Uy
%     dxy = (sqrt(curDx.^2+curDy.^2)~=0);
%     if(dxy)
%         curUx = mesh.UDIC(1,elnFE(i,:));
%         curUy = mesh.UDIC(2,elnFE(i,:));
%         
%         el.n(i,:)  = elnFE(i,:);
%         el.Ux(i,:) = curUx;
%         el.Uy(i,:) = curUy;
%         el.dx(i,:) = curDx;
%         el.dy(i,:) = curDy;
%                idx = idx+1;
%     end
% end
% 
% usf_nod  = unique(el.n(:)); usf_nod(usf_nod==0) = [];
% mesh.UFE = mesh.UDIC(:,usf_nod);
% mesh.dFE = mesh.dDIC(:,usf_nod);
% 
% % Remove 0 disp elements
% msk = ~sum(el.n,2)==0;
% el.n = el.n(msk,:);
% el.Ux = el.Ux(msk,:);
% el.Uy = el.Uy(msk,:);
% el.dx = el.dx(msk,:);
% el.dy = el.dy(msk,:);
% 
% %other info 
% mesh.elFE  = size(el.n,1);
% mesh.nFE   = length(unique(el.n));
% 
% elseif mesh.Operation=='Str'
%     
% %% Mesh DIC for nodes & elements
% count = 0;
% %Number of elements in the FE mesh
% szelnFE = (mesh.winFE(2,2)-mesh.winFE(1,2))*(mesh.winFE(2,1)-mesh.winFE(1,1));
% 
% %For each FE element, find constitutive DIC nodes
% elnFE = zeros(szelnFE,4);
% for c =  mesh.winFE(1,2):1:mesh.winFE(2,2)-1
%     for r =  mesh.winFE(1,1):1:mesh.winFE(2,1)-1
%         count = count + 1;
%         elnFE(count,:) = ...
%             [c+(r)*mesh.winDIC(2), ... %
%             c+1+(r)*mesh.winDIC(2), ...
%             c+1+(r-1)*mesh.winDIC(2), ...
%             c+(r-1)*mesh.winDIC(2)];
%     end
% end
% 
% %% Assemble element nodes for FE
% idx=1;
% for i=1:length(elnFE)
%     %check if the element is not zero displacement everywhere
%     curDx = mesh.dDIC(1,elnFE(i,:));
%     curDy = mesh.dDIC(2,elnFE(i,:));
%     dxy = (sqrt(curDx.^2+curDy.^2)~=0);
%     if(dxy)
%         curUx = mesh.UDIC(1,elnFE(i,:));
%         curUy = mesh.UDIC(2,elnFE(i,:));
%         
%         el.n(i,:) = elnFE(i,:);
%         el.Ux(i,:) = curUx;
%         el.Uy(i,:) = curUy;
%         el.dx(i,:) = curDx;
%         el.dy(i,:) = curDy;
%         idx = idx+1;
%     end
% end
% 
% usf_nod = unique(el.n(:));
% mesh.UFE = mesh.UDIC(:,usf_nod);
% mesh.dFE = mesh.dDIC(:,usf_nod);
% 
% %other info 
% mesh.elFE  = size(el.n,1);
% mesh.nFE   = length(unique(el.n));
% end