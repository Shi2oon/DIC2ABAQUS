function [X,Y,E11,E22,E12,Nodes,Elements,mesh] = Meshing(DATAin)
% if DIM == 2 && strcmp(DATATYPE,'Strain')
    tmp(1:2,:)   = DATAin(:,1:2)'; % x and y
    tmp(3:size(DATAin,2),:)   = DATAin(:,3:end)'; % data
    %Guess the dimensions of the image
    mesh.winodow = [length(unique(DATAin(:,1),'first')),length(unique(DATAin(:,2),'first'))];
    % sort the dataset
    tmp         = sortrows(tmp',[1 2])';
    mesh.xy     = tmp(1:2,:);
    mesh.Data   = tmp(3:end,:);
    mesh.winFE  = [1 1; mesh.winodow(1) mesh.winodow(2)]; % take 2 elements fromc coners

    %% Mesh for nodes & elements
    %Number of elements in the FE mesh
    szelnFE = (mesh.winFE(2,2)-mesh.winFE(1,2))*(mesh.winFE(2,1)-mesh.winFE(1,1));
    count   = 0;
    %For each FE element, find constitutive DIC nodes
    elnFE = zeros(szelnFE,4);
    for c =  mesh.winFE(1,2):1:mesh.winFE(2,2)-1
        for r =  mesh.winFE(1,1):1:mesh.winFE(2,1)-1
            count = count + 1;
            elnFE(count,:) = [c+(r)*mesh.winodow(2),     c+1+(r)*mesh.winodow(2), ...
                              c+1+(r-1)*mesh.winodow(2), c+(r-1)*mesh.winodow(2)];
        end
    end

    %% Assemble element nodes for FE
    for i=1:length(elnFE)
        el.n(i,:)   = elnFE(i,:);                   % nodes numbers
        el.x(i,:)   = mesh.xy(1,elnFE(i,:));        % nodes corrdinates in X
        el.y(i,:)   = mesh.xy(2,elnFE(i,:));        % nodes coordinates in Y
        el.e11(i,:) = mesh.Data(1,elnFE(i,:));      % nodes coordinates in Data 1
        el.e22(i,:) = mesh.Data(2,elnFE(i,:));      % nodes coordinates in Data 2
        if size(DATAin,2)==5
            el.e12(i,:) = mesh.Data(3,elnFE(i,:));	% nodes coordinates in Data 3
        end
    end

    usf_nod  = unique(el.n(:));             usf_nod(usf_nod==0) = [];
%     mesh.UFE = mesh.xy(:,usf_nod);
%     mesh.dFE = mesh.Data(:,usf_nod);
    % Remove 0 disp elements
    msk      = ~sum(el.n,2)==0;
    el.n     = el.n(msk,:);
    el.x     = el.x(msk,:);             X   = el.x';
    el.y     = el.y(msk,:);             Y   = el.y';
    el.e11   = el.e11(msk,:);           E11 = el.e11';
    el.e22   = el.e22(msk,:);           E22 = el.e22';
    if size(DATAin,2)==5
        el.e12   = el.e12(msk,:);           E12 = el.e12';
    else
        E12 = [];
    end

    %% other info 
%     DATAout    = [mesh.UFE' mesh.dFE'];
%     mesh.elFE  = size(el.n,1);
%     mesh.nFE   = length(unique(el.n));
    NumberEl   = 1:length(DATAin);
    Nodes      = [NumberEl(:) DATAin(:,1:2)];
    NumberEl   = 1:length(el.n);
    Elements   = [NumberEl(:) el.n];
    
%     X   = reshape(DATAout(:,1),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
%     X(:,1) = [];                   X(1,:) = [];
%     Y   = reshape(DATAout(:,2),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
%     Y(:,1) = [];                   Y(1,:) = [];
%     E11 = reshape(DATAout(:,3),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
%     E11(:,1) = [];                 E11(1,:) = [];
%     E22 = reshape(DATAout(:,4),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
%     E22(:,1) = [];                 E22(1,:) = [];
%     E12 = reshape(DATAout(:,5),length(unique(DATAout(:,2))),length(unique(DATAout(:,1))));
%     E12(:,1) = [];                 E12(1,:) = [];
end
