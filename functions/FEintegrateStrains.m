function [el,mesh] = FEintegrateStrains(el,mesh)

%%  CHEAT TEST VERSION
% disp 'Attention Integration cheat test version'
% load DispsFEbench.mat;
% for i=1:length(mesh.UDIC)
%     indice = findClosest(mesh.UDIC(:,i)',alldata(:,1:2));
%     mesh.dDIC(1,i) = alldata(indice,3); mesh.dDIC(2,i) = alldata(indice,4);
% end
% for i=1:length(mesh.UDIC)
%     pos = find(el.Ux(:)==mesh.UDIC(1,i) & el.Uy(:)==mesh.UDIC(2,i));
%     if(pos)
%         el.dx(pos) = mesh.dDIC(1,i);
%         el.dy(pos) = mesh.dDIC(2,i);
%     end
% end

%%
% Mmat = zeros(length(el.n)*3,length(unique(el.n(:)))*2);
% Bvec = zeros(length(el.n)*3,1);
% Xvec = zeros(length(unique(el.n(:)))*2,1);
elemselection.n = el.n;
Mmat = zeros(length(elemselection.n)*3,length(unique(elemselection.n(:)))*2);
Bvec = zeros(length(elemselection.n)*3,1);
Xvec = zeros(length(unique(elemselection.n(:)))*2,1);

% [usednod,id] = unique(elemselection.n(:));
% allx =el.Ux(id); ally = el.Ux(id);
% load DispsFEbench.mat;
% for i =1:length(usednod)
%     xcur = allx(i); ycur = ally(i);
%     indice = findClosest([xcur ycur],alldata(:,1:2));
%     x0ini(2*i-1:2*i,1) = [alldata(indice,3);alldata(indice,4)];
% end

%element size
lx = el.Ux(2,1) - el.Ux(1,1);
ly = el.Uy(1,2) - el.Uy(1,1);

% Get node correspondance
% used_nodes(1,:) = unique(el.n(:));
% used_nodes(2,:) = 1:length(used_nodes);
used_nodes(1,:) = unique(elemselection.n(:));
used_nodes(2,:) = 1:length(used_nodes);

% Loop through gauss points
cons_nod = nan(1,4);
for i=1:length(elemselection.n)
    idx = (i-1)*3+1:(i-1)*3+3;
    Bvec(idx) = [el.gpexx(i);el.gpeyy(i);el.gpexy(i)];
    %element consitutive nodes
    for j=1:4
        cons_nod(j) = find(used_nodes(1,:)==elemselection.n(i,j));
    end
    %to fit the convention chosen     
    cons_nod = [cons_nod(4) cons_nod(1) cons_nod(2) cons_nod(3)];
    expos = (2.*cons_nod)-1; eypos=expos+1;
    %add Exx related terms
    Mmat(idx(1),expos) =  Mmat(idx(1),expos)+(1/(2*lx)*[-1 1 1 -1]);
    %add Eyy related terms
    Mmat(idx(2),eypos) =  Mmat(idx(2),eypos)+(1/(2*ly)*[-1 -1 1 1]);
    %add Exy related terms
    Mmat(idx(3),expos) =  Mmat(idx(3),expos)+(1/(2*ly)*[-1 -1 1 1]);
    Mmat(idx(3),eypos) =  Mmat(idx(3),eypos)+(1/(2*lx)*[-1 1 1 -1]);
end

% solve the problem
options =  optimoptions('lsqlin','MaxIter',1E4,'Algorithm','trust-region-reflective','PrecondBandWidth',Inf,'TolFun',1E-14);
[Xvec,resnorm,residual,exitflag,output,lambda] = lsqlin(Mmat,Bvec,[],[],[],[],...
    -2*ones(length(unique(elemselection.n(:)))*2,1),2*ones(length(unique(elemselection.n(:)))*2,1),[],options);
dxval = Xvec(1:2:end);dyval = Xvec(2:2:end);


% gather the results
for i=1:length(dxval)
    el.dx(el.n==used_nodes(1,i)) = dxval(i);
    el.dy(el.n==used_nodes(1,i)) = dyval(i);
end

for i=1:length(mesh.UDIC)
    pos = find(el.Ux(:)==mesh.UDIC(1,i) & el.Uy(:)==mesh.UDIC(2,i));
    if(pos)
        mesh.dDIC(1,i) = el.dx(pos(1));
        mesh.dDIC(2,i) = el.dy(pos(1));
    end
end

% Apply smoothing for noise robustness
Matdx = reshape(mesh.dDIC(1,:),mesh.winDIC(2),mesh.winDIC(1));
Matdy = reshape(mesh.dDIC(2,:),mesh.winDIC(2),mesh.winDIC(1));
% smoothn option
% optsmt.Spacing = [lx,ly];
% smthres = smoothn({Matdx,Matdy},'robust',optsmt);
% Matdx = smthres{1}; Matdy = smthres{2};
% mesh.dDIC(1,:) = Matdx(:)'; mesh.dDIC(2,:) = Matdy(:)';

% % Median filtering option
Matdx_f = medfilt2(Matdx,[4,4],'symmetric');
Matdy_f = medfilt2(Matdy,[4,4],'symmetric');
Matdx_f(isnan(Matdx_f)) = Matdx(isnan(Matdx_f));
Matdy_f(isnan(Matdy_f)) = Matdy(isnan(Matdy_f));
mesh.dDIC(1,:) = Matdx_f(:)';
mesh.dDIC(2,:) = Matdy_f(:)';

% Correctly fill 'el' structure
for i=1:length(mesh.UDIC)
    pos = find(el.Ux(:)==mesh.UDIC(1,i) & el.Uy(:)==mesh.UDIC(2,i));
    if(pos)
        el.dx(pos) = mesh.dDIC(1,i);
        el.dy(pos) = mesh.dDIC(2,i);
    end
end

