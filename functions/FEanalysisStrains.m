function [el] = FEanalysisStrains(el,mat,mesh)
%% for strain
%Set Guass quad
[gaussloc,~] = gaussquad(mat);

%initiate strain, stress and gauss location points
el.gpexx = zeros(length(el.n),1);
el.gpeyy = zeros(length(el.n),1);
el.gpexy = zeros(length(el.n),1);
el.gpsxx = zeros(length(el.n),1);
el.gpsyy = zeros(length(el.n),1);
el.gpsxy = zeros(length(el.n),1);
el.gpx   = zeros(length(el.n),1);
el.gpy   = zeros(length(el.n),1);
el.gpszz = zeros(length(el.n),1);
el.gpsmax= zeros(length(el.n),1);
el.gpsmed= zeros(length(el.n),1);
el.gpsmin= zeros(length(el.n),1);

% disp 'Edit Special Non Linear Felix. Attention';
%loop through elements
for elnum = 1:size(el.n,1)
%   get node positions and displacement values
    Ue = zeros(1,size(el.n,2)*2);
%     de = zeros(1,size(el.n,2)*2);
    
    xel = el.Ux(elnum,:);
    yel = el.Uy(elnum,:);
%     dxel = el.dx(elnum,:);
%     dyel = el.dy(elnum,:);
    
    Ue(1:2:size(el.n,2)*2) = [xel(4) xel(1) xel(2) xel(3)];
    Ue(2:2:size(el.n,2)*2) = [yel(4) yel(1) yel(2) yel(3)];   
%     de(1:2:size(el.n,2)*2) = [dxel(4) dxel(1) dxel(2) dxel(3)];
%     de(2:2:size(el.n,2)*2) = [dyel(4) dyel(1) dyel(2) dyel(3)];
    
    
    %get shape functions
    eta  = gaussloc(1);
    xi   = gaussloc(2);
    [~,dHex] = shapefunct(eta, xi, mat);
    %assemble B matrix
    [B,~] = BJmat(dHex, el, Ue);
    
    %find strain, stress and position at Gauss point
%     strpt = roundn([sum(unique(el.Ux(elnum,:)))/2 sum(unique(el.Uy(elnum,:)))/2],-6);
   
    strpt = [sum(unique(el.Ux(elnum,:)))/2 sum(unique(el.Uy(elnum,:)))/2];
    el.gpx(elnum,:)     = strpt(1);
    el.gpy(elnum,:)     = strpt(2);
    
    indice = findClosest(strpt,mesh.straindata(:,1:2));
    egp = mesh.straindata(indice,3:5)';
%     egp = B*de';
     %% inject HR-EBSD data   
    if  any(ismember(fields(mesh),'OpS'))==1
        if sum(strfind(mesh.OpS,'xED'))~=0 
        indice = findClosest(strpt,mesh.stressdata(:,1:2));
        sgp    = mesh.stressdata(indice,3:5)'; 
        end
    else
        sgp = mat.D*egp;
    end
%     sgp = sign(egp).*(5.082*1e9*abs(egp).^3-3.092*1e7*abs(egp).^2+74949*abs(egp));
    
    %assemble into elements(location of strain, stress at Gauss points)
    el.gpexx(elnum)     = egp(1);
    el.gpeyy(elnum)     = egp(2);
    el.gpexy(elnum)     = egp(3);
    el.gpemax(elnum,1)  = 1/2*(egp(1)+egp(2))+(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2);
    el.gpemin(elnum,1)  = 1/2*(egp(1)+egp(2))-(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2);    
    el.gpsxx(elnum)     = sgp(1);
    el.gpsyy(elnum)     = sgp(2);
    el.gpsxy(elnum)     = sgp(3);
    
    if strcmp(mat.stressstate,'plane_strain')
        el.gpszz(elnum) = mat.nu *(el.gpsxx(elnum) + el.gpsyy(elnum));
    end
 
    Vec=eig([el.gpsxx(elnum) el.gpsxy(elnum) 0; el.gpsxy(elnum) el.gpsyy(elnum) 0;...
        0 0 el.gpszz(elnum)]);
    el.gpsmax(elnum) = Vec (1,:);
    el.gpsmed(elnum) = Vec (2,:);
    el.gpsmin(elnum) = Vec (3,:);
    el.gpsmises(elnum,1)=1/sqrt(2)*sqrt((el.gpsmax(elnum)-el.gpsmed(elnum))^2+...
        (el.gpsmed(elnum)-el.gpsmin(elnum))^2+(el.gpsmin(elnum)-el.gpsmax(elnum))^2);
    el.gpW(elnum,1)     = 1/2* egp'*sgp;
    el.gpsh(elnum,1)=(sgp(1)+sgp(2))/3;    
end

el.Uxdef = el.Ux;% + el.dx*1;
el.Uydef = el.Uy;% + el.dy*1;
