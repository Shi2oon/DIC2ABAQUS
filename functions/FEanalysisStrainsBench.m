function [el] = FEanalysisStrains(el,mat,~)
%%

load AbaqBenchmark;
%initiate global stiffness matrix K
%gl.K = sparse(zeros(mesh.nFE*2));

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

count=1;
%loop through elements
for elnum = 1:size(el.n,1);
    %   get node positions and displacement values
    Ue = zeros(1,size(el.n,2)*2);
    de = zeros(1,size(el.n,2)*2);
    
    %     Ue(1:2:size(el.n,2)*2) = el.Ux(elnum,:);
    %     Ue(2:2:size(el.n,2)*2) = el.Uy(elnum,:);
    %     de(1:2:size(el.n,2)*2) = el.dx(elnum,:);
    %     de(2:2:size(el.n,2)*2) = el.dy(elnum,:);
    
    xel = el.Ux(elnum,:);
    yel = el.Uy(elnum,:);
    dxel = el.dx(elnum,:);
    dyel = el.dy(elnum,:);
    
    Ue(1:2:size(el.n,2)*2) = [xel(4) xel(1) xel(2) xel(3)];
    Ue(2:2:size(el.n,2)*2) = [yel(4) yel(1) yel(2) yel(3)];
    de(1:2:size(el.n,2)*2) = [dxel(4) dxel(1) dxel(2) dxel(3)];
    de(2:2:size(el.n,2)*2) = [dyel(4) dyel(1) dyel(2) dyel(3)];
    
    %get labels of abaqus nodes of this element
%     nodesconst = [];
%     
%     for i=1:4
%         coord = Ue(2*i-1:2*i);
%         posi = find(AbaNodes(:,2)==coord(1)&AbaNodes(:,3)==coord(2));
%         nodesconst(i)=posi(1);
%     end
%     %find element lablel
%     for c = 1:length(AbaElems)
%         if(sum(ismember(AbaElems(c,:),nodesconst))==4)
%             elemindex = c;
%             break;
%         end
%     end
    if(~isnan(AbaLUT(elnum)))
        elemindex=AbaLUT(elnum);
        egp(1,1)   = AbaStrains(elemindex,2);
        egp(2,1)   = AbaStrains(elemindex,3);
        egp(3,1)   = AbaStrains(elemindex,4);
    else
        %get shape functions
        eta  = gaussloc(1);
        xi   = gaussloc(2);
        [~,dHex] = shapefunct(eta, xi, mat);
        %assamble B matrix
        [B,~] = BJmat(dHex, el, Ue);
        %     det(J)
        
        %find strain, stress and position at Gauss point
        egp = B*de';
    end
       
    sgp = mat.D*egp;
    
    %assemble into elements(location of strain, stress at Gauss points)
    el.gpexx(elnum)   = egp(1);
    el.gpeyy(elnum)   = egp(2);
    el.gpexy(elnum)   = egp(3);
    el.gpemax(elnum,1)  = 1/2*(egp(1)+egp(2))+(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2);
    el.gpemin(elnum,1)  = 1/2*(egp(1)+egp(2))-(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2);
    
    el.gpsxx(elnum)   = sgp(1);
    el.gpsyy(elnum)   = sgp(2);
    el.gpsxy(elnum)   = sgp(3);
    
    strpt = roundn([sum(unique(el.Ux(elnum,:)))/2 sum(unique(el.Uy(elnum,:)))/2],-6);
    el.gpx(elnum,:)     = strpt(1);
    el.gpy(elnum,:)     = strpt(2);
    
    if strcmp(mat.stressstate,'plane_strain')
        el.gpszz(elnum) = mat.nu *(el.gpsxx + el.gpsyy);
    end
    %                  el.gpsmax(elnum,gp)  = (sgp(1)+sgp(2))/2+sqrt(((sgp(1)+sgp(2))/2)^2+sgp(3)^2);
    %                  el.gpsmin(elnum,gp)  = (sgp(1)+sgp(2))/2-sqrt(((sgp(1)+sgp(2))/2)^2+sgp(3)^2);
    
    Vec=eig([el.gpsxx(elnum) el.gpsxy(elnum) 0; el.gpsxy(elnum) el.gpsyy(elnum) 0; 0 0 el.gpszz(elnum)]);
    el.gpsmax(elnum) = Vec (1,:);
    el.gpsmed(elnum) = Vec (2,:);
    el.gpsmin(elnum) = Vec (3,:);
    %                  el.gpsmises(elnum,gp)= sqrt(((el.gpsmin(elnum,gp)-el.gpsmax(elnum,gp))^2)+el.gpsmin(elnum,gp)^2+el.gpsmax(elnum,gp)^2)/sqrt(2);
    el.gpsmises(elnum,1)=1/sqrt(2)*sqrt((el.gpsmax(elnum)-el.gpsmed(elnum))^2+(el.gpsmed(elnum)-el.gpsmin(elnum))^2+(el.gpsmin(elnum)-el.gpsmax(elnum))^2);
    %     el.gpx(elnum)     = el.Ux(elnum,:)*H';
    %     el.gpy(elnum)     = el.Uy(elnum,:)*H';
    %     el.gpW(elnum,1)     = 1/2* egp'*mat.D*egp;
    el.gpW(elnum,1)     = 1/2* egp'*sgp;
    el.gpsh(elnum,1)=(sgp(1)+sgp(2))/3;
    %assamble K matrix
    %        Ke = mat.t*B'*mat.D*B*gaussweight(gp)*det(J);
    %        gpindex(1:2:size(el.nFE,2)*2) = (el.nFE(elnum,:)*2)-1;
    %        gpindex(2:2:size(el.nFE,2)*2) = (el.nFE(elnum,:)*2);
    %        gl.K(gpindex, gpindex) = gl.K(gpindex,gpindex) + Ke;
    
end
% el.gpx = mean(el.gpx,2);
% el.gpy = mean(el.gpy,2);
% el.Uxdef = el.Ux;
% el.Uydef = el.Uy ;
el.Uxdef = el.Ux + el.dx*1;
el.Uydef = el.Uy + el.dy*1;
