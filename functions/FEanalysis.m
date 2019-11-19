function [el] = FEanalysis(el,mat)
%% 
%initiate global stiffness matrix K
 %gl.K = sparse(zeros(mesh.nFE*2));

%Set Guass quad
 [gaussloc gaussweight] = gaussquad(mat);

%initiate strain, stress and gauss location points
 el.gpexx = zeros(size(el.n)); 
 el.gpeyy = zeros(size(el.n)); 
 el.gpexy = zeros(size(el.n));
 el.gpsxx = zeros(size(el.n)); 
 el.gpsyy = zeros(size(el.n)); 
 el.gpsxy = zeros(size(el.n));
 el.gpx   = zeros(size(el.n)); 
 el.gpy   = zeros(size(el.n));
 el.gpszz = zeros(size(el.n));
 el.gpsmax= zeros(size(el.n));
 el.gpsmed= zeros(size(el.n));
 el.gpsmin= zeros(size(el.n));

%loop through elements
 for elnum = 1:size(el.n,1);
     %Ue = [el.UxFE(1,:);el.UyFE(1,:)];
     Ue = zeros(1,size(el.n,2)*2);
     de = zeros(1,size(el.n,2)*2);
     Ue(1:2:size(el.n,2)*2) = el.Ux(elnum,:);
     Ue(2:2:size(el.n,2)*2) = el.Uy(elnum,:);
     de(1:2:size(el.n,2)*2) = el.dx(elnum,:);
     de(2:2:size(el.n,2)*2) = el.dy(elnum,:);
    
    %loop through Gauss points
     switch mat.eltype
        case 'Q4'
            for gp = 1:4
                %get shape functions
                 eta  = gaussloc(1,gp);
                 xi   = gaussloc(2,gp);
                 [H dHex] = shapefunct(eta, xi, mat);
                %assamble B matrix
                 [B J] = BJmat(dHex, el, Ue);
                %find strain, stress and position at Gauss point
                 egp = B*de';
                 sgp = mat.D*B*de';
                %assemble into elements(location of strain, stress at Gauss points)        
                 el.gpexx(elnum,gp)   = egp(1); 
                 el.gpeyy(elnum,gp)   = egp(2); 
                 el.gpexy(elnum,gp)   = egp(3);
                 el.gpemax(elnum,gp)  = 1/2*(egp(1)+egp(2))+(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2);
                 el.gpemin(elnum,gp)  = 1/2*(egp(1)+egp(2))-(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2); 
%                  el.gpemax(elnum,gp)  = (egp(1)+egp(2))/2+sqrt(((egp(1)+egp(2))/2)^2+egp(3)^2); %CORRECTION
%                  el.gpemin(elnum,gp)  = (egp(1)+egp(2))/2-sqrt(((egp(1)+egp(2))/2)^2+egp(3)^2); %CORRECTION
                 
                 el.gpsxx(elnum,gp)   = sgp(1); 
                 el.gpsyy(elnum,gp)   = sgp(2); 
                 el.gpsxy(elnum,gp)   = sgp(3);
                 
                 if strcmp(mat.stressstate,'plane_strain')
                     el.gpszz(elnum,gp) = mat.nu *(el.gpsxx + el.gpsyy);
                 end
%                  el.gpsmax(elnum,gp)  = (sgp(1)+sgp(2))/2+sqrt(((sgp(1)+sgp(2))/2)^2+sgp(3)^2);
%                  el.gpsmin(elnum,gp)  = (sgp(1)+sgp(2))/2-sqrt(((sgp(1)+sgp(2))/2)^2+sgp(3)^2);
                 
                 Vec=eig([el.gpsxx(elnum,gp) el.gpsxy(elnum,gp) 0; el.gpsxy(elnum,gp) el.gpsyy(elnum,gp) 0; 0 0 el.gpszz(elnum,gp)]);
                 el.gpsmax(elnum,gp) = Vec (1,:);
                 el.gpsmed(elnum,gp) = Vec (2,:);
                 el.gpsmin(elnum,gp) = Vec (3,:);
%                  el.gpsmises(elnum,gp)= sqrt(((el.gpsmin(elnum,gp)-el.gpsmax(elnum,gp))^2)+el.gpsmin(elnum,gp)^2+el.gpsmax(elnum,gp)^2)/sqrt(2);
                 el.gpsmises(elnum,gp)=1/sqrt(2)*sqrt((el.gpsmax(elnum,gp)-el.gpsmed(elnum,gp))^2+(el.gpsmed(elnum,gp)-el.gpsmin(elnum,gp))^2+(el.gpsmin(elnum,gp)-el.gpsmax(elnum,gp))^2);
                 el.gpx(elnum,gp)     = el.Ux(elnum,:)*H'; 
                 el.gpy(elnum,gp)     = el.Uy(elnum,:)*H';
                 el.gpW(elnum,gp)     = 1/2* egp'*mat.D*egp;
                 el.gpsh(elnum,gp)=(sgp(1)+sgp(2))/3;
                 %assamble K matrix
        %        Ke = mat.t*B'*mat.D*B*gaussweight(gp)*det(J); 
        %        gpindex(1:2:size(el.nFE,2)*2) = (el.nFE(elnum,:)*2)-1;
        %        gpindex(2:2:size(el.nFE,2)*2) = (el.nFE(elnum,:)*2);
        %        gl.K(gpindex, gpindex) = gl.K(gpindex,gpindex) + Ke;
            end
        case {'Q8' ,'Q9'}
            if strcmp(mat.eltype,'Q8')
                maxgp = 8;
            else
                maxgp = 9;
            end
            for gp = 1:maxgp
                %get shape functions
                 eta  = gaussloc(1,gp);
                 xi   = gaussloc(2,gp);
                 [H dHex] = shapefunct(eta, xi, mat);
                %assamble B matrix
                 [B J] = BJmat(dHex, el, Ue);
                %find strain, stress and position at Gauss point
                 egp = B*de';
                 sgp = mat.D*B*de';
                %assemble into elements(location of strain, stress at Gauss points)        
                 el.gpx(elnum,gp)     = el.Ux(elnum,:)*H'; 
                 el.gpy(elnum,gp)     = el.Uy(elnum,:)*H';
                 el.gpexx(elnum,gp)   = egp(1); 
                 el.gpeyy(elnum,gp)   = egp(2); 
                 el.gpexy(elnum,gp)   = egp(3);
                 el.gpemax(elnum,gp)  = 1/2*(egp(1)+egp(2))+(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2);
                 el.gpemin(elnum,gp)  = 1/2*(egp(1)+egp(2))-(1/2*(egp(1)+egp(2))^2+egp(3)^2)^(1/2);
                 el.gpsxx(elnum,gp)   = sgp(1); 
                 el.gpsyy(elnum,gp)   = sgp(2); 
                 el.gpsxy(elnum,gp)   = sgp(3);
%                  el.gpsmax(elnum,gp)  = 1/2*(sgp(1)+sgp(2))+(1/2*(sgp(1)+sgp(2))^2+sgp(3)^2)^(1/2);
%                  el.gpsmin(elnum,gp)  = 1/2*(sgp(1)+sgp(2))-(1/2*(sgp(1)+sgp(2))^2+sgp(3)^2)^(1/2);
%                  el.gpsmises(elnum,gp)= (el.gpsmin(elnum,gp)^2+el.gpsmax(elnum,gp)^2)-el.gpsmin(elnum,gp)*el.gpsmax(elnum,gp);
%                  el.gpsmises(elnum,gp)= sqrt(((el.gpsmin(elnum,gp)-el.gpsmax(elnum,gp))^2)+el.gpsmin(elnum,gp)^2+el.gpsmax(elnum,gp)^2)/sqrt(2);

                 Vec=eig([el.gpsxx(elnum,gp) el.gpsxy(elnum,gp) 0; el.gpsxy(elnum,gp) el.gpsyy(elnum,gp) 0; 0 0 el.gpszz(elnum,gp)]);
                 el.gpsmax(elnum,gp) = Vec (1,:);
                 el.gpsmed(elnum,gp) = Vec (2,:);
                 el.gpsmin(elnum,gp) = Vec (3,:);
                 el.gpsmises(elnum,gp)=1/sqrt(2)*sqrt((el.gpsmax(elnum,gp)-el.gpsmed(elnum,gp))^2+(el.gpsmed(elnum,gp)-el.gpsmin(elnum,gp))^2+(el.gpsmin(elnum,gp)-el.gpsmax(elnum,gp))^2);
                 
                 el.gpW(elnum,gp)     = 1/2* egp'*mat.D*egp;
                 el.gpsh(elnum,gp)=(sgp(1)+sgp(2))/3;
                 %assamble K matrix
        %        Ke = mat.t*B'*mat.D*B*gaussweight(gp)*det(J); 
        %        gpindex(1:2:size(el.nFE,2)*2) = (el.nFE(elnum,:)*2)-1;
        %        gpindex(2:2:size(el.nFE,2)*2) = (el.nFE(elnum,:)*2);
        %        gl.K(gpindex, gpindex) = gl.K(gpindex,gpindex) + Ke;
            end
     end
end
 
 el.Uxdef = el.Ux + el.dx*1;
 el.Uydef = el.Uy + el.dy*1; 
 