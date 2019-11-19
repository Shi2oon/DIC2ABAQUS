function [Jint] = getq(el, Jint, mesh)

%Assemble (Jint)ergral elements from FE results
 Jint.n     = el.n(Jint.el,:);
 Jint.Ux    = el.Ux(Jint.el,:);
 Jint.Uy    = el.Uy(Jint.el,:);
 Jint.dx    = el.dx(Jint.el,:);
 Jint.dy    = el.dy(Jint.el,:);
 Jint.gpx   = el.gpx(Jint.el,1);
 Jint.gpy   = el.gpy(Jint.el,1);
 Jint.gpexx = el.gpexx(Jint.el,1);
 Jint.gpeyy = el.gpeyy(Jint.el,1);
 Jint.gpexy = el.gpexy(Jint.el,1);
 Jint.gpsxx = el.gpsxx(Jint.el,1);
 Jint.gpsyy = el.gpsyy(Jint.el,1);
 Jint.gpsxy = el.gpsxy(Jint.el,1);
 Jint.gpW   = el.gpW(Jint.el,1);

%get q from surface interpolation using two squares (outer&inner) 
%Z defines q at inner and outer J-intergral path

% disp 'ICI aussi changement'
nodes = reshape(1:length(mesh.UDIC),mesh.winDIC(2),mesh.winDIC(1));
x = mesh.UDIC(1,[nodes(Jint.nout(1,2),Jint.nout(2,1)),nodes(Jint.nin(1,2),Jint.nin(2,1)),...
          nodes(Jint.nin(1,1),Jint.nin(2,2)),nodes(Jint.nout(1,1),Jint.nout(2,2))]);
y = mesh.UDIC(2,[nodes(Jint.nout(1,2),Jint.nout(2,1)),nodes(Jint.nin(1,2),Jint.nin(2,1)),...
          nodes(Jint.nin(1,1),Jint.nin(2,2)),nodes(Jint.nout(1,1),Jint.nout(2,2))]);
      
points = [x(1) y(1) 0;x(1) y(4) 0;x(4) y(4) 0;x(4) y(1) 0;x(2) y(2) 1;...
    x(2) y(3) 1;x(3) y(3) 1;x(3) y(2) 1;];
F = scatteredInterpolant(points(:,1), points(:,2), points(:,3));

Jint.nQ = F(Jint.Ux,Jint.Uy);
Jint.nQ = roundn(Jint.nQ,-5);


% disp 'debug'
% x = [x(1) -9.92 -9.92 x(4)];
% y = [y(1) 0 0 y(4)];
      
% [a,b]=meshgrid(x,y);
% sizeFE = abs(mesh.winFE(1,:)-mesh.winFE(2,:))+1;
% nodes = reshape(1:length(mesh.UFE),sizeFE(2),sizeFE(1));
% 
% % x = mesh.UFE(1,[nodes(Jint.nout(1,1),Jint.nout(2,1)),nodes(Jint.nin(1,1),Jint.nin(2,1)),...
% %           nodes(Jint.nin(1,2),Jint.nin(2,2)),nodes(Jint.nout(1,2),Jint.nout(2,2))]);
% % y = mesh.UFE(2,[nodes(Jint.nout(1,1),Jint.nout(2,1)),nodes(Jint.nin(1,1),Jint.nin(2,1)),...
% %           nodes(Jint.nin(1,2),Jint.nin(2,2)),nodes(Jint.nout(1,2),Jint.nout(2,2))]);
% 
% x = mesh.UFE(1,[nodes(Jint.nout(1,2),Jint.nout(2,1)),nodes(Jint.nin(1,2),Jint.nin(2,1)),...
%           nodes(Jint.nin(1,1),Jint.nin(2,2)),nodes(Jint.nout(1,1),Jint.nout(2,2))]);
% y = mesh.UFE(2,[nodes(Jint.nout(1,2),Jint.nout(2,1)),nodes(Jint.nin(1,2),Jint.nin(2,1)),...
%           nodes(Jint.nin(1,1),Jint.nin(2,2)),nodes(Jint.nout(1,1),Jint.nout(2,2))]);
%       
% [a,b]=meshgrid(x,y);
% Z = [0 0 0 0; 0 1 1 0; 0 1 1 0;0 0 0 0];
% 
% Jint.nQ=griddata(a,b,Z,Jint.Ux,Jint.Uy);