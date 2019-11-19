function [QualityFact] = FitWilliamsResidual(mesh,estcrktip,mat)

% Read input displacement field
mask = ~isnan(mesh.dDIC(1,:));
Xpos = mesh.UDIC(1,mask); Ypos = mesh.UDIC(2,mask);
dXpos = mesh.dDIC(1,mask); dYpos = mesh.dDIC(2,mask);

% [posX,posY] = meshgrid(unique(Xpos),unique(Ypos));
posx = Xpos(:);
posy = Ypos(:);

% % Set posX and posY
% n = 30;
% lins = linspace(-1,1,n);
% [posX,posY] = meshgrid(lins,lins);
% posx = posX(:);
% posy = posY(:);

%Set constants
E = mat.E;
v = mat.nu;

% Displacement functions
% Material parameters
bulk = (3-v)/(1+v);
G =  E/(2*(1+v));
% Polar co-ordinates
r = @(crcktipx,crcktipy) sqrt((posx+crcktipx).^2+(posy+crcktipy).^2);
th = @(crcktipx,crcktipy,crckang) atan2((posy+crcktipy),(posx+crcktipx))+crckang;
% Analytical displacement fields
uxKI = @(KI,crcktipx,crcktipy,crckang) ...
    KI/(2*G)*sqrt(r(crcktipx,crcktipy)/(2*pi)).*...
    cos(th(crcktipx,crcktipy,crckang)/2).*(bulk-1+2*sin(th(crcktipx,crcktipy,crckang)/2).^2);
uyKI = @(KI,crcktipx,crcktipy,crckang) ...
    KI/(2*G)*sqrt(r(crcktipx,crcktipy)/(2*pi)).*...
    sin(th(crcktipx,crcktipy,crckang)/2).*(bulk+1-2*cos(th(crcktipx,crcktipy,crckang)/2).^2);
% uxKII = @(KII,crcktipx,crcktipy,crckang) ...
%     KII/(2*G)*sqrt(r(crcktipx,crcktipy)/(2*pi)).*...
%     sin(th(crcktipx,crcktipy,crckang)/2).*(bulk+1+2*cos(th(crcktipx,crcktipy,crckang)/2).^2);
% uyKII = @(KII,crcktipx,crcktipy,crckang) ...
%     -KII/(2*G)*sqrt(r(crcktipx,crcktipy)/(2*pi)).*...
%     cos(th(crcktipx,crcktipy,crckang)/2).*(bulk-1-2*sin(th(crcktipx,crcktipy,crckang)/2).^2);
% Rigid body movement + rotation displacement field
uxrig = @(rigx,rot) rigx-tan(rot)*posy;
uyrig = @(rigy,rot) rigy+tan(rot)*posx;


% Analytical displacement + rigid body movement + rotation fields
% ux = @(KI,KII,crcktipx,crcktipy,crckang,urigx,urot) ...
%     uxKI(KI,crcktipx,crcktipy,crckang)+uxKII(KII,crcktipx,crcktipy,crckang)+...
%     uxrig(urigx,urot);
% uy = @(KI,KII,crcktipx,crcktipy,crckang,urigy,urot) ...
%     uyKI(KI,crcktipx,crcktipy,crckang)+uyKII(KII,crcktipx,crcktipy,crckang)+...
%     uyrig(urigy,urot);
ux = @(KI,crcktipx,crcktipy,crckang,urigx,urot) ...
    uxKI(KI,crcktipx,crcktipy,crckang)+uxrig(urigx,urot);
uy = @(KI,crcktipx,crcktipy,crckang,urigy,urot) ...
    uyKI(KI,crcktipx,crcktipy,crckang)+uyrig(urigy,urot);

% Set values to find

% KI = 100;
% KII = 50;
% crcktipx = 0.2;
% crcktipy = -0.1;
% crckang = -pi/8;
% rigx = 1;
% rigy = 1;
% rot = pi/6;
% % Add some noise to data
% percent = 0.02;
% noise = (rand(size(posx))-rand(size(posx)))*range(posx)/n*percent;
% % Analitical field + rigid body movment + rigid body rotation;
% uX0 = ux(KI,KII,crcktipx,crcktipy,crckang,rigx,rot)+noise;
% uY0 = uy(KI,KII,crcktipx,crcktipy,crckang,rigy,rot)+noise;

uX0 = dXpos(:); uY0 = dYpos(:);

%% Optimisation solver
options = optimoptions(@lsqnonlin,'Display','off','MaxFunEvals',1E5,'MaxIter',500);
% algorithm = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

% Rigid body movement and rotation
% xrig = [rigx,rigy,rot];
rig0 = [0,0,0];
urigfun = @(rig) [uxrig(rig(1),rig(3))-uX0;
    uyrig(rig(2),rig(3))-uY0];
Xrig = lsqnonlin(urigfun,rig0,[],[],options);
% printmat([xrig',rig0',Xrig'],'Rigid','rigx rigy rot','known initial calculated');

% Analitical field + rigid body movment + rigid body rotation;
% xall =  [KI,crcktipx,crcktipy,crckang,rigx,rigy,rot];
% Initial guess
xall0 = [1,estcrktip(1),estcrktip(2),0,Xrig(1),Xrig(2),Xrig(3)];
uallfun = @(xall) [ux(xall(1),xall(2),xall(3),xall(4),xall(5),xall(7))-uX0...
    uy(xall(1),xall(2),xall(3),xall(4),xall(6),xall(7))-uY0];
% uallfun = @(xall) [ux(2,estcrktip(1),estcrktip(2),0,xall(5),xall(7))-uX0...
%     uy(2,estcrktip(1),estcrktip(2),0,xall(6),xall(7))-uY0];
Xall = lsqnonlin(uallfun,xall0,[],[],options);
% Show solution
% printmat([xall',xall0',Xall'],'ans','KI crcktipx crcktipy crckang rigx rigy rot','known initial calculated');

FittedResidual = [ux(Xall(1),Xall(2),Xall(3),Xall(4),Xall(5),Xall(7))-uX0 ...
    uy(Xall(1),Xall(2),Xall(3),Xall(4),Xall(6),Xall(7))-uY0];
QualityFact = sum(abs(FittedResidual(:)));

dtocrktip = pdist([estcrktip;Xall(2) Xall(3)],'euclidean');

QualityFact = QualityFact+1E3*dtocrktip;

