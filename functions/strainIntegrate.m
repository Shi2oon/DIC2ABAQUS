% clear
% clc
% 
% 
% load origdisps.mat

function [ReconsUx,ReconsUy] = strainIntegrate(Exx,Eyy,Exy,IniRBM,IniRot)

% [dUxdx,dUxdy] = gradient(xrdUxq);
% [dUydx,dUydy] = gradient(xrdUyq);

% Exx = dUxdx; Eyy = dUydy;
% Exy = dUxdy+dUydx;

ReconsUx =nan(size(Exx)); ReconsUydx = nan(size(Exx));
ReconsUy =nan(size(Exx)); ReconsUxdy = nan(size(Exx));
Ini_x = IniRBM(1);
Ini_y = IniRBM(2);

% Integrate Exx in the x direction
Uxint = zeros(size(Exx));
for i=1:size(Exx,1)
    Uxint(i,:) = intgrad1(Exx(i,:));
end

% Integrate Eyy in the y direction
Uyint = zeros(size(Eyy));
for i=1:size(Eyy,2)
    Uyint(:,i) = intgrad1(Eyy(:,i));
end

ReconsUx(1,:) = Uxint(1,:)+Ini_x;
ReconsUy(:,1) = Uyint(:,1)+Ini_y;

% Second BC at (1,1) dUx/dy connu
% IniD = 0.005128985613646;
IniD = IniRot;
ReconsUx(2,1) = ReconsUx(1,1)+IniD;
ReconsUy(1,2) = ReconsUy(1,1)+Exy(1,1)-IniD;

% So we can reconstruct Ux ligne 2 et Uy colone 2
ReconsUx(2,:) = Uxint(2,:) + ReconsUx(2,1);
ReconsUy(:,2) = Uyint(:,2) + ReconsUy(1,2);

%Comput dUx/dy complete vector
ReconsUydx(:,1) = (ReconsUy(:,2)-ReconsUy(:,1));
ReconsUxdy(1,:) = (ReconsUx(2,:)-ReconsUx(1,:));
ReconsUydx(1,2:end) = Exy(1,2:end)-ReconsUxdy(1,2:end);
ReconsUxdy(2:end,:) = Exy(2:end,:)-ReconsUydx(2:end,:);

% integrate it
intXdata = intgrad1(ReconsUxdy(:,1))+ReconsUx(1,1);
intYdata = intgrad1(ReconsUydx(1,:))+ReconsUy(1,1);
ReconsUx(:,1) = intXdata;
ReconsUy(1,:) = intYdata;

for i=3:size(ReconsUx,1)
    ReconsUx(i,:) = Uxint(i,:) + ReconsUx(i,1); 
end
for i=3:size(ReconsUy,2)
    ReconsUy(:,i) = Uyint(:,i) + ReconsUy(1,i); 
end

%Reconstruct Uy the same way

return



