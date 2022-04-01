function crop = DIC_Rot2Crop(Maps)
Maps.U = {      Maps.Ux         zeros(size(Maps.Ux))	zeros(size(Maps.Ux));...
          zeros(size(Maps.Ux))         Maps.Uy        zeros(size(Maps.Ux));...
          zeros(size(Maps.Ux)) zeros(size(Maps.Ux))   zeros(size(Maps.Ux))};
%
%% DEFINE HORIZONTAL LINE
Maps.X = Maps.X1;     Maps.Y = Maps.Y1;
H1 = figure; H1=subplot(1,1,1);
Xvec = Maps.X(1,:);     Yvec = Maps.Y(:,1);
imagesc(Xvec,Yvec,Maps.U{2,2}); colormap jet; 
axis image;         axis xy;        H1.YDir='reverse';   
colormap jet;       H1.XDir='reverse'; 
axis equal; title('pick the crack start from the tip')
xlim([min(Xvec) max(Xvec)]);    ylim([min(Yvec) max(Yvec)])
pos = get(gcf,'position'); set(gcf,'position',[100 100 pos(3:4)*2]) 
[ptX,ptY] = ginput(2);hold on
plot(ptX,ptY);  hold off; close all

%%
theta = atan((ptY(2)-ptY(1))/(ptX(2)-ptX(1)));
% if atand((ptY(2)-ptY(1))/(ptX(2)-ptX(1)))<5
%     %}
%     theta = 0;
% end
X = Maps.X(:);              Y = Maps.Y(:);
a = DirectionCosine(theta);      % edited by Abdo 22.10.19
crop.theta = rad2deg(theta);

%% Rotate the locations of the data points
% Define a new set of points, with respect to the original axes, but whose
% locations are rotated by an angle theta about the origin.
Xnew   = a(1,1).*X + a(1,2)*Y;            Ynew = a(2,1).*X + a(2,2)*Y; 
Xnew   = double(Xnew);                    Ynew = double(Ynew);
ptXnew = a(1,1).*ptX + a(1,2)*ptY;      ptYnew = a(2,1).*ptX + a(2,2)*ptY;
Maps.stepsize = unique(diff(X));    crop.stepsize = Maps.stepsize(2);
% create a coarse grid of points to visualise data
i = (min(Xnew):crop.stepsize:max(Xnew));
j = (min(Ynew):crop.stepsize:max(Ynew));
% xq and yq are the co-ordinates of a grid parellel to the crack,
% defined with respect to the original frame of reference.
[xq,yq] = meshgrid(i,j);
xq      = double(xq);                       yq = double(yq);
% interpolate data onto this new grid of xq,yq to allow a contour plot to
% be made
[fi,fj] = size(Maps.U);
for i = 1:fi
    for j = 1:fj
        Maps.rot.U{i,j}     = griddata(Xnew,Ynew,Maps.U{i,j}(:),xq,yq);
    end
end

%% Perform tensor transformation to new, rotated, co-ordinate system
% R = @(theta)[cosd(theta) sind(theta) 0;-sind(theta) cosd(theta) 0;0 0 1];
% A0 = R(theta)*A0*R(theta)';
% Displacement Tensor:
[Maps.txf.U]    = componentTensorTransform(Maps.rot.U,a);

%% CROP THE DATA
[selXcrop,selYcrop] = DIC_selectCrop(xq,yq,Maps.txf.U{2,2},ptXnew,ptYnew);
% saveas(gcf,[Maps.results  '\Cropped Zone.tif']); 
% saveas(gcf,[Maps.results  '\Cropped Zone.fig']);close all
% Generate crop Mask
cropMask                   = xq;
cropMask(xq<min(selXcrop)) = 0;
cropMask(xq>max(selXcrop)) = 0;
cropMask(yq<min(selYcrop)) = 0;
cropMask(yq>max(selYcrop)) = 0;
cropMask(cropMask~=0)      = 1;
 
[Ycrop, Xcrop] = find(cropMask~=0);
Xcrop          = [min(Xcrop) max(Xcrop)];
Ycrop          = [min(Ycrop) max(Ycrop)];

% Crop the data
% Crop the x-coordinates
Maps.crop.X = xq(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
% Crop the y-coordinates
Maps.crop.Y = yq(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop)); 

for i = 1:3
    for j = 1:3
        Maps.crop.U{i,j}    = Crop(Maps.txf.U{i,j},Xcrop,Ycrop);
    end
end
crop.Ux = Maps.crop.U{1,1};     
crop.Uy = Maps.crop.U{2,2};

% Shift the X and Y axes to give values starting at zero
crop.X1 = Maps.crop.X - min(Maps.crop.X(:));
crop.Y1 = Maps.crop.Y - min(Maps.crop.Y(:));
end