function [Crop] = Cropping10(X,Y,Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10)
% Plots Z (defined at co-ordinates [X],[Y]) and the associated horizontal
% line defined by lineX and lineY
% Asks for user input to crop a region
close all;                  fig=subplot(1,1,1);
imagesc(X(1,:),Y(:,1),Z1);   title('1^{st} Z data: Select Area to Crop');  
axis image;                 set(gca,'Ydir','normal');   %axis off;  
colorbar;   %colormap jet;                           
set(gcf,'position',[30 50 1300 950])
xlabel('X [Raw Data Units]');          ylabel('Y [Raw Data Units]');

[Xcrop,Ycrop] = ginput(2);
Xcrop = [min(Xcrop);max(Xcrop)];
Ycrop = [min(Ycrop);max(Ycrop)];
hold on
plot([Xcrop(1) Xcrop(2) Xcrop(2) Xcrop(1) Xcrop(1)],...
    [Ycrop(1) Ycrop(1) Ycrop(2) Ycrop(2) Ycrop(1)],'color','k')
hold off

%% Data
xLin          = X(1,:);                     yLin         = Y(:,1); 
[~, Xcrop(1)] = min(abs(xLin-Xcrop(1)));   [~, Xcrop(2)] = min(abs(xLin-Xcrop(2)));         
[~, Ycrop(1)] = min(abs(yLin-Ycrop(1)));   [~, Ycrop(2)] = min(abs(yLin-Ycrop(2)));  

%% Z data
for i=1:nargin-2
    eval(sprintf('Crop.Z%d = Z%d(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));',i,i));
end

%% XY, steps and stifness
Crop.X   = X(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.Y   = Y(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
Crop.X   = Crop.X - min(min(Crop.X));  	Crop.Y   = Crop.Y - min(min(Crop.Y));
if (Crop.X(1) - Crop.X(end))>0;         Crop.X   = flip(Crop.X,2);         end
if (Crop.Y(1) - Crop.Y(end))>0;         Crop.Y   = flip(Crop.Y,1);         end

%%
close all;              fig = subplot(1,1,1); 
imagesc(Crop.X(1,:),Crop.Y(:,1),Crop.Z1)
axis image;             set(gca,'Ydir','normal');   % axis off;  
colormap jet;           colorbar;                            
xlabel('X [Raw Data Units]');          ylabel('Y [Raw Data Units]');
set(gcf,'position',[30 50 1300 950]);    
end