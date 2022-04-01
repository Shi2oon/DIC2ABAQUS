function [Xcrop,Ycrop] = DIC_selectCrop(X,Y,Z,lineX,lineY)
% Plots Z (defined at co-ordinates [X],[Y]) and the associated horizontal
% line defined by lineX and lineY

% Asks for user input to crop a region

Xvec = X(1,:);
Yvec = Y(:,1);
fig = figure;   fig=subplot(1,1,1);
hold on

    imagesc(Xvec,Yvec,Z)
    colormap(jet(256));                 %set(gca,'CLim',[13 15.5]); 

axis image;         axis xy;        fig.YDir='reverse';   
colormap jet;       fig.XDir='reverse'; 
% set(gca,'YDir','normal')
% axis equal
xlim([min(Xvec) max(Xvec)])
ylim([min(Yvec) max(Yvec)])
colorbar
plot(lineX,lineY,'color','k')

hold off
% set(fig,'Name','Crop required data','NumberTitle','off');
pos = get(gcf,'position');          set(gcf,'position',[100 100 pos(3:4)*2]) 
% Xcrop =     [20.9812; 33.5856];  Ycrop =  [-2.4920;13.1597];
% Xcrop =[ 20.0456; 14.3934]; Ycrop = [-1.3322; -5.8397];
[Xcrop,Ycrop] = ginput(2);
Xcrop = [min(Xcrop);max(Xcrop)];
Ycrop = [min(Ycrop);max(Ycrop)];
hold on
plot([Xcrop(1) Xcrop(2) Xcrop(2) Xcrop(1) Xcrop(1)],...
    [Ycrop(1) Ycrop(1) Ycrop(2) Ycrop(2) Ycrop(1)],'color','r')
hold off
% set(fig,'Name','Cropped region','NumberTitle','off');
title('Cropped region shown in red');

%% placeing the crack at the middle
% if abs(mean(lineY)-Ycrop(1)) ~= abs(mean(lineY)-Ycrop(2))
%         addi  = (abs(mean(lineY)-Ycrop(1))+abs(mean(lineY)-Ycrop(2)))/2;
%         Ycrop = [mean(lineY)-addi, mean(lineY)+addi];
% end

hold on
plot([Xcrop(1) Xcrop(2) Xcrop(2) Xcrop(1) Xcrop(1)],...
    [Ycrop(1) Ycrop(1) Ycrop(2) Ycrop(2) Ycrop(1)],'color','k')
hold off
title('Crack Centered-Cropped region shown in black');
