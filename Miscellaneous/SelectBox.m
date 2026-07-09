function ebsdBox = SelectBox(ebsd)
%% GRADIENT
clc;        close all;      warning off
plot(ebsd,ebsd.Eequi); caxis([-1e-1 1e-1]); colormap jet
% [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);
% plot(grains.boundary); 

%% define a sub region
counterss=1;        isDispOk = 'Y'; 
while isDispOk ~= 'N'
    uiwait(msgbox('Click two points on the box corners.','Select Corners','modal'));  
    [xdata{counterss}, ydata{counterss}] = ginput(2);     hold on;
    xmin = min(xdata{counterss});              xmax = max(xdata{counterss});
    ymin = min(ydata{counterss});              ymax = max(ydata{counterss});

    region = [xmin ymin xmax-xmin ymax-ymin];   % marke the sub region
    h = rectangle('position',region,'edgecolor','r','LineStyle','--','linewidth',2);
    
    opts.Interpreter = 'tex'; % Include the desired Default answer
    opts.Default     = 'Y';     % Use the TeX interpreter to format the question
    quest            = '(Y) Chose Another Bix, (C) Remove Previous Selection, (N) Done with Box Selection';
    isDispOk         = questdlg(quest,'Boundary Condition','Y','N','C', opts);    
    
    if isDispOk == 'Y' 
        % select EBSD data within region
        condition           = inpolygon(ebsd,region);  % select indices by polygon
        ebsdBox{counterss}  = ebsd(condition);  
        text(mean(xdata{counterss}),mean(ydata{counterss}),...
            num2str(counterss),'Color','red','FontSize',14)
        counterss           = counterss+1;
        figure; 
    elseif isDispOk == 'C' 
        delete(h)
    end 
end