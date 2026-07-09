function [Maps] = BoxHrEBSD(xdata,ydata,fname,DirSave,ebsd)
close all
load(fname,'Maps');

% Dir.xEBSD  =  [erase(fname,'.ctf') '_XEBSD.mat'];         load(Dir.xEBSD);

for iv = 1:2
    if     iv == 1 % the whole mape
        % do nothing
        xmin          = min(xdata);                        xmax = max(xdata);
        ymin          = min(ydata);                        ymax = max(ydata);
        region        = [xmin ymin xmax-xmin ymax-ymin];   % marke the sub region
    elseif iv==2   % the selcted bit
        %% find location in the data set
        [~, index(1)] = min(abs(Maps.X-xdata(1)));      xdata(1) = Maps.X(index(1));
        [~, index(2)] = min(abs(Maps.X-xdata(2)));     	xdata(2) = Maps.X(index(2));
        [~, indey(1)] = min(abs(Maps.Y-ydata(1)));     	ydata(1) = Maps.Y(indey(1));
        [~, indey(2)] = min(abs(Maps.Y-ydata(2)));     	ydata(2) = Maps.Y(indey(2));
        
        %% draw the box
        condition = inpolygon(ebsd,region);  % select indices by polygon
        ebsd      = ebsd(condition);
        [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',5*degree); 
        ebsd(grains (grains.grainSize<=5))  = [];
        grains    = grains ('indexed');       
        grains.boundary = grains.boundary('indexed');
        grains    = smooth(grains ,5);
        
        %% get data on the box
        ctx         = 0;
        for i       = min(index):1:max(index)
            ctx     = ctx+1;                    cty = 0;    
            for ii  = min(indey):1:max(indey)
                cty = cty+1;                     SqData.S33(cty,ctx) = Maps.S33(ii,i);
        SqData.GND(cty,ctx) = GND(ii,i);         SqData.S11(cty,ctx) = Maps.S11(ii,i);
        SqData.S12(cty,ctx) = Maps.S12(ii,i);    SqData.S13(cty,ctx) = Maps.S13(ii,i);
        SqData.S22(cty,ctx) = Maps.S22(ii,i);    SqData.S23(cty,ctx) = Maps.S23(ii,i);
            end
        end
        % reverse to enter the loop
        Maps.GND = SqData.GND;    Maps.S11 = SqData.S11;	Maps.S12 = SqData.S12;   
        Maps.S13 = SqData.S13;    Maps.S22 = SqData.S22;    Maps.S23 = SqData.S23;                
        Maps.S33 = SqData.S33;    xplot    = Maps.X(min(index):max(index));  
        Maps.X   = xplot;         yplot    = Maps.Y(min(indey):max(indey));  
        Maps.Y   = yplot;
    end
    
    %% plot desired area
    newMtexFigure('nrows',3,'ncols',3);         %counT  = 1; 
    for i = 1:3
        for ii = 1:3
            if ii >= i
                nextAxis(i,ii)
                %subplot(3,3,counT)
                plot(grains.boundary,'micronBar','off');      hold on
                eval(sprintf('imagesc(flip(Maps.X),Maps.Y,Maps.S%d%d);',i,ii));
                if i==3 && ii==3
                    plot(grains.boundary);
                    if iv ==1
                        rectangle('position',region,'edgecolor','r',...
                                    'LineStyle','--','linewidth',2);
                        text(min(xdata),mean(ydata),'Study Area',...
                                    'Color','red','FontSize',14)
                    end
                else
                    plot(grains.boundary,'micronBar','off'); 
                end
                colormap(jet(256));     set(gca,'CLim',[-1.5 1.5]);
                hold off;               title(['\sigma_' num2str(i) '_' num2str(ii)]);
                mtexColorbar;           mtexColorbar;               %axis off;
            end                        %;counT = counT+1;
        end
    end
    setColorRange('equal');                                 CLim(gcm,'equal'); 
    mtexColorbar('all','title',['\sigma (GPa), Step Size = ' num2str(Maps.stepsize) ' \mum']);
    
    % GNDs
    nextAxis(3,1)
        plot(grains.boundary,'micronBar','off');      hold on
        imagesc(flip(x),y,GND); 
        plot(grains.boundary,'micronBar','off');      hold off
        colormap(jet(256));                 
        set(gca,'ColorScale','log');         % this works only starting with Matlab 2018a
        set(gca,'CLim',[10^13 10^15.5]);     c = colorbar;    
        title('log(\rho_G_N_D_s), m/m^{3}'); % labelling

    %% save
    if     iv == 1     
        saveas(gcf,fullfile(DirSave,'All Area.png'));  close all
    elseif iv == 2
        DSave = fullfile(DirSave,'Box.png');     
        saveas(gcf,DSave);                       close all
    end
end
