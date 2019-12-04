function [alldata,DataFile,saf] = Locate_Crack(datum,DataFile,Operation,mechDat)
% a function to take a dic file and then give the user the librity to
% define the crack tip and mask which will be used in ababqus code

% suggsted imporvement: data location for mask and y values
if Operation == 'DIC'
    try
        filed = [DataFile '.dat'];
        data  = importdata(filed);           datum = data.data;
    catch
        filed = [DataFile '.mat'];
        datum = load(filed);               datum = datum.alldata;
    end
    switch mechDat.input_unit
        case 'mm'
            % convert mm --> m
            datum = datum.* (mechDat.pixel_size * 10^-3); 
        case 'um'
            % convert um --> m
            datum = datum.* (mechDat.pixel_size * 10^-6); % convert pixels --> m
        case 'm'
            datum = datum.* (mechDat.pixel_size); 
        otherwise
            error('Invalid input_unit')
    end
    mechDat.input_unit = 'm';
    DataFile = fileparts(DataFile);
end
%% Units + 
[datum,mechDat.input_unit, saf] = unist4Abaqus(datum,mechDat.input_unit);

%%   
datum = reshapeData(datum);
close all;      
if      mechDat.type  == 'A';       fig=subplot(1,1,1); 
% imagesc(datum.X(1,:),datum.Y(:,1),real(log10(mechDat.Maps.GND)));   set(gca,'CLim',[14 15.5]); 
imagesc(datum.X(1,:),datum.Y(:,1),mechDat.Maps.S11);                set(gca,'CLim',[-1.5 1.5]); 
fig.XDir='reverse';             fig.YDir='reverse';                 c = colorbar; 	
        c.Label.String = '\rho_G_N_D_s [log10(m/m^{3}])';
else;   imagesc(datum.X(1,:),datum.Y(:,1),datum.Uy); 
        c=colorbar; c.Label.String = ['U_Y [' mechDat.input_unit ']'];     
end  
set(gca,'Ydir','normal');	axis image;
title('Answer in the command line');
xlabel(['X [' mechDat.input_unit ' ]'],'FontSize',20,'FontName','Times New Roman');          
ylabel(['Y [' mechDat.input_unit ' ]'],'FontSize',20,'FontName','Times New Roman');

set(gcf,'position',[30 50 1300 950]); 

%% Crop and rotate
[datum] = Crack_align(datum); % rotate data
    opts.Interpreter = 'tex'; % Include the desired Default answer
    opts.Default     = 'N';     % Use the TeX interpreter to format the question
    quest            = 'Do you want to crop the map?';
    answer           = questdlg(quest,'Boundary Condition','Y','N', opts);
if answer == 'Y' || answer == 'y'
    [datum] = Cropping10(datum.X,datum.Y,datum.Uy, datum.Ux);
    saveas(gcf,[DataFile '\Cropped Uy map.fig']);     
    saveas(gcf,[DataFile '\Cropped Uy map_.png']);    
end

%% get dim
% ATTENTION : crackpoints should be defined with the crack tip in first position.....
title('U_Y :: Select the crack tip start from crack tip');   
[xo,yo] = ginput(2);
title('U_Y :: Select the Crack mask, start from crack tip');
[xm,ym] = ginput(2);

%% get excat from data in
xLin       = datum.X(1,:);
[~, index] = min(abs(xLin-xo(1)));      xo(1) = xLin(index);
[~, index] = min(abs(xLin-xo(2)));      xo(2) = xLin(index);
[~, index] = min(abs(xLin-xm(1)));      xm(1) = xLin(index);
[~, index] = min(abs(xLin-xm(2)));      xm(2) = xLin(index);

yLin       = datum.Y(:,1);
[~, index] = min(abs(yLin-yo(1)));      yo(1) = yLin(index);
[~, index] = min(abs(yLin-yo(2)));      yo(2) = yLin(index);
[~, index] = min(abs(yLin-ym(1)));      ym(1) = yLin(index);    msk.ds1 = index;
[~, index] = min(abs(yLin-ym(2)));      ym(2) = yLin(index);    msk.ds2 = index;

%% adjust to min and max
if xo(2)<=min(datum.X(1,:)) || abs(xo(2)-min(datum.X(1,:)))<abs(xo(2)-max(datum.X(1,:)))
    xo(2)=min(datum.X(1,:)); end
if xo(2)>=max(datum.X(1,:)) || abs(xo(2)-min(datum.X(1,:)))>abs(xo(2)-max(datum.X(1,:)))
    xo(2)=max(datum.X(1,:)); end
if xm(2)<=min(datum.X(1,:)) || abs(xm(2)-min(datum.X(1,:)))<abs(xm(2)-max(datum.X(1,:)))
    xm(2)=min(datum.X(1,:)); end
if xm(2)>=max(datum.X(1,:)) || abs(xm(2)-min(datum.X(1,:)))>abs(xm(2)-max(datum.X(1,:))) 
    xm(2)=max(datum.X(1,:)); end
hold on; plot(xm,ym,'.w');
rectangle('Position',[min(xo(1),xo(2)),min(yo(1),yo(2)),abs(xo(2)-xo(1)),...
    abs(yo(2)-yo(1))],'EdgeColor','W'); drawnow;
rectangle('Position',[min(xm(1),xm(2)),min(ym(1),ym(2)),abs(xm(2)-xm(1)),...
    abs(ym(2)-ym(1))],'EdgeColor','w','FaceColor','k');drawnow;
line(xo,yo,'Color','w','LineStyle','-.','linewidth',2)

title({'Crack Position (x,y) from tip to the end is at ',[ num2str(xo(1)) ' , ' ...
       num2str(yo(1)) ' and ' num2str(xo(2)) ' , ' num2str(yo(2))], ...
       ['Marker position ' num2str(xm(2)) ' , ' num2str(ym(2)) ' and ' ...
       num2str(xm(1)) ' , ' num2str(ym(1))], ['Picture Frame is at ' ...
       num2str(min(datum.X(1,:))) ' , ' num2str(min(datum.Y(:,1))) ' and ' ...
       num2str(max(datum.X(1,:))) ' , ' num2str(max(datum.Y(:,1)))],''});
% fprintf(['Crack Position (x,y) from tip to the end is at ',[ num2str(xo(1)) ' , ' ...
%        num2str(yo(1)) ' and ' num2str(xo(2)) ' , ' num2str(yo(2))], ...
%        ['Marker position ' num2str(xm(2)) ' , ' num2str(ym(2)) ' and ' ...
%        num2str(xm(1)) ' , ' num2str(ym(1))], ['Picture Frame is at ' ...
%        num2str(min(datum.X(1,:))) ' , ' num2str(min(datum.Y(:,1))) ' and ' ...
%        num2str(max(datum.X(1,:))) ' , ' num2str(max(datum.Y(:,1)))],'']);
saveas(gcf, [DataFile '\DIC2ABAQUS Coodrinate.png']); 
msk.xo = xo;      msk.yo = yo;      msk.xm = xm;      msk.ym = ym;

%% save cropped data
if answer == 'Y' || answer == 'y'
alldata = table(datum.X(:), datum.Y(:), datum.Z1(:), datum.Z2(:), 'VariableNames',{'x',...
        'y', 'displacement_x', 'displacement_y'});
    SaveD = [DataFile '\Cropped Data Uxy.dat'];
writetable(alldata, SaveD, 'Delimiter',' '); clc
else
alldata = table(datum.X(:), datum.Y(:), datum.Ux(:), datum.Uy(:), 'VariableNames',{'x',...
        'y', 'displacement_x', 'displacement_y'});
    SaveD = [DataFile '\Data Uxy.dat'];
writetable(alldata, SaveD, 'Delimiter',' '); clc
end

%% writing up
DataFile = PrintCode(mechDat, msk,SaveD,min(size(datum.X)));
% Out = importdata('P:\Abdo\ABAQUS\test.txt'); % output file 