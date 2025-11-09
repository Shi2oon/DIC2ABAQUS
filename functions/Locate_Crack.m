function [datum,saf, msk,SaveD] = ...
    Locate_Crack(datum)
% a function to take a dic file and then give the user the librity to
% define the crack tip and mask which will be used in ababqus code
[~,~, saf] = unist4Abaqus(datum,datum.input_unit);% Units +

%% restructure grids
Fx = scatteredInterpolant(datum.X(:),datum.Y(:),datum.Ux(:),'natural');
Fy = scatteredInterpolant(datum.X(:),datum.Y(:),datum.Uy(:),'natural');
xLin  = linspace(min(datum.X(:)),max(datum.X(:)),length(unique(datum.X)));
yLin  = linspace(min(datum.Y(:)),max(datum.Y(:)),length(unique(datum.Y)));
[datum.X,datum.Y] = meshgrid(xLin,yLin);
datum.Ux = Fx(datum.X(:),datum.Y(:));
datum.Uy = Fy(datum.X(:),datum.Y(:));
datum.Ux = reshape(datum.Ux,length(unique(datum.Y)),length(unique(datum.X)));
datum.Uy = reshape(datum.Uy,length(unique(datum.Y)),length(unique(datum.X)));

%%
close all;
% if strcmpi(Operation, 'lazy')
%     imagesc(datum.X(1,:),datum.Y(:,1),datum.E12);
%     c=colorbar; c.Label.String = ['E_{12} [' input_unit ']'];
%     caxis([-5e-3 5e-3])
% else
U  = (datum.Ux.^2+datum.Uy.^2).^0.5;
% [~,~,~,Exy] = xgLOBAL(datum.M4.mesh,datum.M4.E11,datum.M4.E12,datum.M4.X,datum.M4.Y);
imagesc(datum.X(1,:),datum.Y(:,1),U);
c=colorbar; c.Label.String = ['U_{Mag.} [' datum.input_unit ']'];
% end
set(gca,'Ydir','normal');	axis image;
title('Answer in the command line');
xlabel(['X [' datum.input_unit ' ]'],'FontSize',20,'FontName','Times New Roman');
ylabel(['Y [' datum.input_unit ' ]'],'FontSize',20,'FontName','Times New Roman');
set(gcf,'WindowStyle','normal')
set(gcf,'position',[30 50 1300 950]);

%% Crop and rotate
[datum] = Crack_align(datum);   % rotate data
opts.Interpreter = 'tex';       % Include the desired Default answer
opts.Default     = 'N';         % Use the TeX interpreter to format the question
quest            = 'Do you want to crop the map?';
answer           = questdlg(quest,'Boundary Condition','Y','N', opts);
if strcmpi(answer,'Y') % crop data
    [dat] = Cropping10(datum.X,datum.Y,datum.Ux, datum.Uy);
    dat = uniformMesh(dat);
    datum.X = dat.X;    datum.Y = dat.Y; 
    datum.Ux = dat.Ux;  datum.Uy = dat.Uy;
end
%{
 % change mesh
    opts.Interpreter = 'tex';       % Include the desired Default answer
    opts.Default     = 'N';         % Use the TeX interpreter to format the question
    quest            = 'Do you want to change the mesh to irregular mesh?';
    answer           = questdlg(quest,'Boundary Condition','Y','N', opts);
if strcmpi(answer,'Y')
%     [Xq,Yq] = Hyperbolic_Tangent_Clustering([min(datum.X(:)) ...
%         max(datum.X(:))],[min(datum.Y(:)) max(datum.Y(:))],size(datum.X));
      [Xq,Yq] = addRandomMesh(unique(datum.X),unique(datum.Y));
%     plot(Xq,Yq,'k');hold on; plot(Xq.',Yq.','k');hold off; axis equal;
    datum.Ux = interp2(datum.X,datum.Y,datum.Ux,Xq,Yq,'nearest');
    datum.Uy = interp2(datum.X,datum.Y,datum.Uy,Xq,Yq,'nearest');
    datum.X = Xq;      datum.Y = Yq;
    imagesc(Xq(1,:),Yq(:,1),(datum.Ux.^2+datum.Uy.^2).^0.5);
    c=colorbar; c.Label.String = ['U_{Mag.} [' input_unit ']'];
    set(gca,'Ydir','normal');	axis image;
    title('Answer in the command line');
    xlabel(['X [' input_unit ' ]'],'FontSize',20,'FontName','Times New Roman');
    ylabel(['Y [' input_unit ' ]'],'FontSize',20,'FontName','Times New Roman');
    set(gcf,'WindowStyle','normal')
    set(gcf,'position',[30 50 1300 950]);
end
%}
%% get dim
if ~isfield(datum,"msk") || strcmpi(answer,'Y')
    % ATTENTION : crackpoints should be defined with the crack tip in first position.....
    title('Select the crack tip, start from crack tip',FontSize=20);
    [xo,yo] = ginput(2);
    % yo = [yo(1); yo(1)];     %xm = [xo(1); xm(2)]; if the crack is on x axis
    xom = xo; yom = yo;
    [~, index] = min(abs(xLin-xom(2)));      xom(2) = xLin(index);
    line(xom,yom,'Color','w','LineStyle','-.','linewidth',2)
    title('Select the Crack mask, start from crack tip',FontSize=20);
    [xm,ym] = ginput(2);
else
    xo = datum.msk.xo;    yo = datum.msk.yo;
    xm = datum.msk.xm;    ym = datum.msk.ym;
end
%% get excat from data in
xLin       = datum.X(1,:);
% [~, index] = min(abs(xLin-xo(1)));      xo(1) = xLin(index);
[~, index] = min(abs(xLin-xo(2)));      xo(2) = xLin(index);
[~, index] = min(abs(xLin-xm(1)));      xm(1) = xLin(index);
[~, index] = min(abs(xLin-xm(2)));      xm(2) = xLin(index);

yLin       = datum.Y(:,1);
% [~, index] = min(abs(yLin-yo(1)));      yo(1) = yLin(index);
% [~, index] = min(abs(yLin-yo(2)));      yo(2) = yLin(index);
[~, index] = min(abs(yLin-ym(1)));      ym(1) = yLin(index);
msk.ds1 = index;
[~, index] = min(abs(yLin-ym(2)));      ym(2) = yLin(index);
msk.ds2 = index;
%{
%% adjust to min and max
if xo(2)<=min(datum.X(1,:)) || abs(xo(2)-min(datum.X(1,:)))<abs(xo(2)-max(datum.X(1,:)))
    xo(2)=min(datum.X(1,:)); end
if xo(2)>=max(datum.X(1,:)) || abs(xo(2)-min(datum.X(1,:)))>abs(xo(2)-max(datum.X(1,:)))
    xo(2)=max(datum.X(1,:)); end
if xm(2)<=min(datum.X(1,:)) || abs(xm(2)-min(datum.X(1,:)))<abs(xm(2)-max(datum.X(1,:)))
    xm(2)=min(datum.X(1,:)); end
if xm(2)>=max(datum.X(1,:)) || abs(xm(2)-min(datum.X(1,:)))>abs(xm(2)-max(datum.X(1,:)))
    xm(2)=max(datum.X(1,:)); end
%}
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
saveas(gcf, [datum.results '_Corr.tif']); close all

msk.xo = xo;            msk.yo = yo;
msk.xm = xm;            msk.ym = ym;
datum.xo = xo;        datum.yo = yo;
datum.xm = xm;        datum.ym = ym;

%% save cropped data
datum.Ux(isnan(datum.Ux))=0;            datum.Uy(isnan(datum.Uy))=0;
alldata = [datum.X(:), datum.Y(:), datum.Ux(:), datum.Uy(:)];
alldata = table(alldata(:,1), alldata(:,2), alldata(:,3), alldata(:,4),...
    'VariableNames',{'x','y', 'displacement_x', 'displacement_y'});
SaveD = [datum.results '_Data.dat'];
writetable(alldata, SaveD, 'Delimiter',' ');
end