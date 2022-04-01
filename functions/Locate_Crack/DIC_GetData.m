function Maps=DIC_GetData(Maps)
% file name is fname
close all
Maps = DIC_Rot2Crop(Maps);
%{
%% crack coordinates
close all;
% Maps.Ux(Maps.Ux==0)=NaN;    Maps.Uy(Maps.Uy==0)=NaN;    Maps.Uz(Maps.Uz==0)=NaN;
imagesc(Maps.X1(1,:),Maps.Y1(:,1),Maps.Uy);
c=colorbar; c.Label.String = ['U_{y} [mm]'];
%     caxis([-5e-3 5e-3]);
set(gca,'Ydir','normal');	axis image;colormap jet
title('Answer in the command line');
xlabel(['X [mm]'],'FontSize',20);
ylabel(['Y [mm]'],'FontSize',20);
set(gcf,'WindowStyle','normal')
set(gcf,'position',[30 50 1300 950]);
title('Select the Crack, start from crack tip');
[Maps.xo,Maps.yo] = ginput(2);
title('Select the Crack mask, start from crack tip');
[Maps.xm,Maps.ym] = ginput(2);
% Maps.ym=[Maps.yo(1)-Maps.stepsize;Maps.yo(1)+Maps.stepsize];
% Maps.xm=[Maps.xo(1)-Maps.stepsize;Maps.xo(2)+Maps.stepsize];
% [Maps.xm,Maps.ym] = ginput(2);
close
%% Get stiffness tensor
    % the stiffness tensor is defined with respect to the measurement x and y, once x
    % once x and y are rotated the same rotation angle need to be applied to the
    % the stiffness tensor to maintain the Euler angles x and y reference.
    % this can be see clearly when the measurement x an y were changed for (001)
    % Single Si Crystal as the [-110] and [110] pole kept moving.
    % See variation with direction in fig. 3 of https://doi.org/10.1109/JMEMS.2009.2039697
%     Maps.Stiffness  = S2DRot(Maps.Stiffness,Maps.theta); % rotate the stiffness
%}
end