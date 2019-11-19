function CroppedPlot2(Maps,gl,fontsize)
set(0,'defaultAxesFontSize',fontsize);    set(0,'DefaultLineMarkerSize',fontsize/4*3)  
%% Arrange data and limits
        S11  = Maps.S11;              S12  = Maps.S12;
        S13  = Maps.S13;              S22  = Maps.S22;
        S23  = gl.SMISES;	
        unit = '[GPa]';                     LAB  = '\sigma';

Wo   = Maps.Wo;
GNDs = real(log10(Maps.GND));  
% Maxiall = max([E11(:); E12(:); E13(:); E22(:); E23(:); E33(:)]);
% Miniall = min([E11(:); E12(:); E13(:); E22(:); E23(:); E33(:)]);
Xvec = Maps.X(1,:)*1e3;       
Yvec = Maps.Y(:,1)*1e3;

%% Remove extreme values from plot
factor = 10;
S11(abs(S11)>factor*nanmean(abs(S11(:)))) = NaN;
S12(abs(S12)>factor*nanmean(abs(S12(:)))) = NaN;
S13(abs(S13)>factor*nanmean(abs(S13(:)))) = NaN;
S22(abs(S22)>factor*nanmean(abs(S22(:)))) = NaN;
S23(abs(S23)>factor*nanmean(abs(S23(:)))) = NaN;
Wo(abs(Wo)>factor*nanmean(abs(Wo(:)))) = NaN;
GNDs(abs(GNDs)>factor*nanmean(abs(GNDs(:)))) = NaN;

%% Ploting
Maxiall = 1.5;        Miniall=-1.5;
close all
figh = figure(1);               colormap jet
h1 = subplot(3,3,1);            imagesc(Xvec,Yvec,S11);         hold on
contour(Maps.X*1e3,Maps.Y*1e3,S11,10,'LineWidth',1,'LineColor',[0 0 0]);
set(gca,'Ydir','normal');       axis equal;     axis tight;     axis off
title([ LAB '_1_1']);           hold off                                  

h2 = subplot(3,3,2);            imagesc(Xvec,Yvec,S12);         hold on
contour(Maps.X*1e3,Maps.Y*1e3,S12,10,'LineWidth',1,'LineColor',[0 0 0]);
set(gca,'Ydir','normal');       axis equal;     axis tight;     axis off
title([ LAB '_1_2']);           hold off                                

h3 = subplot(3,3,3);            imagesc(Xvec,Yvec,S13);         hold on
contour(Maps.X*1e3,Maps.Y*1e3,S13,10,'LineWidth',1,'LineColor',[0 0 0]);
set(gca,'Ydir','normal');       axis equal;     axis tight;            
c = colorbar;                   c.Label.String = unit;      %labelling  
xlabel('x[\mum]');          	ylabel('y[\mum]');
title([ LAB '_1_3']);           hold off

h5 = subplot(3,3,5);            imagesc(Xvec,Yvec,S22);         hold on
contour(Maps.X*1e3,Maps.Y*1e3,S22,10,'LineWidth',1,'LineColor',[0 0 0]);
set(gca,'Ydir','normal');       axis equal;     axis tight;     axis off           
title([ LAB '_2_2']);           hold off

h6 = subplot(3,3,6);            imagesc(gl.Ux(1,:)*1e6,gl.Uy(:,1)*1e6,S23);   hold on
contour(gl.Ux*1e6,gl.Uy*1e6,S23,10,'LineWidth',1,'LineColor',[0 0 0]);
set(gca,'Ydir','normal');       axis equal;     axis tight;     axis off                     
title([ LAB '_V_M']);           hold off

h4 = subplot(3,3,4);            imagesc(Xvec,Yvec,Wo);          hold on
contour(Maps.X*1e3,Maps.Y*1e3,Wo,10,'LineWidth',1,'LineColor',[0 0 0]);
set(gca,'Ydir','normal');       axis equal;     axis tight;     axis off
caxis([0 1.5]);                 title('W');     hold off
% colorbar;                       

h7 = subplot(3,3,7);            imagesc(Xvec,Yvec,GNDs);        hold on
set(gca,'Ydir','normal');       axis equal;     axis tight;     axis off
set(gca,'CLim',[14 15.5]);      c = colorbar;                  	
title('\rho_G_N_D_s');          c.Label.String = 'log10(m/m^{3})';%labelling       
hold off

h8 = subplot(3,3,8);            
pcolor(gl.Ux*1e6,gl.Uy*1e6,gl.dy*1e6);        shading interp;         hold on
contour(gl.Ux*1e6,gl.Uy*1e6,gl.dy*1e6,10,'LineWidth',1,'LineColor',[0 0 0]);
colormap(jet);                  title('U_y');
c = colorbar;                   c.Label.String = '[\mum]';%labelling
axis equal;                     axis tight; axis off
Lim = c.Limits;                 hold off                 

h9 = subplot(3,3,9);            
pcolor(gl.Ux*1e6,gl.Uy*1e6,gl.dx*1e6);        shading interp;         hold on
contour(gl.Ux*1e6,gl.Uy*1e6,gl.dx*1e6,10,'LineWidth',1,'LineColor',[0 0 0]);
colormap(jet);                  title('U_x');
c = colorbar;                   c.Label.String = '[\mum]';%labelling
caxis(Lim);                     hold off
axis equal;                     axis tight; axis off

% pos = get(figh,'position');
% set(figh,'position',[pos(1:2)/4 pos(3:4)*2])
set(figh,'position',[30 50 1100 1550])

set([h1 h2 h3 h5 h6 h4],'clim',[Miniall Maxiall]);
end

