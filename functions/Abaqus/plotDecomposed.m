function plotDecomposed(M4)
tmp = sortrows([M4.X(:) M4.Y(:) M4.Z(:) M4.Ux(:) M4.Uy(:)  M4.Uz(:)],[3,1,2]);
[~,dataum ] = reshapeData(tmp);
X1  = squeeze(dataum.X1(:,:,1));    Y1 = squeeze(dataum.Y1(:,:,1));
Uy  = squeeze(dataum.Uy(:,:,1));    Uz = squeeze(dataum.Uz(:,:,1));
Ux = squeeze(dataum.Ux(:,:,1));     Um  = sqrt(Ux.^2+Uy.^2+Uz.^2);
Um1 = ((1/2*(Ux+flipud(Ux))).^2 + (1/2*(Uy-flipud(Uy))).^2).^0.5;
Um2 = ((1/2*(Ux-flipud(Ux))).^2+  (1/2*(Uy+flipud(Uy))).^2).^0.5;
Um3 = 1/2*(Uz-flipud(Uz));
        
%%
s1=subplot(2,4,1);  	contourf(X1,Y1,Ux,'LineStyle','none'); 	
title('U_x','fontsize',20);
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(2,4,2);  	contourf(X1,Y1,Uy,'LineStyle','none'); 	
title('U_y','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(2,4,3);  	contourf(X1,Y1,Uz,'LineStyle','none'); 	
title('U_z','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;

s5=subplot(2,4,5);  	contourf(X1,Y1,Um1,'LineStyle','none'); 	
title('U^{I}_m','fontsize',20);
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off; 
s6=subplot(2,4,6);  	contourf(X1,Y1,Um2,'LineStyle','none'); 	
title('U^{II}_m','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;     colorbar off;
s7=subplot(2,4,7);  	contourf(X1,Y1,Um3,'LineStyle','none'); 	
title('U^{III}_m','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;
s8=subplot(1,4,4);contourf(X1,Y1,Um,'LineStyle','none')
axis image; axis off; colormap jet;  box off; %set(gca,'Ydir','reverse')
c  =colorbar;	     colorbar off;%cU(4,:) = c.Limits;
addScale([1 4 4],[X1(:) Y1(:)]);title('U_m','fontsize',20);
%
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'southoutside','position', [0.3513 0.0833 0.3 0.03] );
h.Label.String = [ 'U [\mum]']; 
set([s1 s2 s3 s5 s6 s7 s8],"clim",caxis);
%}
set(gcf,'position',[1 51 1900 1000]);  