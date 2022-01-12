function plotDecomposedStrain(uXXd,uYYd,uZZd,uXYd,uXZd,uYZd,Maps)
figure; 
s1=subplot(3,3,1);  	contourf(Maps.X,Maps.Y,Maps.E11,'LineStyle','none');
title([char(949) '_{xx}'],'fontsize',19);
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(3,3,2);  	contourf(Maps.X,Maps.Y,Maps.E12,'LineStyle','none');
title([char(949) '_{xy}'],'fontsize',19);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(3,3,3);  	contourf(Maps.X,Maps.Y,Maps.E31,'LineStyle','none');
title([char(949) '_{xz}'],'fontsize',19);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;
s5=subplot(3,3,5);  	contourf(Maps.X,Maps.Y,Maps.E22,'LineStyle','none');
title([char(949) '_{yy}'],'fontsize',19);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(4,:) = c.Limits;     colorbar off;
s6=subplot(3,3,6);  	contourf(Maps.X,Maps.Y,Maps.E32,'LineStyle','none');
title([char(949) '_{yz}'],'fontsize',19);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;    colorbar off; 
s9=subplot(3,3,9);  	contourf(Maps.X,Maps.Y,Maps.E33,'LineStyle','none');
title([char(949) '_{zz}'],'fontsize',19);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off;
addScale([3 3 9],[Maps.X(:) Maps.Y(:)]);

EId   = sqrt(0.5*((uXXd(:,:,1)-uYYd(:,:,1)).^2+(uYYd(:,:,1)-uZZd(:,:,1)).^2+...
                  (uZZd(:,:,1)-uXXd(:,:,1)).^2+ ...
                  (uXYd(:,:,1).^2+uYZd(:,:,1).^2+uXZd(:,:,1).^2).*6));
EIId  = sqrt(0.5*((uXXd(:,:,2)-uYYd(:,:,2)).^2+(uYYd(:,:,2)-uZZd(:,:,2)).^2+...
                  (uZZd(:,:,2)-uXXd(:,:,2)).^2+ ...
                  (uXYd(:,:,2).^2+uYZd(:,:,2).^2+uXZd(:,:,2).^2).*6));
EIIId = sqrt(0.5*((uXXd(:,:,3)-uYYd(:,:,3)).^2+(uYYd(:,:,3)-uZZd(:,:,3)).^2+...
                  (uZZd(:,:,3)-uXXd(:,:,3)).^2+ ...
                  (uXYd(:,:,3).^2+uYZd(:,:,3).^2+uXZd(:,:,3).^2).*6));
 
s4=subplot(3,3,4);  	contourf(Maps.X,Maps.Y,EId,'LineStyle','none'); 	
title([char(949) '^{I}_M'],'fontsize',19);
axis image; axis off;  box off; colormap jet;
c  =colorbar;	cU(7,:) = c.Limits;     colorbar off; 
s7=subplot(3,3,7);  	contourf(Maps.X,Maps.Y,EIId,'LineStyle','none'); 	
title([char(949) '^{II}_M'],'fontsize',19);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(8,:) = c.Limits;     colorbar off;
s8=subplot(3,3,8);  	contourf(Maps.X,Maps.Y,EIIId,'LineStyle','none'); 	
title([char(949) '^{III}_M'],'fontsize',19);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(9,:) = c.Limits;    colorbar off; 
%
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'westoutside','position', [0.9011 0.1211 0.0121 0.7533] );
h.Label.String = [char(949)]; 
h.Label.FontSize = 30; 
set([s1 s2 s3 s5 s6 s7 s8 s9],"clim",caxis); 
%}
set(gcf,'position',[1 -41 1900 1000]); 
end