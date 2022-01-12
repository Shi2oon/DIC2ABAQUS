function plotDecomposed_v2(M4)
tmp = sortrows([M4.X(:) M4.Y(:) M4.Z(:) M4.Ux(:) M4.Uy(:)  M4.Uz(:)],[3,1,2]);
[~,dataum ] = reshapeData(tmp);
X1  = squeeze(dataum.X1(:,:,1));    Y1 = squeeze(dataum.Y1(:,:,1));
Uy  = squeeze(dataum.Uy(:,:,1));    Uz = squeeze(dataum.Uz(:,:,1));
Ux = squeeze(dataum.Ux(:,:,1));     Um  = sqrt(Ux.^2+Uy.^2+Uz.^2);
Um3 = 1/2*(Uz-flipud(Uz));
        % in case it is zero as Abaqus won't work
Um4 = 1/2*(Uz+flipud(Uz)) + ones(size(Ux))*1e-12;
        
%%
s1=subplot(2,3,1);  	contourf(X1,Y1,Ux,'LineStyle','none'); 	
title('U_x','fontsize',20);
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(1,:) = c.Limits;     colorbar off; 
s2=subplot(2,3,2);  	contourf(X1,Y1,Uy,'LineStyle','none'); 	
title('U_y','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(2,:) = c.Limits;     colorbar off;
s3=subplot(2,3,3);  	contourf(X1,Y1,Uz,'LineStyle','none'); 	
title('U_z','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(3,:) = c.Limits;     colorbar off;

s5=subplot(2,3,6);  	contourf(X1,Y1,Um3,'LineStyle','none'); 	
title('U^{III}_{Asy.}','fontsize',20);
axis image; axis off; colormap jet; box off; %set(gca,'Ydir','reverse')
c  =colorbar;	cU(5,:) = c.Limits;     colorbar off;
s6=subplot(2,3,5);  	contourf(X1,Y1,Um4,'LineStyle','none'); 	
axis image; axis off; colormap jet; box off; 
c  =colorbar;	cU(6,:) = c.Limits;     colorbar off; 
title('U^{III}_{Sy.}','fontsize',20);
addScale([2 3 6],[X1(:) Y1(:)]);
%
cbax  = axes('visible', 'off');         cU(abs(cU)==1)=0;
caxis(cbax,[min(cU(:)) max(cU(:))]);
h = colorbar(cbax, 'location', 'southoutside','position', [0.1356 0.2432 0.2151 0.0322] );
h.Label.String = [ 'U [\mum]']; 
set([s1 s2 s3 s5 s6],"clim",caxis);
%}
set(gcf,'position',[1 51 1460 950]);  