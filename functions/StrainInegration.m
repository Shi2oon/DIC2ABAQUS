function [el,mesh,gl]=StrainInegration(mesh,el,Dir)
close all
    [el,mesh] = FEintegrateStrains(el,mesh);
%Interpolates FE results into a global domain at nodal co-ordinates
    [el,mesh] = RotRemoval('true',mesh,el);
    [gl]      = makeglobal(el, mesh);

%% plotting
close all
h1=subplot(1,2,2); imagesc(gl.Ux(1,:),gl.Uy(:,1),gl.dy); 
title('Integrated U_Y');    colorbar; C = caxis; colormap jet
set(gca,'Ydir','normal');           axis image; 
% h1.XDir='reverse';                  h1.YDir='reverse';
xlabel('X[m]','FontSize',20,'FontName','Times New Roman');          
ylabel('Y[m]','FontSize',20,'FontName','Times New Roman');

h2=subplot(1,2,1); imagesc(gl.Ux(1,:),gl.Uy(:,1),gl.dx); 
title('Integrated U_X');    colorbar; caxis(C);
set(gca,'Ydir','normal');           axis image; 
% h2.XDir='reverse';                  h2.YDir='reverse';
xlabel('X[m]','FontSize',20,'FontName','Times New Roman');          
ylabel('Y[m]','FontSize',20,'FontName','Times New Roman'); 
set(gcf,'position',[30 50 1300 950]);   
       
saveas(gcf,fullfile(Dir.results,[Dir.fillname '_integrated_displacements.png']));       
saveas(gcf,fullfile(Dir.results,[Dir.fillname '_integrated_displacements.fig']));    close all;

if  any(ismember(fields(mesh),'OpS'))==1
    if sum(strfind(mesh.OpS,'xED'))~=0
        CroppedPlot2(Dir.Maps,gl,20)
        pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '_All.png']);
        saveas(gcf,pyxe_D_path); 
        pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '_All.fig']);
        saveas(gcf,pyxe_D_path); close all        
    end
end
end
