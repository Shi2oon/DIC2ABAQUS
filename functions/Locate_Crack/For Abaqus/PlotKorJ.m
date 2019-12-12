function [J,Keff,KI,KII]=PlotKorJ(saveto,E,offset,byCode)
% fole ris the folder where all the data is vaed
% E is the elasric modulus or the stifness tensor in Pa
% offset: is when abaqus consider your data in meters or standard SI units 
% and your data is actually in mm (offset = 1000) or nm (offest = 1e6) or
% in Km (offset = 1e-3) and so on

%% input values
if size(E) ~=1      % elastic modulus in Pa
E11 = E(1,1);  E12 = E(1,2);  nu = E12/(E11+ E12);  E = E11*(1-2*nu)*(1+nu)/(1-nu);
end

if byCode == 'Y' % if the processing went through without proplems
    dir
    dataum = importdata(fullfile(saveto, 'KJ_Output.txt'));
    J   = dataum.data(1,:)*offset*1e-3;
    KI  = abs(dataum.data(2,:)*sqrt(offset)*1e-6);
    KII = abs(dataum.data(3,:)*sqrt(offset)*1e-6);
    Keff= sqrt(abs(J)*E)*1e-6;
elseif byCode == 'N' % if there was a proplem and done alone
    filrname = dir(fullfile(saveto,'/*.dat'));
    [J,KI, KII, JK] = readDATAbaqus(fullfile(saveto,filrname.name));
    J   = J*offset*1e-3;
    KI  = KI*sqrt(offset)*1e-6;
    KII = KII*sqrt(offset)*1e-6;
    Keff= sqrt(abs(JK*offset)*E)*1e-6;
end
    
%% the Plot
contrs = length(J);                     contrs = contrs - round(contrs/4);
Jtrue  = ((mean(J(contrs:end))+max(J(contrs:end)))/2);
Jdiv   = std(J(contrs:end));
Ktrue  = ((mean(Keff(contrs:end))+max(Keff(contrs:end)))/2);
Kdiv   = std(Keff(contrs:end));

plot(Keff,'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
set(gcf,'position',[600,100,950,650])
xlabel('Contour Number');       ylabel('K (MPa\surdm)')
title ({['Stress Intensity Factor (K_{eff}) = ' num2str(round(Ktrue,2)) ' ± ' ...
    num2str(round(Kdiv,3)) ' MPa\surdm' ];''});
saveas(gcf, [saveto '\K.fig']); saveas(gcf, [saveto '\K.png']); close all

plot(J,'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
set(gcf,'position',[600,100,950,650])
xlabel('Contour Number');       ylabel('J (KJ/m^2)')
title ({['J-integral = ' num2str(round(Jtrue,2)) ' ± ' ...
    num2str(round(Jdiv,4)) ' KJ/m^2' ];''});
saveas(gcf, [saveto '\J.fig']);     saveas(gcf, [saveto '\J.png']); close all

%% for K1 and K2
KtrueI  = ((mean(KI(contrs:end))+max(KI(contrs:end)))/2);
KdivI   = std(KI(contrs:end));
KtrueII = ((mean(KII(contrs:end))+max(KII(contrs:end)))/2);
KdivII  = std(KII(contrs:end));

plot(KI,'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r'); hold on;
plot(KII,'b--o','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b'); hold off
set(gcf,'position',[600,100,950,650])
xlabel('Contour Number');       ylabel('K (MPa m^{1/2})')
title ({['Stress Intensity Factor (K_{I}) = ' num2str(round(KtrueI,2)) ' ± ' ...
    num2str(round(KdivI,3)) ' MPa\surdm' ];...
    ['Stress Intensity Factor (K_{II}) = ' num2str(round(KtrueII,2)) ' ± ' ...
    num2str(round(KdivII,3)) ' MPa\surdm' ]});
legend('K_{I}','K_{II}');   
saveas(gcf, [saveto '\KI and KII.fig']);
saveas(gcf, [saveto '\KI and KII.png']);    close all

%% plot ALL
plot(Keff,'k--o','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k'); hold on;
plot(KI,'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r'); 
plot(KII,'b--o','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b');hold off;

set(gcf,'position',[600,100,950,650])
xlabel('Contour Number');       ylabel('K (MPa m^{1/2})')
title ({['Stress Intensity Factor (K_{eff}) = ' num2str(round(Ktrue,2)) ' ± ' ...
    num2str(round(Kdiv,3)) ' MPa\surdm' ]; ...
    ['Stress Intensity Factor (K_{I}) = ' num2str(round(KtrueI,2)) ' ± ' ...
    num2str(round(KdivI,3)) ' MPa\surdm' ];...
    ['Stress Intensity Factor (K_{II}) = ' num2str(round(KtrueII,2)) ' ± ' ...
    num2str(round(KdivII,3)) ' MPa\surdm' ]});
legend('K_{eff}','K_{I}','K_{II}');   
saveas(gcf, [saveto '\Keff, KI and KII.fig']);
saveas(gcf, [saveto '\Keff, KI and KII.png']); close all

