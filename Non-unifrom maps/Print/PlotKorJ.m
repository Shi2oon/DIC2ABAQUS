function [J,K,KI,KII,Direction]=PlotKorJ(saveto,E,offset,pp)
% fole ris the folder where all the data is vaed
% E is the elasric modulus or the stifness tensor in Pa
% offset: is when abaqus consider your data in meters or standard SI units
% and your data is actually in mm (offset = 1000) or nm (offest = 1e6) or
% in Km (offset = 1e-3) and so on
set(0,'defaultAxesFontSize',25); 	set(0,'DefaultLineMarkerSize',15)
%% input values
if size(E) ~=1      % elastic modulus in Pa
    E11 = E(1,1);                   E12 = E(1,2);
    nu  = E12/(E11+ E12);            E  = E11*(1-2*nu)*(1+nu)/(1-nu);
end
if E<1e6
    if ~exist('pp','var')
        disp('Check the Youngs Modulus units, should be in Pa');
    end
    E = E*1e9;
end
%
% if byCode == 'Y' % if the processing went through without proplems
%     dataum = importdata(fullfile(saveto, 'KJ_Output.txt'));
%     J.Raw      = dataum.data(1,:)*offset*1e-3;
%     KI.RAw     = abs(dataum.data(2,:)*sqrt(offset)*1e-6);
%     KII.Raw    = abs(dataum.data(3,:)*sqrt(offset)*1e-6);
%     K.Raw   = sqrt(abs(J.Raw)*E)*1e-6;

[J.Raw,KI.Raw, KII.Raw, J.K.Raw,Direction.Raw] = ...
    readDATAbaqus([saveto '.dat']);
[saveto,Ond] = fileparts(saveto);
J.Raw     = J.Raw(:)./offset;             % in J/m^2
J.K.Raw   = J.K.Raw(:)./offset;           % in J/m^2
KI.Raw    = KI.Raw(:)./offset^1.5*1e-6;     % in MPa
KII.Raw   = KII.Raw(:)./offset^1.5*1e-6;    % in MPa
K.J_K.Raw = sqrt(J.K.Raw(:)*E)*1e-6;            % in MPa
K.J.Raw   = sqrt(abs(J.Raw(:))*E)*1e-6;         % in MPa
LENK = min(length(KII.Raw), length(KI.Raw));
K.I_II.Raw = sqrt(KII.Raw(1:LENK).^2+KI.Raw(1:LENK).^2);

%% remove outliers
contrs   = length(KI.Raw);        contrs = contrs - round(contrs*0.4);
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<1;   dic = 1; end
J.true   = round(mean(rmoutliers(J.Raw(contrs:end))),dic);
J.div    = round(std(rmoutliers(J.Raw(contrs:end)),1),dic);
J.K.true = round(mean(rmoutliers(J.K.Raw(contrs:end))),dic);
J.K.div  = round(std(rmoutliers(J.K.Raw(contrs:end)),1),dic);

K.J_K.true   = round(mean(rmoutliers(K.J_K.Raw(contrs:end))),dic);
K.J_K.div    = round(std(rmoutliers(K.J_K.Raw(contrs:end)),1),dic);
K.J.true     = round(mean(rmoutliers(K.J.Raw(contrs:end))),dic);
K.J.div      = round(std(rmoutliers(K.J.Raw(contrs:end)),1),dic);
K.I_II.true  = round(mean(rmoutliers(K.I_II.Raw(contrs:end))),dic);
K.I_II.div   = round(std(rmoutliers(K.I_II.Raw(contrs:end)),1),dic);

KI.true  = round(mean(rmoutliers(KI.Raw(contrs:end))),dic);
KI.div   = round(std(rmoutliers(KI.Raw(contrs:end)),1),dic);
KII.true = round(mean(rmoutliers(KII.Raw(contrs:end))),dic);
KII.div  = round(std(rmoutliers(KII.Raw(contrs:end)),1),dic);

if ~isempty(Direction.Raw)
    Direction.true = round(mean(rmoutliers(Direction.Raw(contrs:end))),1);
    Direction.div  = round(std(rmoutliers(Direction.Raw(contrs:end)),1),1);
end
%% Plotting
if ~exist('pp','var'); close all;
    %% for J
    plot(J.Raw,'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');  hold on
    plot(J.K.Raw,'b--s','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b');hold off
    xlabel('Contour Number');       ylabel('J [J/m^2]')
    if min([J.K.Raw(:); J.Raw(:)])>0;     ylim([0 inf]);      end
    title ({['J_{integral} = ' num2str(J.true)   ' ± ' num2str(J.div)   ' J/m^2' ];...
        ['J_{SIF} = '      num2str(J.K.true) ' ± ' num2str(J.K.div) ' J/m^2']});
    legend('J_{integral}','J_{SIFs}','location','best','box','off');
    set(gcf,'WindowStyle','normal');
    set(gcf,'position',[600,100,950,850]);
    axis tight;     xlim([0 length(J.Raw)+2]);
    box off; saveas(gcf, [saveto '\' Ond '_J.fig']);
    saveas(gcf, [saveto '\' Ond '_J.tif']); close;
    
    %% for K
    plot(K.J_K.Raw,'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r'); hold on
    plot(K.J.Raw,'b--s','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b');
    plot(K.I_II.Raw,'k--<','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
    xlabel('Contour Number');       ylabel('K (MPa\surdm)'); hold off;
    title ({['K_{SIF} = '  num2str(K.J_K.true)  ' ± ' num2str(K.J_K.div)  ' MPa\surdm' ]; ...
        ['K_{J} = '    num2str(K.J.true)    ' ± ' num2str(K.J.div)    ' MPa\surdm' ]; ...
        ['K_{I-II} = ' num2str(K.I_II.true) ' ± ' num2str(K.I_II.div) ' MPa\surdm' ]});
    set(gcf,'WindowStyle','normal');    ylim([0 inf]);
    legend('K_{eff-SIFs}','K_{eff-J-integral}','K_{eff-I,II}','location','best','box','off');
    set(gcf,'position',[600,20,950,950]);   xlim([0 length(J.Raw)+2])
    box off; saveas(gcf, [saveto '\' Ond '_Keffs.fig']);
    saveas(gcf, [saveto '\' Ond '_Keffs.tif']); close;
    
    %% for K1 and K2
    plot(KI.Raw,'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');  hold on;
    plot(KII.Raw,'b--s','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b'); hold off
    xlabel('Contour Number');       ylabel('K (MPa m^{1/2})')
    title ({['K_{I} = '  num2str(KI.true)   ' ± ' num2str(KI.div)  ' MPa\surdm' ];...
        ['K_{II} = ' num2str(KII.true)  ' ± ' num2str(KII.div) ' MPa\surdm' ]});
    legend('K_{I}','K_{II}','location','best','box','off');
    set(gcf,'WindowStyle','normal');
    if KII.true>0 && KI.true>0;     ylim([0 inf]);      end
    set(gcf,'position',[600,100,950,850]);  xlim([0 length(J.Raw)+2])
    box off; saveas(gcf, [saveto '\' Ond '_KI and KII.fig']);
    saveas(gcf, [saveto '\' Ond '_KI and KII.tif']);    close
    
    %% plot ALL K
    plot(K.J_K.Raw,'k--o','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');hold on;
    plot(KI.Raw,'r--s','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
    plot(KII.Raw,'b--<','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b'); hold off;
    
    xlabel('Contour Number');       ylabel('K_{abs} (MPa m^{0.5})')
    title ({['K_{eff} = ' num2str(K.J_K.true) ' ± ' num2str(K.J_K.div) ' MPa\surdm' ];...
        ['K_{I} = '   num2str(KI.true)    ' ± ' num2str(KI.div)    ' MPa\surdm' ];...
        ['K_{II} = '  num2str(KII.true)   ' ± ' num2str(KII.div)   ' MPa\surdm' ]});
    legend('K_{eff}','K_{I}','K_{II}','location','best','box','off');
    set(gcf,'WindowStyle','normal'); set(gcf,'position',[600,50,1000,920]);
    axis tight;     xlim([0 length(J.Raw)+2]);
    if min([KII.Raw(:); KI.Raw(:)])>0;     ylim([0 inf]);      end
    box off; saveas(gcf, [saveto '\' Ond '_Keff, KI and KII.fig']);
    saveas(gcf, [saveto '\' Ond '_Keff, KI and KII.tif']); close
    
    %%
    fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
    yyaxis left;    hold on;
    plot(KI.Raw,'k--o','MarkerEdgeColor','k','LineWidth',4);
    plot(KII.Raw,'k--s','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
    ylabel('K (MPa m^{0.5})'); hold off
    if min([KII.Raw(:); KI.Raw(:)])>0;     ylim([0 inf]);      end
    yyaxis right;
    plot(J.K.Raw,'r--<','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
    ylabel('J [J/m^2]');        ylim([0 inf]);
    xlabel('Contour Number');
    legend(['K_{I} = '        num2str(KI.true)  ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
        ['K_{II} = '       num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
        ['J_{integral} = ' num2str(J.K.true) ' ± ' num2str(J.K.div) ' J/m^2'],...
        'location','northoutside','box','off');
    set(gcf,'WindowStyle','normal');
    set(gcf,'position',[60,10,850,990]);  xlim([0 length(J.Raw)+2])
    box off; saveas(gcf, [saveto '\' Ond '_KI, KII and J.fig']);
    saveas(gcf, [saveto '\' Ond '_KI, KII and J.tif']);    %close
    
    %%
    if ~isempty(Direction.Raw)
        Direction.Raw = Direction.Raw(1:length(KII.Raw));
        fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
        yyaxis left;    hold on;
        plot(KI.Raw,'k--o','MarkerEdgeColor','k','LineWidth',4);
        plot(KII.Raw,'k--s','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
        ylabel('K (MPa m^{0.5})'); hold off
        if min([KII.Raw(:); KI.Raw(:)])>0;     ylim([0 inf]);      end
        yyaxis right;
        plot(Direction.Raw,'r--<','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
        ylabel('\theta^{o}'); ylim([-90 90]);
        xlabel('Contour Number');
        legend(['K_{I} = '     num2str(KI.true)  ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
            ['K_{II} = '    num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
            ['Direction = ' num2str(round(mean(rmoutliers(Direction.Raw(contrs:end))),dic))...
            ' ± ' num2str(round(std(rmoutliers(Direction.Raw(contrs:end)),1),dic)) '^o'],...
            'location','northoutside','box','off');
        set(gcf,'WindowStyle','normal');
        set(gcf,'position',[60,10,850,990]);  xlim([0 length(J.Raw)+2])
        box off; saveas(gcf, [saveto '\' Ond '_KI, KII and D.fig']);
        saveas(gcf, [saveto '\' Ond '_KI, KII and D.tif']);    %close
    end
    
else
    if pp==99
        fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
        yyaxis left;    hold on;
        plot(KI.Raw,'k--o','MarkerEdgeColor','k','LineWidth',4);
        plot(KII.Raw,'k--s','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
        ylabel('K (MPa m^{0.5})'); hold off
        if min([KII.Raw(:); KI.Raw(:)])>0;     ylim([0 inf]);      end
        yyaxis right;
        plot(J.Raw,'r--<','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
        ylabel('J [J/m^2]');        ylim([0 inf]);
        xlabel('Contour Number');
        legend(['K_{I} = '        num2str(KI.true)  ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
            ['K_{II} = '       num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
            ['J_{integral} = ' num2str(J.true)   ' ± ' num2str(J.div)   ' J/m^2'],...
            'location','northoutside','box','off');
        set(gcf,'WindowStyle','normal');
        set(gcf,'position',[60,10,850,1050]);  xlim([0 length(J.Raw)+2])
        box off; saveas(gcf, [saveto '\' Ond '_KI, KII and J.fig']);
        saveas(gcf, [saveto '\' Ond '_KI, KII and J.tif']);    close
        
        if ~isempty(Direction.Raw)
            Direction.Raw = Direction.Raw(1:length(KII.Raw));
            fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
            yyaxis left;    hold on;
            plot(KI.Raw,'b--o','MarkerEdgeColor','b','LineWidth',1.5,'MarkerFaceColor','b');
            plot(KII.Raw,'g--s','MarkerEdgeColor','g','LineWidth',1.5,'MarkerFaceColor','g');
            ylabel('K (MPa m^{0.5})'); hold off
            if min([KII.Raw(:); KI.Raw(:)])>0;         end
            yyaxis right;
            plot(Direction.Raw,'r--<','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
            ylabel('\theta^{o}'); ylim([-90 90]);
            xlabel('Contour Number');
            legend(['K_{I} = '     num2str(KI.true)  ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
                ['K_{II} = '    num2str(KII.true) ' ± ' num2str(KII.div) ' MPa\surdm' ],...
                ['Direction = ' num2str(round(mean(rmoutliers(Direction.Raw(contrs:end))),...
                dic)) ' ± ' num2str(round(std(rmoutliers(Direction.Raw(contrs:end)),1),...
                dic)) '^o'],'location','northoutside','box','off');
            set(gcf,'WindowStyle','normal');
            set(gcf,'position',[60,10,850,990]);  xlim([0 length(J.Raw)+2])
            box off; saveas(gcf, [saveto '\' Ond '_KI, KII and D.fig']);
            saveas(gcf, [saveto '\' Ond '_KI, KII and D.tif']);    close
        end
    end
end

