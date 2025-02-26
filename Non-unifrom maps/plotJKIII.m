function plotJKIII(KI,KII,KIII,J,stepsize,input_unit)
if isempty(KIII)
    KII.Raw = abs(KII.Raw);
    KII.true = abs(KII.true);
    Kd = [KI.Raw(:); KII.Raw(:)];
else
    Kd = [KI.Raw(:); KII.Raw(:); KIII.Raw(:)];
end
% if min(Kd(:))>0
%     Kd(end+1,1) = 0;
% end
contrs   = length(J.Raw);        contrs = contrs - round(contrs*0.4);
set(0,'defaultAxesFontSize',22);       set(0,'DefaultLineMarkerSize',16)
Contour = (1:length(J.Raw))*stepsize;
fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
yyaxis left;    hold on;
    patch([Contour(contrs) Contour(end) Contour(end) Contour(contrs)], ...
        [min(Kd(:)) min(Kd(:)), max(Kd(:)) max(Kd(:))]...
       , [0.95 0.5 0.9],'edgecolor','none','facealpha','0.3')
plot(Contour,KI.Raw,'k--o','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
plot(Contour,KII.Raw,'k--s','MarkerEdgeColor','k','LineWidth',1.5,'MarkerFaceColor','k');
if ~isempty(KIII)
    plot(Contour,KIII.Raw,'k--d','MarkerEdgeColor','k','LineWidth',...
        1.5','MarkerFaceColor','k');
end
ylabel('K (MPa m^{0.5})'); hold off
if min(Kd(:))>0;     ylim([0 max(Kd(:))+min(Kd(:))/3]);      end
yyaxis right;
plot(Contour,J.K.Raw,'r--<','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
ylabel('J (J/m^2)');        
ylim([0 inf]);
if isempty(KIII)
    legend( 'Convergent EDI', ...
            ['K_{I} = '     num2str(KI.true)   ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
            ['K_{II} = '       num2str(KII.true)  ' ± ' num2str(KII.div) ' MPa\surdm' ],...
            ['J_{integral} = ' num2str(J.K.true)    ' ± ' num2str(J.K.div)   ' J/m^2'],...
            'location','northoutside','box','off','NumColumns',2,'fontsize',16);
else
legend( 'Convergent EDI',...
        ['K_{I} = '     num2str(KI.true)   ' ± ' num2str(KI.div)  ' MPa\surdm' ],...
        ['K_{II} = '       num2str(KII.true)  ' ± ' num2str(KII.div) ' MPa\surdm' ],...
        ['K_{III} = '      num2str(KIII.true) ' ± ' num2str(KIII.div) ' MPa\surdm' ],...
        ['J_{integral} = ' num2str(J.K.true)    ' ± ' num2str(J.K.div)   ' J/m^2'],...
        'location','northoutside','box','off','NumColumns',2,'fontsize',16);
end
set(gcf,'position',[60,-70,776,1000]);grid on;  box off;
ax1 = gca;  axPos = ax1.Position;
% Change the position of ax1 to make room for extra axes
% format is [left bottom width height], so moving up and making shorter here...
ax1.Position = axPos + [0 0.2 0 -0.15];
% Exactly the same as for plots (above), axes LineWidth can be changed inline or after
ax1.LineWidth = 1;
% Add two more axes objects, with small multiplier for height, and offset for bottom
ax2 = axes('position', (axPos .* [1 1 1 1e-3]) + [0 0.08 0 0], 'color', 'none', 'linewidth', 1);
% You can change the limits of the new axes using XLim
ax2.XLim = [0 length(Contour)+1];     ax1.XLim = [0 max(Contour)+stepsize];
% You can label the axes using XLabel.String
if strcmpi(input_unit,'um')
    input_unit = '\mum';
end
ax1.XLabel.String = ['Contour Length (' input_unit ')'];
ax2.XLabel.String = 'Contour Number';
end