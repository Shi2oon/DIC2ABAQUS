function [Ang,Jv,JE,Mv,ME,KIv,KIE,KIIv,KIIE]=DualJ(saveTo,OfBy,Ang,ToS)
% a code to calculate J, K-I, K-II for a specific gemoetry using declerared
% angles (Ang) as vector for virtual crack extension (q).
% The code takes:
% SaveTo: The .inp file which was used at least once for 1 q direction to
%           get the desired results.
% Stifness: anistropic stifness matrix or young modulus
% OfBy: for the units conversion to meter, 1 for meter, 1e-3 for mm, 1e-6 for um
% Ang: Range of desired angles for q, the range is + from 0 to 180 and
%       negative for 181 (-179) to 359 (-1)
% ToS: used if decleared angeles are reversed in Y direction.

% Written by Abdalrhaman Koko (abdo.koko@materials.ox.ac.uk) as part of 
% DPhil project
% last edit: 17 Dec. 2020 (Abdo)

warning off; iv=0; 
if exist([fileparts(saveTo) '\0_Data.mat'],'file')
    load([fileparts(saveTo) '\0_Data.mat'],'Ang','Jv','JE','KIv',...
         'KIE','KIIv','KIIE','Mv','ME');
else
    %%
for IO = Ang
    %%
    fclose all; clc;
    iv = 1+iv;
  if ~exist([fileparts(saveTo) '\Number_' num2str(iv) '.dat'],'file')
    fid = fopen([erase(saveTo, '.inp') '.inp'], 'r');
    C = textscan(fid, '%s', 'Delimiter', '\n');     fclose(fid);
    S = regexp(C{1}{end-6}, ',', 'split');
    setN = erase(S{1},'_PickedSet');    setN = str2double(setN);
    if isnan(setN)
        setN = erase(S{1},'_PICKEDSET');    setN = str2double(setN);
    end
    C{1}{end-6} = ['_PickedSet' num2str(setN) ', _PickedSet' num2str(setN+1)...
        ', ' num2str(cosd(IO)) ', ' num2str(sind(IO))  ', 0.'];
    C{1}{end-1} = ['_PickedSet' num2str(setN) ', _PickedSet' num2str(setN+1)...
        ', ' num2str(cosd(IO)) ', ' num2str(sind(IO))  ', 0.'];
%     C{1}(end-5:end-1)=[];
    finalform = C{1}(1:end);
    fileID = fopen([fileparts(saveTo) '\Job-' num2str(iv) '.inp'],'w');
    for i=1:size(finalform,1);          stri = finalform(i);
        if ~cellfun('isempty',stri);    fprintf(fileID,'%s\n',char(stri)); end
    end
    fclose(fileID); 
%     if rem(iv,10)==0
%         folder = Js_Code10(iv,fileparts(saveTo));
%     else
        folder = Js_Code(iv,fileparts(saveTo));
%     end
    disp(['Going to ' num2str(iv+1)]);
  else
        folder = [fileparts(saveTo) '\Number_' num2str(iv)];
  end
    try
       	[J{iv},~,KI{iv},KII{iv}] = PlotKorJ(folder,0,OfBy,1);
        Jv(iv)  = J{iv}.true;                   JE(iv)  = J{iv}.div;
    catch
        Jv(iv)  = NaN;                          JE(iv)  = NaN;
    end
    try
        Jv(iv)  = J{iv}.true;                   JE(iv)  = J{iv}.div;        % J inegral
        Mv(iv)  = J{iv}.K.true;                 ME(iv)  = J{iv}.K.div;        % M inegral
        KIv(iv) = KI{iv}.true;                  KIE(iv) = KI{iv}.div;
        KIIv(iv)= KII{iv}.true;                 KIIE(iv)= KII{iv}.div;
    catch
        Mv(iv)  = NaN;                          ME(iv)  = NaN;
        KIv(iv) = NaN;                          KIE(iv) = NaN;
        KIIv(iv)= NaN;                          KIIE(iv)= NaN;
    end
    
    if Ang(iv)<0; Ang(iv) = 360+Ang(iv); end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete %%%%%%%%%%%%%%%%%%%%%%%%%%%
       delete([fileparts(saveTo) '\Job-' num2str(iv) '.inp']);
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.com']);
   if rem(iv,10)~=0
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.inp']);
   else
       try
            PlotKorJ([fileparts(saveTo) '\Number_' num2str(iv)],0,OfBy,99);
       end
   end
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.log']);
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.msg']);
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.odb']);
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.prt']);
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.sim']);
       delete([fileparts(saveTo) '\Number_' num2str(iv) '.sta']);
       clearvars -except J Jv JE iv IO Stiffness Ang saveTo OfBy ToS ...
                         KIv KIE KIIv KIIE setN Mv ME KJv KJE KI_IIv KI_IIE
end
    %% check if angles are actually for y reversed plot 
    if exist('ToS','var');      Ang=ToS;    end
    TableS = table([1:length(Ang)]',Ang(:),Jv(:),JE(:),Mv(:),ME(:),...
                    KIv(:),KIE(:),KIIv(:),KIIE(:),'VariableNames',...
                    {'No','Angle','J_m','STD_J_m','M_int','M_int_STD',...
                   'KI_MPa','STD_KI_MPa','KII_MPa','STD_KII_MPa'});
     writetable(TableS,[fileparts(saveTo) '\0_Data.xlsx'],'Sheet','All_Angles')
     save([fileparts(saveTo) '\0_Data.mat'],'Ang','Jv','JE','KIv',...
         'KIE','KIIv','KIIE','Mv','ME');
end

%%
if exist('ToS','var'); Ang=ToS; end
kol = isoutlier(Jv,'movmedian',15);         Jv(kol)     = NaN;	
Jv(isnan(Ang))   = [];                   	JE(isnan(Ang))   = [];  
kol = isoutlier(Mv,'movmedian',15);         Mv(kol)     = NaN; 
Mv(isnan(Ang))   = [];                      ME(isnan(Ang))   = [];  
kol = isoutlier(KIv,'movmedian',15);        KIv(kol)    = NaN; 
KIv(isnan(Ang))  = [];                      KIE(isnan(Ang))  = [];  
kol = isoutlier(KIIv,'movmedian',15);       KIIv(kol)   = NaN; 
KIIv(isnan(Ang)) = [];                      KIIE(isnan(Ang)) = [];  
Ang(isnan(Ang)) = [];

%%
Para = {'J','KI','KII'};
Unit = {'J [J/m^2]','KI [MPa\surdm]','KII [MPa\surdm]'};

for iK=1:length(Para)
    Tol = 1;
    eval(sprintf('meanJv = sort(%sv);',Para{iK}));%lowest to highest
    meanJv(isnan(meanJv))=[];
    eval(sprintf('STDmin = %sE(ind2sub(size(%sv),find(%sv==meanJv(Tol))));'    ...
        ,Para{iK},Para{iK},Para{iK}));
    eval(sprintf('STDmax = %sE(ind2sub(size(%sv),find(%sv==meanJv(end-Tol))));'...
        ,Para{iK},Para{iK},Para{iK}));
    while std([meanJv(1)   meanJv(Tol+1)])<=STDmin(1) &&...
          std([meanJv(end) meanJv(end-Tol)])<=STDmax(1)
        Tol=Tol+1;
    end
    tol(iK)=Tol;
end
Tol = round(mean(rmoutliers(tol)),0);
if Tol > ceil(length(Ang)*0.015)
    Tol = floor(mean([Tol length(Ang)*0.015]));    
end

%% plot
close all; clc; set(0,'defaultAxesFontSize',35); 	set(0,'DefaultLineMarkerSize',15) 
for iK=1:length(Para)
    for OO=2 %1 for min and 2 for max
eval(sprintf('errorbar(Ang,%sv,%sE,"ok","MarkerSize",10,"MarkerEdgeColor","k","MarkerFaceColor","k");',...
    Para{iK},Para{iK}));
box off;    xlabel('q^{\circ}');  ylabel(Unit{iK})
set(gcf,'position',[20,50,1000,920]);  xlim([min(Ang)-5 max(Ang)+5]); 
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
saveas(gcf, [fileparts(saveTo) '\0_' Para{iK} '.fig']);      
saveas(gcf, [fileparts(saveTo) '\0_' Para{iK} '.tif']); close

%% plot KI,KII, J
% eval(sprintf('meanJv = nanmean(abs(%sv));',Para{iK}));
eval(sprintf('meanJv = sort(%sv);',Para{iK}));%lowest to highest
meanJv(isnan(meanJv))=[];  
if OO == 1
%     eval(sprintf('MinJv = min(abs(%sv)); ind = ind2sub(size(%sv),find(abs(%sv)<=(MinJv+meanJv*Tol)));',...
%                  Para{iK},Para{iK},Para{iK}));
    for Ito=1:Tol
        eval(sprintf('O = ind2sub(size(%sv),find(%sv==meanJv(Ito)));',Para{iK},Para{iK}));
        ind(1:length(O),Ito)=O;
    end
elseif OO == 2
%     eval(sprintf('MaxJv = max(abs(%sv)); ind = ind2sub(size(%sv),find(abs(%sv)>=(MaxJv-meanJv*Tol)));',...
%                  Para{iK},Para{iK},Para{iK}));
	for Ito=1:Tol
        eval(sprintf('O = ind2sub(size(%sv),find(%sv==meanJv(end+1-Ito)));',Para{iK},Para{iK}));
        ind(1:length(O),Ito)=O;
    end 
end
ind=ind(:);     ind(ind==0)=[];     ind = unique(ind);
    eval(sprintf('Angles{OO}.%s = Ang(ind); Angles{OO}.ind%s = ind;',Para{iK},Para{iK}));   
close all;  th=deg2rad(Ang(:));
R = 1;	x = R*cos(th);	y = R*sin(th);	plot(x,y,'--k','linewidth',0.8);    
hold on;    %th = linspace( deg2rad(Ang(ind(1))), deg2rad(Ang(ind(end))), 100);	
for iS=1:length(ind)
    if iS==1
        plot([0 cosd(Ang(ind(iS)))],[0 sind(Ang(ind(iS)))],'-r','linewidth',5); 
    end
    plot([0 cosd(Ang(ind(iS)))],[0 sind(Ang(ind(iS)))],'-r','linewidth',...
        5,'HandleVisibility','off');
end
OKp = Ang; OKp(isnan(OKp))=[];
plot([0 cosd(OKp(1))],[0 sind(OKp(1))],'--k','linewidth',0.8,'HandleVisibility','off');
plot([0 cosd(OKp(end))],[0 sind(OKp(end))],'--k','linewidth',0.8,'HandleVisibility','off');
plot([-1.1 1.1],[0 0],'-b','linewidth',1.5);
plot([0 0],[-1.1 1.1],'-b','linewidth',1.5);
axis image; hold off; set(gcf,'position',[20,50,1000,920]); axis image; box off
if OO == 1
    [~,objh]=legend('q_{Available Range}',['q_{Min. ' Para{iK} '}'],...
                    'Axis','location','best','box','off','fontsize',35); axis off
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12); %// set marker size as desired
    set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
    y=title(num2str(Ang(ind))); y.FontSize = 14;
    saveas(gcf, [fileparts(saveTo) '\0_' Para{iK} '_Min_Range.fig']);      
    saveas(gcf, [fileparts(saveTo) '\0_' Para{iK} '_Min_Range.tif']); close
elseif OO == 2
    [~,objh]=legend('q_{Available Range}',['q_{Max. ' Para{iK} '}'],...
                    'Axis','location','best','box','off','fontsize',35); axis off
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12); %// set marker size as desire
    set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
    y=title(num2str(Ang(ind))); y.FontSize = 14;
    saveas(gcf, [fileparts(saveTo) '\0_' Para{iK} '_Max_Range.fig']);      
    saveas(gcf, [fileparts(saveTo) '\0_' Para{iK} '_Max_Range.tif']); close
end

clearvars -except Tol J Jv JE Mv ME iK OO Ang saveTo KIv KIE KIIv ...
    KIIE Para Unit Angles
    end
end

%% plot KI
for OO=2%1 for min and 2 for max
close all;  th=deg2rad(Ang(:));
R = 1;	x = R*cos(th);	y = R*sin(th);	plot(x,y,'--k','linewidth',0.8);    
hold on; % plotiing J
for iS=1:length(Angles{OO}.indJ)
    if iS==1
        plot([0 cosd(Ang(Angles{OO}.indJ(iS)))],[0 sind(Ang(Angles{OO}.indJ(iS)))],'-r','linewidth',5); 
    end
    plot([0 cosd(Ang(Angles{OO}.indJ(iS)))],[0 sind(Ang(Angles{OO}.indJ(iS)))],'-r',...
        'linewidth',5,'HandleVisibility','off');
end
hold on; % plotiing KI
for iS=1:length(Angles{OO}.indKI)
    if iS==1
        plot([0 cosd(Ang(Angles{OO}.indKI(iS)))],[0 sind(Ang(Angles{OO}.indKI(iS)))],'-b','linewidth',5); 
    end
    plot([0 cosd(Ang(Angles{OO}.indKI(iS)))],[0 sind(Ang(Angles{OO}.indKI(iS)))],'-b',...
        'linewidth',5,'HandleVisibility','off');
end
hold on;   % plotiing KII
for iS=1:length(Angles{OO}.indKII)
    if iS==1
        plot([0 cosd(Ang(Angles{OO}.indKII(iS)))],[0 sind(Ang(Angles{OO}.indKII(iS)))],'-g','linewidth',5); 
    end
    plot([0 cosd(Ang(Angles{OO}.indKII(iS)))],[0 sind(Ang(Angles{OO}.indKII(iS)))],'-g',...
        'linewidth',5,'HandleVisibility','off');
end
OKp = Ang; OKp(isnan(OKp))=[];
plot([0 cosd(OKp(1))],[0 sind(OKp(1))],'--k','linewidth',0.8,'HandleVisibility','off');
plot([0 cosd(OKp(end))],[0 sind(OKp(end))],'--k','linewidth',0.8,'HandleVisibility','off');
plot([-1.1 1.1],[0 0],'-k','linewidth',1.5);
plot([0 0],[-1.1 1.1],'-k','linewidth',1.5);
axis image; hold off; set(gcf,'position',[20,50,1000,920]); axis image; box off

    if OO==1
        [~,objh]=legend('q_{Available Range}','q_{Min. J}','q_{Min. KI}',...
                'q_{Min. KII}','Axis','location','best','box','off','fontsize',35); axis off
        objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
        set(objhl, 'Markersize', 12); %// set marker size as desire
        saveas(gcf, [fileparts(saveTo) '\0_All_Min_Range.fig']);      
        saveas(gcf, [fileparts(saveTo) '\0_All_Min_Range.tif']); close
    elseif OO==2
        [~,objh]=legend('q_{Range}','q_{Max. J}','q_{Max. KI}',...
                'q_{Max. KII}','location','best','box','off','fontsize',35); axis off
        objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
        set(objhl, 'Markersize', 12); %// set marker size as desire
        saveas(gcf, [fileparts(saveTo) '\0_All_Max_Range.fig']);      
        saveas(gcf, [fileparts(saveTo) '\0_All_Max_Range.tif']); close
    end
    clearvars -except Tol J Jv JE Mv ME iK OO Ang saveTo KIv KIE KIIv KIIE ...
        Para Unit meanJv meanKIv meanKIIv Angles
end
%%
fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
yyaxis left;     
errorbar(Ang,KIv,KIE,'bo','MarkerEdgeColor','b','LineWidth',1.5,...
         'MarkerFaceColor','b','MarkerSize',10);   hold on;
errorbar(Ang,KIIv,KIIE,'gs','MarkerEdgeColor','g','LineWidth',1.5,...
        'MarkerFaceColor','g','MarkerSize',10); hold off
ylabel('K (MPa m^{0.5})')
yyaxis right;   
errorbar(Ang,Jv,JE,'r<','MarkerEdgeColor','r','LineWidth',1.5,...
         'MarkerFaceColor','r','MarkerSize',10);  
ylabel('J [J/m^2]')
box off;    xlabel('q^{\circ}');        
[~,objh]=legend('K_{I}','K_{II}','J_{integral}','location','northoutside','box',...
    'off','Orientation','horizontal','fontsize',35); 
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 16); %// set marker size as desire
set(gcf,'WindowStyle','normal');        box off; 
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
set(gcf,'position',[20,20,1200,900]);  xlim([min(Ang)-5 max(Ang)+5]);
saveas(gcf, [fileparts(saveTo) '\0_All_Range_J_KI_II.fig']);
saveas(gcf, [fileparts(saveTo) '\0_All_Range_J_KI_II.tif']);    close 

%%
fig=figure;set(fig,'defaultAxesColorOrder',[[0 0 0]; [1 0 0]]);
yyaxis left;     
errorbar(Ang,KIv,KIE,'ko','MarkerEdgeColor','k','LineWidth',1.5,...
         'MarkerFaceColor','k','MarkerSize',10);   hold on;
errorbar(Ang,KIIv,KIIE,'ko','MarkerEdgeColor','k','LineWidth',1.5,...
        'MarkerSize',10); hold off
ylabel('K (MPa m^{0.5})')
yyaxis right;   
errorbar(Ang,Jv,JE,'r<','MarkerEdgeColor','r','LineWidth',1.5,...
         'MarkerFaceColor','r','MarkerSize',10);  hold on;
errorbar(Ang,Mv,ME,'rs','MarkerEdgeColor','r','LineWidth',1.5,...
         'MarkerSize',10);  hold off;
ylabel('J [J/m^2]')
box off;    xlabel('q^{\circ}');        
[~,objh]=legend('K_{I}','K_{II}','J_{integral}','J_{+interaction}','location',...
    'northoutside','box','off','Orientation','horizontal','fontsize',35); 
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 16); %// set marker size as desire
set(gcf,'WindowStyle','normal');        box off; 
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
set(gcf,'position',[20,20,1200,900]);  xlim([min(Ang)-5 max(Ang)+5]);
saveas(gcf, [fileparts(saveTo) '\0_All_Range_J_M_KI_II.fig']);
saveas(gcf, [fileparts(saveTo) '\0_All_Range_J_M_KI_II.tif']);    close 

%%   
errorbar(Ang,Jv,JE,'r<','MarkerEdgeColor','r','LineWidth',1.5,...
         'MarkerFaceColor','r','MarkerSize',10);   hold on;
errorbar(Ang,Mv,ME,'b>','MarkerEdgeColor','b','LineWidth',1.5,...
         'MarkerFaceColor','b','MarkerSize',10);   hold off;
ylabel('J [J/m^2]')
box off;    xlabel('q^{\circ}');        
[~,objh]=legend('J_{integral}','J_{+interaction}','location',...
    'northoutside','box','off','Orientation','horizontal','fontsize',35); 
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 16); %// set marker size as desire
set(gcf,'WindowStyle','normal'); box off;
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
set(gcf,'position',[20,20,1200,900]);  xlim([min(Ang)-5 max(Ang)+5]);
saveas(gcf, [fileparts(saveTo) '\0_All_Range_J_M.fig']);
saveas(gcf, [fileparts(saveTo) '\0_All_Range_J_M.tif']);    close 

%%
% Maximum Circumferential stress MSC 
subplot(2,2,1); OO=2;
MSC = KIv.*sind(Ang)+KIIv.*(3*cosd(Ang)-1);
plot(Ang,MSC,'.r','markersize',22); hold on
plot(Ang(Angles{OO}.indJ),MSC(Angles{OO}.indJ),'kd','markersize',12,'markerface','k','HandleVisibility','off');
[~,ind_MSC] = min(abs(MSC-0)); % MSC = 0
plot(Ang(ind_MSC),MSC(ind_MSC),'kp','markersize',12,'markerface','w','HandleVisibility','off');
hold off; title('Maximum Circumferential Stress (MSC)','fontsize',25); box off
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

subplot(2,2,2);
KMv = (KIv.^2+KIIv.^2).^0.5;
MSCe = (KIv./KMv).*cos(Ang./2).^3-1.5*(KIIv./KMv).*cos(Ang./2).*sin(Ang);
plot(Ang,MSCe,'.b','markersize',22);hold on
plot(Ang(Angles{OO}.indJ),MSCe(Angles{OO}.indJ),'kd','markersize',12,'markerface','k','HandleVisibility','off');
[~,ind_MSCe] = min(abs(MSCe-1));
plot(Ang(ind_MSCe),MSCe(ind_MSCe),'kp','markersize',12,'markerface','w','HandleVisibility','off'); 
hold off; legend('MSCE','location','best','fontsize',25,'box','off');
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

subplot(2,2,3); 
% minimum strain energy density (MSED) 
% kappa = 3 - (4 .* nu); % [/]'plane_strain'
nu=0.3;
kappa = (3 - nu)./(1 + nu); % [/]'plane_stress'
a11 = (1+cosd(Ang)).*(kappa-cosd(Ang))./16;
a12 = (sind(Ang)/16).*(2*cosd(Ang)-(kappa-1));
a22 = ((kappa+1).*(1-cosd(Ang))+(1+cosd(Ang).*(3*cosd(Ang)-1)))./16;
MSED = 1/(kappa-1)*(a11.*(KIv./KMv).^2+2*a12.*(KIv.*KIIv./KMv.^2)+a22.*(KIIv./KMv).^2);
plot(Ang,MSED,'.r','markersize',22); hold on
plot(Ang(Angles{OO}.indJ),MSED(Angles{OO}.indJ),'kd','markersize',12,'markerface','k','HandleVisibility','off');
MSED(1:ceil(length(Ang)/6))   = NaN;
MSED(end:ceil(length(Ang)/6)) = NaN;
kol = isoutlier(MSED,'movmedian',1);      MSED(kol)     = NaN;    
[~,ind_MSED] = min(MSED);
plot(Ang(ind_MSED),MSED(ind_MSED),'kp','markersize',12,'markerface','w','HandleVisibility','off'); 
hold off; title('Minimum Strain Energy Density (MSED)','fontsize',25); box off
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

subplot(2,2,2); 
% maximum energy release rate (MERR) 
MERR = (kappa+1)/8*(KIv.^2.*(sind(Ang/2)+sind(3*Ang/2))+4*KIv.*KIIv.*cosd(3*Ang/2)-...
        KIIv.^2.*(3*sind(3*Ang/2)-5*sind(Ang/2)));
plot(Ang,MERR,'.r','markersize',22);     hold on
plot(Ang(Angles{OO}.indJ),MERR(Angles{OO}.indJ),'kd','markersize',12,'markerface','k','HandleVisibility','off');
[~,ind_MERR] = min(abs(MERR-0));
plot(Ang(ind_MERR),MERR(ind_MERR),'kp','markersize',12,'markerface','w','HandleVisibility','off'); 
hold off; %legend('MERR','Direction','Max J','location','best','fontsize',25,'box','off');
box off; title('Maximum Energy Release Rate (MERR-Nuismer)','fontsize',25)
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

subplot(2,2,4);
gamma_MERRe = rad2deg(atan(KIIv./KIv));
MERRe = (kappa+1)/8*(2*sind(3*Ang/2+2*gamma_MERRe)+((1+4*sind(gamma_MERRe).^2).*sind(Ang/2)-sind(3*Ang/2)));
plot(Ang,MERRe,'.r','markersize',22); hold on
plot(Ang(Angles{OO}.indJ),MERRe(Angles{OO}.indJ),'kd','markersize',12,'markerface','k'); box off
[~,ind_MERRe] = min(abs(MERRe-0));
plot(Ang(ind_MERRe),MERRe(ind_MERRe),'kp','markersize',12,'markerface','w');
hold off; legend('MERR-Chang','Max J','Direction','location','best','fontsize',25,'box','off');
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

set(gcf,'position',[1 41 1920 1400]); 
saveas(gcf, [fileparts(saveTo) '\0_theta_Cir.fig']);
saveas(gcf, [fileparts(saveTo) '\0_theta_Cir.tif']);    close 

%%
% Maximum tangential stress criterion
r=1;    Thet_MTSC = rad2deg(acos((3*KIIv.^2+(KIv.^4+8*KIv.^2.*...
                                KIIv.^2).^0.5)./(KIv.^2+9*KIIv.^2))); % MSC
ToRth = (cosd(0.5*Thet_MTSC)/sqrt(2*pi*r)).*(KIv.*sind(Thet_MTSC)+...
         KIIv.*(3*cosd(Thet_MTSC)-1));
Theta_2D_MERR = rad2deg(atan(2*KIv.*KIIv)./(KIv.^2+KIIv.^2));
% Maximum energy release rate criterion, direcion that maximise J
subplot(2,2,1);  plot(Ang,Theta_2D_MERR,'.r','markersize',22,'HandleVisibility','off'); hold on
                 plot(Ang(Angles{OO}.indJ),Theta_2D_MERR(Angles{OO}.indJ),...
                     'dk','markersize',12,'markerface','k');    
                 [~,ind_Theta_2D_MERR] = min(abs(Theta_2D_MERR-0));
                 plot(Ang(ind_Theta_2D_MERR),Theta_2D_MERR(ind_Theta_2D_MERR),'pw','markersize',12,...
                     'markeredgecolor','k','markerface','w');    hold off
	xlabel('q^\circ');  ylabel('\theta^{MERR}');                box off
[~,objh] = legend('Max J','Criteria', 'location','best','fontsize',20,'box','off'); 
objhl = findobj(objh, 'type', 'line');      set(objhl, 'Markersize', 15); 
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

subplot(2,2,2);  plot(Ang,Thet_MTSC,'.r','markersize',22,'HandleVisibility','off');          hold on
                 plot(Ang(Angles{OO}.indJ),Thet_MTSC(Angles{OO}.indJ),...
                      'dk','markersize',12,'markerface','k');   
                 [~,ind_Thet_MTSC] = min(abs(Thet_MTSC-0));
                 plot(Ang(ind_Thet_MTSC),Thet_MTSC(ind_Thet_MTSC),'pw','markersize',12,...
                     'markeredgecolor','k','markerface','w');    hold off
	xlabel('q^\circ');  ylabel('\theta^{MSC}');             	box off
[~,objh] = legend('Max J','Criteria', 'location','best','fontsize',20,'box','off'); 
objhl = findobj(objh, 'type', 'line');      set(objhl, 'Markersize', 15); 
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

subplot(2,2,3);  plot(Ang,gamma_MERRe,'.r','markersize',22,'HandleVisibility','off');          hold on
                 plot(Ang(Angles{OO}.indJ),gamma_MERRe(Angles{OO}.indJ),...
                     'dk','markersize',12,'markerface','k');   
                 [~,ind_gamma_MERRe] = min(abs(gamma_MERRe-0));
                 plot(Ang(ind_gamma_MERRe),gamma_MERRe(ind_gamma_MERRe),'pw','markersize',12,...
                     'markeredgecolor','k','markerface','w');    hold off
	xlabel('q^\circ');  ylabel('tan^{-1}(K_{II}/K_I)');         box off
[~,objh] = legend('Max J','Criteria', 'location','best','fontsize',20,'box','off'); 
objhl = findobj(objh, 'type', 'line');      set(objhl, 'Markersize', 15); 
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);

subplot(2,2,4);  plot(abs(KIv./KMv),abs(KIIv./KMv),...
                      '.r','markersize',22,'HandleVisibility','off'); 
                 hold on;   plot(abs(KIv(Angles{OO}.indJ)./KMv(Angles{OO}.indJ)),...
                                 abs(KIIv(Angles{OO}.indJ)./KMv(Angles{OO}.indJ)),...
                                 'dk','markersize',12,'markerface','k'); hold off
xlabel('K_{I}/K_{IC}');         ylabel('K_{II}/K_{IC}');    axis image;box off
% whichisbigger = max([abs(KIv./KMv) abs(KIIv./KMv)]);
xlim([0 1]);                    ylim([0 1]); 
xticks([0 0.25 0.5 0.75 1]);    yticks([0 0.25 0.5 0.75 1]);
[~,objh] = legend('Max J', 'location','southwest','fontsize',20,'box','off'); 
objhl = findobj(objh, 'type', 'line');      set(objhl, 'Markersize', 15); 

set(gcf,'position',[1 41 1920 1400]); 
saveas(gcf, [fileparts(saveTo) '\0_theta.fig']);
saveas(gcf, [fileparts(saveTo) '\0_theta.tif']);    close 
save([fileparts(saveTo) '\0_Data.mat'],'Angles','OO','saveTo','-append');

%%
Y_Name ={'MSC','MSCe','MSED','MERR','MERRe','Theta_2D_MERR','Thet_MTSC','gamma_MERRe'};
Angv = Ang;     AngE = zeros(size(Ang));
ind_Ang = Angles{OO}.indJ;
X_Data = {'Ang','J','M','KI','KII'};
for iV = 1:length(Y_Name) %creieria
    for iX = 1:length(X_Data) % data
        eval(sprintf('Datav(iV,iX) = %sv(ind_%s); DataE(iV,iX) = %sv(ind_%s);',...
          X_Data{iX},Y_Name{iV},X_Data{iX},Y_Name{iV}));
    end
end
TableS = table(Y_Name',Datav(:,1),Datav(:,2),DataE(:,2),Datav(:,3),...
            DataE(:,3),Datav(:,4),DataE(:,4),Datav(:,5),DataE(:,5),...
            'VariableNames',{'Criterion','Angle','J','STD_J','M',...
            'STD_M','KI_MPa','STD_KI_MPa','KII_MPa','STD_KII_MPa'});
writetable(TableS,[fileparts(saveTo) '\0_Data.xlsx'],'Sheet','Criteria')

%% plot
close all; clc; set(0,'defaultAxesFontSize',35); 	set(0,'DefaultLineMarkerSize',15) 
for iK=1:length(Y_Name)
eval(sprintf('ind = ind_%s;',Y_Name{iK}));
ind=ind(:);     ind(ind==0)=[];     ind = unique(ind);
  
close all;  th=deg2rad(Ang(:));
R = 1;	x = R*cos(th);	y = R*sin(th);	plot(x,y,'--k','linewidth',0.8);    
hold on;    %th = linspace( deg2rad(Ang(ind(1))), deg2rad(Ang(ind(end))), 100);	
for iS=1:length(ind)
    if iS==1
        plot([0 cosd(Ang(ind(iS)))],[0 sind(Ang(ind(iS)))],'-r','linewidth',5); 
    end
    plot([0 cosd(Ang(ind(iS)))],[0 sind(Ang(ind(iS)))],'-r','linewidth',...
        5,'HandleVisibility','off');
end
plot([0 cosd(Ang(1))],[0 sind(Ang(1))],'--k','linewidth',0.8,'HandleVisibility','off');
plot([0 cosd(Ang(end))],[0 sind(Ang(end))],'--k','linewidth',0.8,'HandleVisibility','off');
plot([-1.1 1.1],[0 0],'-b','linewidth',1.5);
plot([0 0],[-1.1 1.1],'-b','linewidth',1.5);
axis image; hold off; set(gcf,'position',[20,50,1000,920]); axis image; box off
    [~,objh]=legend('q_{Available Range}',['q_{' Y_Name{iK} '}'],...
                    'Axis','location','best','box','off','fontsize',35); axis off
    objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
    set(objhl, 'Markersize', 12); %// set marker size as desire
    set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
    y=title(num2str(Ang(ind))); y.FontSize = 14;
    saveas(gcf, [fileparts(saveTo) '\0_Angle_' Y_Name{iK} '.fig']);      
    saveas(gcf, [fileparts(saveTo) '\0_Angle_' Y_Name{iK} '.tif']); close
end

%% my method
try
for iV=1:floor(length(Jv)/2)
    try
        j1 = Jv(iV+90);        
        j2 = Jv(iV);
        func=@(theta) j1*sin(theta)+j2*cos(theta); %find the root
        func=@(theta) j1*cos(theta)+j2*cos(theta);
        x1x2dr(:,iV) = [Ang(iV+90) Ang(iV) ...
            rad2deg(fzero(func,[deg2rad(min(Ang)) deg2rad(max(Ang))]))];
%         if Ang(iV+90) == Ang(Angles{OO}.indJ)
    end
end
% count=1;
% for iV=1:length(Angles{OO}.indJ)
%     if ~isempty(find(Angles{OO}.J(iV)== x1x2(1,:)))
%         Ok = find(Angles{OO}.J(iV)== abs(x1x2(1,:)));
%         oAn(1,count) = Ok(1);
%         oAn(2,count) = Angles{OO}.J(iV);
%         count=count+1;
%     end
% end
plot(x1x2dr(1,:), x1x2dr(3,:),'.r','markersize',22); hold on
plot(x1x2dr(2,:), x1x2dr(3,:),'.b','markersize',22); hold off
% try
%     plot(oAn(2,:),Crack_direction(oAn(1,:)),'dk','markersize',12,'markerface','k');   hold off
% end
xlabel('q^\circ');  ylabel('\theta^{*}');             	box off
set(gca,'Xtick',[min(Ang):30:0, 30:30:max(Ang)]);
[~,objh] = legend('x_1','x_2','Max J','location','best','fontsize',25); 
objhl = findobj(objh, 'type', 'line');      set(objhl, 'Markersize', 16);
title({'J_k^* = R_{kk} J_k,      J_2=0';''},'fontsize',25)
set(gcf,'position',[432 111 900 817]); 
saveas(gcf, [fileparts(saveTo) '\0_theta_my.fig']);
saveas(gcf, [fileparts(saveTo) '\0_theta_my.tif']);    %close 

TableS = table(x1x2dr','VariableNames',{'x1_x2_dir'});
writetable(TableS,[fileparts(saveTo) '\0_Data.xlsx'],'Sheet','My_Criteria')
end
end

%%
function [OtF] = Js_Code(NumB,folder)
fileid = fullfile(folder,'Abaqus_script.py');
fileID = fopen(fileid,'w');
fprintf(fileID,'# -*- coding: mbcs -*- \n');
fprintf(fileID,'from part import * \n');
fprintf(fileID,'from material import * \n');
fprintf(fileID,'from section import * \n');
fprintf(fileID,'from assembly import * \n');
fprintf(fileID,'from step import * \n');
fprintf(fileID,'from interaction import * \n');
fprintf(fileID,'from load import * \n');
fprintf(fileID,'from mesh import * \n');
fprintf(fileID,'from optimization import * \n');
fprintf(fileID,'from job import * \n');
fprintf(fileID,'from sketch import * \n');
fprintf(fileID,'from visualization import * \n');
fprintf(fileID,'from connectorBehavior import * \n');
fprintf(fileID,'mdb.ModelFromInputFile(inputFileName= \n');
SaveTmp = pythonFileName([folder '/Job-' num2str(NumB) '.inp']);
fprintf(fileID,'    "%s",  \n',SaveTmp);
fprintf(fileID,'    name="Job-%d") \n',NumB);
fprintf(fileID,'from part import *\n');
fprintf(fileID,'from material import *\n');
fprintf(fileID,'from section import *\n');
fprintf(fileID,'from assembly import *\n');
fprintf(fileID,'from step import *\n');
fprintf(fileID,'from interaction import *\n');
fprintf(fileID,'from load import *\n');
fprintf(fileID,'from mesh import *\n');
fprintf(fileID,'from optimization import *\n');
fprintf(fileID,'from job import *\n');
fprintf(fileID,'from sketch import *\n');
fprintf(fileID,'from visualization import *\n');
fprintf(fileID,'from connectorBehavior import *\n');
fprintf(fileID,'mdb.Job(atTime=None, contactPrint=OFF, description="", echoPrint=OFF, \n');
fprintf(fileID,'    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, \n');
fprintf(fileID,'    memory=90, memoryUnits=PERCENTAGE, model="Job-%d", modelPrint=OFF, \n',NumB);
fprintf(fileID,'    multiprocessingMode=DEFAULT, name="Number_%d", nodalOutputPrecision=SINGLE, \n',NumB);
fprintf(fileID,'    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch="", type=\n');
fprintf(fileID,'    ANALYSIS, userSubroutine="", waitHours=0, waitMinutes=0)\n');
fprintf(fileID,'mdb.jobs["Number_%d"].submit(consistencyChecking=OFF)\n',NumB);
%{
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STARTED, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "clientHost": "OUMS-TJM10", "handle": 0, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(WARNING, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "message": "THE REQUEST FOR MISES OUTPUT WILL BE REPLACED BY A REQUEST FOR S OUTPUT", \n');
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(ODB_FILE, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "file": "A:\\OneDrive - Nexus365\\Work\\ABAQUS\\Python and ABAQUS\\Dual\\Number_%d.odb", \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(COMPLETED, {"phase": BATCHPRE_PHASE, \n',NumB);
fprintf(fileID,'    "message": "Analysis phase complete", "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STARTED, {"phase": STANDARD_PHASE, \n',NumB);
fprintf(fileID,'    "clientHost": "OUMS-TJM10", "handle": 14604, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STEP, {"phase": STANDARD_PHASE, "stepId": 1, \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(ODB_FRAME, {"phase": STANDARD_PHASE, "step": 0, \n',NumB);
fprintf(fileID,'    "frame": 0, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(MEMORY_ESTIMATE, {"phase": STANDARD_PHASE, \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d", "memory": 35.0})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(ODB_FRAME, {"phase": STANDARD_PHASE, "step": 0, \n',NumB);
fprintf(fileID,'    "frame": 1, "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(STATUS, {"totalTime": 1.0, "attempts": 1, \n',NumB);
fprintf(fileID,'    "timeIncrement": 1.0, "increment": 1, "stepTime": 1.0, "step": 1, \n');
fprintf(fileID,'    "jobName": "Number_%d", "severe": 0, "iterations": 1, \n',NumB);
fprintf(fileID,'    "phase": STANDARD_PHASE, "equilibrium": 1})\n');
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(END_STEP, {"phase": STANDARD_PHASE, "stepId": 1, \n',NumB);
fprintf(fileID,'    "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'mdb.jobs["Number_%d"]._Message(COMPLETED, {"phase": STANDARD_PHASE, \n',NumB);
fprintf(fileID,'    "message": "Analysis phase complete", "jobName": "Number_%d"})\n',NumB);
fprintf(fileID,'# Save by scro3511 on 2020_05_10-17.16.29; build 2016 2015_09_24-21.31.09 126547\n');
%}
fclose(fileID);

%% Excute
PWD =pwd;           cd(folder)
system(['abaqus cae ','noGUI','=Abaqus_script.py']); % Windows system
cd(PWD);
OtF = [folder '\Number_' num2str(NumB)];
end