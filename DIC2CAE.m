function [J,KI,KII,KIII,Direction] = DIC2CAE(Maps)
% variable M4 is for the displacement components
%%
warning on; addpath([pwd '\functions'])
nn =2;
if ~isfield(Maps,"Uz")
    Maps.Uz = ones(size(Maps.Ux))*1e-12;
    Maps.Z = zeros(size(Maps.X));
    nn =1;
end
if size(Maps.Ux,2)*size(Maps.Ux,2) ~= size(Maps.Uz,2)*size(Maps.Uz,2)
    Maps.Uz = squeeze(Maps.Uz(2:end,2:end));
end
if ~isfield(Maps,"stepsize")
    Maps.stepsize = mean(unique(diff(unique(Maps.Y(:)))));
end

% this can be used to fix the number of converged J values through
% different sets
if ~isfield(Maps,"Alot")
    Maps.Alot=1;
end

%%
if strcmpi(Maps.type, 'A')
    [Maps.E,Maps.nu,Maps.G,Maps.Co] = effectiveE_v(Maps.Stiffness); % in Pa
else
    Maps.G = Maps.E/(2*(1 + Maps.nu));
end

if strcmpi(Maps.stressstat, 'plane_strain')
    Maps.E = Maps.E/(1-Maps.nu^2);% for HR-EBSD plane strain conditions
    Maps.G = Maps.E/(2*(1 + Maps.nu));
end

names = {'KI','KIII'};
for iO=1:nn
    Dirxyz = Maps;
    Dirxyz.unique = names{iO};
    if iO == 1      % Mode I/II
        Dirxyz.Ux = Maps.Ux;
        Dirxyz.Uy = Maps.Uy;
    elseif iO == 2 % Mode III
        Dirxyz.Ux = Maps.Uz;
        % in case it is zero as Abaqus won't work
        Dirxyz.Uy = ones(size(Maps.Ux))*1e-12;
    end
    [Dirxyz,UnitOffset,Dirxyz.msk,SaveD] = ...
        Locate_Crack(Dirxyz);
    if ~isfield(Maps,"msk")
        Maps.msk = Dirxyz.msk;
    end
    % prepare and run abaqus cae
    [Abaqus,~] = PrintRunCode(Dirxyz, ...
        Dirxyz.msk,SaveD,ceil(min(size(Dirxyz.X))*0.5-2),UnitOffset);

    if iO == 1      % Mode I
        [Jd,K,KI,KII,Direction] = PlotKorJ(Abaqus,Maps.E,UnitOffset,Maps.Alot);
        if ~isempty(Direction.Raw)
            fprintf('\nRecommended J-integral direction is %d ± %d\t',...
                round(Direction.true,1), Direction.div)
            Ans = questdlg_timer(60,['J-integral direction is ' ...
                num2str(Direction.true) ' ± ' num2str(Direction.div) ...
                ', Do you want to adjust?'],...
                [ num2str(Direction.true) ' ± ' num2str(Direction.div)],...
                'Y','N','N');
            if Ans == 'Y'
                if (abs(Direction.true(1))-abs(Direction.div(1)))>5 && ...
                        (abs(Direction.true(1))/abs(Direction.div(1)))>5
                    Ans = 'Y';
                    fprintf(' is being corrected\n')
                else
                    Ans = 'N';
                    fprintf(' not corrected due to huge uncertanity in the direction\n')
                end
            end
            if strcmpi(Ans,'Y')
                [Abaqus] = Adjust4Direction(Abaqus,Direction.true);
                OldDirection = Direction;
                [Jd,~,KI,KII,Direction] = PlotKorJ(Abaqus,Maps.E,UnitOffset,1);
                Direction.OldRaw(iO,1:length(OldDirection.Raw))=OldDirection.Raw;
                Direction.Oldtrue(iO)   = OldDirection.true;
                Direction.Olddiv(iO)    = OldDirection.div;
            end
        end
        if mean(isnan(KI.Raw))==1
            KI = K.J;
        end
        loT(iO) = length(KI.Raw);

    elseif iO==2 % fix KIII to shear rather than modulus
        [~,K,~,KIII,~] = PlotKorJ(Abaqus,Maps.E,UnitOffset,Maps.Alot);
        % correct from in-plane to out-of-plane shear
        if strcmpi(Ans,'Y')
            [Abaqus] = Adjust4Direction(Abaqus,Direction.Oldtrue);
            [~,~,~,KIII,~] = PlotKorJ(Abaqus,Maps.E,UnitOffset,Maps.Alot);
        end
        if mean(isnan(KIII.Raw))==1
            KIII = K.J;
        end
        KIII.Raw = KIII.Raw*2*Maps.G/Maps.E;
        Jd.Raw = (KIII.Raw*1e6).^2/(2*Maps.G);
        loT(iO)  = length(KIII.Raw);
    end
    % J when calculating the SIF (more accurate)
    JKRaw(iO,1:length(Jd.K.Raw)) = Jd.Raw;
    JRaw(iO,1:length(Jd.Raw)) = Jd.Raw; % J from J analysis
end
J.JKIII = JKRaw;
J.JIII = JRaw;

%% Cut to the same contour convergence (IoT value)
if nn==1
    J.Raw = J.JIII(:,1:min(loT));
    J.K.Raw = J.JKIII(:,1:min(loT));
else
    J.Raw = sum(J.JIII(:,1:min(loT)));
    J.K.Raw = sum(J.JKIII(:,1:min(loT)));
    KIII.Raw = abs(KIII.Raw(1:min(loT)));
end

J.Raw = J.Raw(1:min(loT));
J.K.Raw = J.K.Raw(1:min(loT));
KI.Raw = KI.Raw(1:min(loT));
KII.Raw = abs(KII.Raw(1:min(loT)));

%%
contrs   = length(J.Raw);        contrs = contrs - round(contrs*0.4);
dic = real(ceil(-log10(nanmean(rmoutliers(J.Raw(contrs:end))))))+2;
if dic<2;       dic = 2;    end
J.true   = round(mean(rmoutliers(J.Raw(contrs:end))),dic);
J.div    = round(std(rmoutliers(J.Raw(contrs:end)),1),dic);
J.K.true   = round(mean(rmoutliers(J.K.Raw(contrs:end))),dic);
J.K.div    = round(std(rmoutliers(J.K.Raw(contrs:end)),1),dic);

KI.true  = round(mean(rmoutliers(KI.Raw(contrs:end))),dic);
KI.div   = round(std(rmoutliers(KI.Raw(contrs:end)),1),dic);
if nn==2
    KIII.true = round(mean(rmoutliers(KIII.Raw(contrs:end))),dic);
    KIII.div  = round(std(rmoutliers(KIII.Raw(contrs:end)),1),dic);
else
    KIII=[];
end
KII.true = round(mean(rmoutliers(KII.Raw(contrs:end))),dic);
KII.div  = round(std(rmoutliers(KII.Raw(contrs:end)),1),dic);

if ~isempty(Direction.Raw)
    Direction.true= round(mean(rmoutliers(Direction.Raw(contrs:end,1))),1);
    Direction.div = round(std(rmoutliers(Direction.Raw(contrs:end,1)),1),1);
    if strcmpi(Ans,'Y')
        Direction.Oldtrue= round(mean(rmoutliers(Direction.OldRaw(contrs:end,1))),1);
        Direction.Olddiv = round(std(rmoutliers(Direction.OldRaw(contrs:end,1)),1),1);
    end
end
%
%%
plotJKIII(KI,KII,KIII,J,Maps.stepsize,Maps.input_unit)
saveas(gcf, [Maps.results '_J_KI_II_III.fig']);
saveas(gcf, [Maps.results '_J_KI_II_III.tif']);    %close all

save([Maps.results '_DIC2CAE.mat'],'Maps','J','KI','KII',...
    'KIII','Direction');
end
