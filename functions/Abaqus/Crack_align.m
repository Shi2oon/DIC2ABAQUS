function [alldat] = Crack_align(alldat)

    %% decide to continue or not
    opts.Interpreter = 'tex'; % Include the desired Default answer
    opts.Default     = 'Y';     % Use the TeX interpreter to format the question
    quest            = 'Is the crack on the x axis?';
    reply           = questdlg(quest,'Boundary Condition','Y','N', opts);
    
if isempty(reply)
    reply = 'Y';
end
if(reply=='N')
    alldata = [alldat.X1(:) alldat.Y1(:) alldat.Ux(:) alldat.Uy(:)];
    Tri(1:2,:)  = flipud(alldata(:,1:2)');
    Tri(3:4,:)  = flipud(alldata(:,3:4)');

    Tri = sortrows(Tri',[1 2])';
    alldata   = Tri';
    [~,alldat]   = reshapeData(alldata);
    
    close all;    
    
imagesc(alldat.X1(1,:),alldat.Y1(:,1),alldat.Uy);      
set(gca,'Ydir','normal');	axis image;
title('Answer in the command line');
xlabel('X','FontSize',20,'FontName','Times New Roman');          
ylabel('Y','FontSize',20,'FontName','Times New Roman');
set(gcf,'position',[30 50 1300 950]);
end