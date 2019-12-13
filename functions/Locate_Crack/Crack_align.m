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
    alldata = [alldat.X(:) alldat.Y(:) alldat.Ux(:) alldat.Uy(:)];
    Tri(1:2,:)  = flipud(alldata(:,1:2)');
    Tri(3:4,:)  = flipud(alldata(:,3:4)');

    Tri = sortrows(Tri',[1 2])';
    alldata   = Tri';
    alldat   = reshapeData(alldata);
end