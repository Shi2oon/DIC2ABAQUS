function [el,mesh] = RotRemoval(centre,mesh,el)
% RBMCORR2D removes the 2D (planar) rigid body motions (translation and rotation from 
%     a areal vector dataset (e.g. from DIC). The algorithm for rotation 
%     removal is based on Ken Shoemake's Euler angle correction.
%
% Inputs:
%     dat2D_col       4 column data as [coordinates displacements], as from
%                     DaVis - originally 'cdata'
%     dataSize        2 element vector giving size of dataset [X Y]
%     
% Outputs:
%     dataOut         4 column data corrected for rigid body movement
%     eulerAngles     calculated Euler angles 
%     rotCentre       calculated rotation centre
%     
% Code originally developed by M. Mostafavi, with additions by M. Jordan.
% Based on Ken Shoemake's Euler angle extraction.
%
% Last edit: April 2014 (M. Jordan)
% (C) 2014, University of Oxford

data2D_col = [mesh.UDIC' mesh.dDIC'];

%% 1. Clean the dataset

cdata=data2D_col(~any(isnan(data2D_col),2),:);

if(sum(sum(sum(isnan(cdata)))))~=0
    disp('WARNING: Data contains vectors with partial NaN entries');
end

%% 2. Subtract RBD from displacements
C0=zeros(1,4);
RBD = nanmean(cdata(:,3:4));      %row vector of 2 RBD coordinates
C0(1,3:4)= RBD;
disp(['RBD /pixel = ' num2str(RBD,6)])
cdata = bsxfun(@minus,cdata, C0);
%bsxfun performs same operation as following two lines


%% 3. Check for rotation correction request

if ~strcmpi(centre,'norot')
    %% 4. Finding the centre of rotation and setting as coordinate origin
    % Use one of 4 methods to calculate position of rigid body motion and set
    % as coordinate origin.
    %
    % X0 coords = index coords of rotn centre

    switch lower(centre)
        case 'minimum'
            %For using centre as minimum displacement vector
            [~,X0L]=min(sum(abs(cdata(:,3:4)')));     %minimum vector of cdata
            rotCentre=cdata(X0L,1:2);
            disp('##ISSUE: Improve minimum finding algorithm##')
        case 'true'
            %For using centre after RBD correction
            rotCentre = (nanmean(cdata(~isnan(sum(cdata(:,3:4),2)),1:2)));
    end

    %Set rotn centre as coord origin
    X0 = [rotCentre 0 0];
    cdata = bsxfun(@minus,cdata,X0);

    %xi0 is the location of points before loading and xip is the location after
    %displacement

    xi0=cdata(:,1:2);                   % Reference coords
    xip=cdata(:,1:2)+cdata(:,3:4);      % Comparison dataset coords

    %% 4.1 Calculating rotation matrix, and checking the rank 
    % Rank should equal 2 to avoid singularity
    rank_xi0 = rank(xi0');
    rank_xip = rank(xip');

    if rank_xi0 ~= 2 || rank_xip ~= 2
        disp(['##ERROR:xi0 or xip rank deficient. rank_xi0 = ' ...
            num2str(rank_xi0) ', rank_xip = ' num2str(rank_xip)]); ERROR
    end

    R = xi0\xip;  % backslash is MATLAB matrix division - solves for rotation efficiently

    %% 4.2 Extracting Euler angles from the rotation matrix
    % Theoretical rotation matrix based on the calculated angles:
    % R_theo=   [cos(t), sin(t)
    %            -sin(t), cos(t)];
    theta = atan2(R(1,2),R(1,1));       %initial estimate for theta

    % disp(['t = ' num2str(t,'% 1.3f ')])
    % disp(['Euler angle estimate /degree: (gamma) = (' num2str(radtodeg(t),'% 1.3f ') ').'])

    %% 4.2.1 Initial estimate of Rotation matrix 
    %OC1
    %% 4.2.2 Refine estimate of Euler angles and Rot matrix
    disp '>>Solver details (refining theta):'
    disp '***********'
    options2 = optimoptions('lsqnonlin','MaxFunEval',1E9,'MaxIter',1E12,'TolFun',...
        1e-10,'Display','none'); %'Algorithm','sqp'
    thetaSolver2D = @(theta) ([cos(theta) sin(theta);-sin(theta) cos(theta)] - R);
    [theta,resnorm] = lsqnonlin(thetaSolver2D,theta,-2*pi, 2*pi,options2);%theta);
    disp '***********'
    disp '>>End Solver details'
    disp ' '

    c=cos(theta); s=sin(theta);
    R_sol = [c s;...
            -s c];

    disp(['Euler angle /degree: (gamma) = (' num2str(radtodeg(theta),'% 1.6f ') ').'])
    %OC2
    %% 4.3 Generate theoretical displacement field
    xip_theo=(R_sol*xi0')';
    disp_theo=xi0-xip_theo;       %= orig coords - rotn coords

    %% 4.4 Subtract theoretical displacements
    % Ux,Uy and Uz are displacements. Parameters without a prefix are
    % original and rb means rotated back. big_new_data_natural is the final
    % matrix similar to data corrected for rotation referenced to input x-y
    % frame.

    % rUi = disp_theo(:,1:2)  {theoretical displacments} 
    % Ui = cdata(:,3:4)       {measured displacements}
    % rbUi = Ui - rUi         {deformation displacements}

    rbUi = cdata(:,3:4) - disp_theo(:,1:2); 
    new_data=[cdata(:,1:2) rbUi];      %[CleanedCoords DeformDisps] N.B. cdata has been recentred c.f. tdata
    
    %% Outputs prints
    disp(' ')
    disp(['RBD /pixel = ' num2str(RBD,6)])
    disp(['rotCentre /pixel = ' num2str(rotCentre,6) ])
    disp(['Euler angle /degree: (gamma) = (' num2str(radtodeg(theta),'% 1.6f ') ').'])
    
else
    new_data = cdata;
    theta    = nan;
    rotCentre = nan(1,2);
    disp ' '
%     disp(['RBD /pixel = ' num2str(RBD,6)])
    disp('No correction for rotation requested');
end

%% 5. Recreate original and cleaned displacement fields
big_new_data_natural=data2D_col;

j=1;
for i=1:size(big_new_data_natural,1)
    if ~isnan(big_new_data_natural(i,4))
        big_new_data_natural(i,:)=new_data(j,:);
        j=j+1;    
    end
end
%% 6. Define outputs
% eulerAngles = theta;
mesh.dDIC = big_new_data_natural(:,3:4)';
% Correctly fill 'el' structure
for i=1:length(mesh.UDIC)
    pos = find(el.Ux(:)==mesh.UDIC(1,i) & el.Uy(:)==mesh.UDIC(2,i));
    if(pos)
        el.dx(pos) = mesh.dDIC(1,i);
        el.dy(pos) = mesh.dDIC(2,i);
    end
end

disp(' ')
%OC3
end
