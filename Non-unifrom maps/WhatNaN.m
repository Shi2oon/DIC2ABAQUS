function [Elements,E11,E22,E12,X,Y,answer] = WhatNaN(Elements,E11,E22,E12,X,Y)
[f(:,1),f(:,2)]=ind2sub(size(E11(1,:)),find(isnan(E11(1,:))));
if length(f)/length(E11) == 0
    clearvars f
    [f(:,1),f(:,2)]=ind2sub(size(E11(1,:)),find(E11(1,:)==0));
end
if length(f)/length(E11) > 0.01
[LV1] = RemoveOutNaN(E11); %CHECK FOR OUTLIERS
[f1(:,1),f1(:,2)] = ind2sub(size(LV1),find(LV1==0)); % check how many we found
    % check if the opercentage is acceptable
	if length(f1)/length(LV1) > 0.15 && length(f1)/length(LV1) < 0.85 ...
            &&  abs(length(f)/length(E11)-length(f1)/length(LV1)) < 0.15
        % do nothing and use calculated outliner
    else
        LV1 = ones(size(LV1));
        LV1(f) = 0;
	end
	fprintf('\nOutliers percentage is %1.2f (> 0.5 allowed). I will now trim them out\n',...
        length(f)/length(E11)*100);
	E11(:,LV1==0) = []; 		E22(:,LV1==0) = []; 
	E12(:,LV1==0) = []; 		X(:,LV1==0) = []; 
	Y(:,LV1==0) = [];			Elements = Elements';   
	Elements(:,LV1==0) = [];    Elements = Elements';  
    E11(isnan(E11)) = 0;        E22(isnan(E22)) = 0;        E12(isnan(E12)) = 0;
	answer = 'Y';
else
    answer = 'N';
    E11(isnan(E11)) = 0;        E22(isnan(E22)) = 0;        E12(isnan(E12)) = 0;
end
% fill(X,Y,'w'); axis image; fill(X,Y,E12,'LineStyle','none'); axis image; shg
end

function [LV1] = RemoveOutNaN(tmp)
tmp(tmp==0) = NaN;    
mdn = nanmedian(tmp(:));
tmp(isnan(tmp)) = mdn;
LV1 = isoutlier(tmp,'median',2); 
LV1 = max(LV1);
end