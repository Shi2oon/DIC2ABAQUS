function [indice] = findClosest(target,array)
% find the closest X and Y couple to the tqrget value in a given array.

diff_array = abs(array-repmat(target,length(array),1));
diff_array = sum(diff_array,2);

indice = find(diff_array == min(diff_array));