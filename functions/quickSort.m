function sortedArray = quickSort(array)
 
    if numel(array) <= 1 %If the array has 1 element  can't be sorted       
        sortedArray = array;
        return
    end
 
    pivot = array(:,end);
    array(:,end) = [];

    maskle = array(2,:) <= pivot(2,:);
    k=find(maskle~=0);
    less = array(:,k);
    
    maskgr = array(2,:) > pivot(2,:);
    k=find(maskgr~=0);
    greater = array(:,k);

    sortedArray = [quickSort(less) pivot quickSort(greater)];
 
end