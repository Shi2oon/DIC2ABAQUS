function [Size] = GetSize(Aray)

nb=Aray(1);
count = 1;

for i = 2 : length(Aray)
    if(Aray(i) == nb)
        count = count + 1; 
%     else
%         break;
    end
end

 Size(1)=count;
 Size(2)=length(Aray)/count;