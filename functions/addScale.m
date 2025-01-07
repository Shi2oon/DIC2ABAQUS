function addScale(No,alldata)
% funciton to add a measurment line to the graph, you can either input
% number of the suplots or the exact suplot where you want to add the scale
% bar
    if length(No) == 1      && No == 2
        subplot(1,2,1); AddScalePar(alldata(:,1),alldata(:,2))
        subplot(1,2,2); AddScalePar(alldata(:,1),alldata(:,2))
    elseif length(No) == 1  && No == 3
        subplot(1,3,1); AddScalePar(alldata(:,1),alldata(:,2))
        subplot(1,3,2); AddScalePar(alldata(:,1),alldata(:,2))
        subplot(1,3,3); AddScalePar(alldata(:,1),alldata(:,2))
    elseif length(No) == 1	&& No == 1
        AddScalePar(alldata(:,1),alldata(:,2));
    elseif length(No) == 3
        subplot(No(1),No(2),No(3)); AddScalePar(alldata(:,1),alldata(:,2));
    end
end

function AddScalePar(X,Y)
X = unique(X);         Y = unique(Y);
hold on; line([X(ceil(end-length(X)*0.05)),X(end-ceil(length(X)*0.15))],...
            [Y(ceil(length(Y)*0.05)),Y(ceil(length(Y)*0.05))],'Color','k','LineWidth',5)
ht = text(double(X(end-ceil(length(X)*0.24))),double(Y(ceil(length(Y)*0.14)))...
    ,[num2str(round(abs(X(ceil(length(X)*0.05))-X(ceil(length(X)*0.15))),1)) '\mum']);
set(ht,'FontSize',16,'FontWeight','Bold')
set(gca,'Visible','off')
hold off
end