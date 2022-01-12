function dat = unifromMesh(datum)
StepSize = [diff(unique(datum.X1));diff(unique(datum.Y1))];
StepSize = min(round(StepSize,4));
    xLin  = min(round(datum.X1(:),4)):StepSize:max(round(datum.X1(:),4));
    yLin  = min(round(datum.Y1(:),4)):StepSize:max(round(datum.Y1(:),4));
        [X,Y] = meshgrid(xLin,yLin); 
    X = X(:);       Y = Y(:);
    
	F11 = scatteredInterpolant(datum.X1(:),datum.Y1(:),datum.Ux(:),'natural');
    F22 = scatteredInterpolant(datum.X1(:),datum.Y1(:),datum.Uy(:),'natural');
    newM(:,1) = X(:);       
    newM(:,2) = Y(:);
    newM(:,3) = F11(X,Y); %Evaluate FE data on grid of experimental results
	newM(:,4) = F22(X,Y); %Evaluate FE data on grid of experimental results

[~,dat] = reshapeData(newM); 
end