function datum = unifromMesh(datum)
StepSize = [diff(unique(datum.X));diff(unique(datum.Y))];
StepSize = min(round(StepSize,4));
    xLin  = min(round(datum.X(:),4)):StepSize:max(round(datum.X(:),4));
    yLin  = min(round(datum.Y(:),4)):StepSize:max(round(datum.Y(:),4));
        [X,Y] = meshgrid(xLin,yLin); 
    X = X(:);       Y = Y(:);
    
	F11 = scatteredInterpolant(datum.X(:),datum.Y(:),datum.Ux(:),'natural');
    F22 = scatteredInterpolant(datum.X(:),datum.Y(:),datum.Uy(:),'natural');
    newM(:,1) = X(:);       
    newM(:,2) = Y(:);
    newM(:,3) = F11(X,Y); %Evaluate FE data on grid of experimental results
	newM(:,4) = F22(X,Y); %Evaluate FE data on grid of experimental results

[~,dat] = reshapeData(newM); 
datum.X = dat.X1;       datum.Y = dat.Y1;
datum.Ux = dat.Ux;      datum.Uy = dat.Uy;
end