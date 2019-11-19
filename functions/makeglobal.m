function [gl] = makeglobal(el, mesh)
if mesh.Operation=='Dis'
    %find non zero displacment nodes
    logic = ((sqrt(mesh.dDIC(1,:).^2+mesh.dDIC(2,:).^2))~=0);
    logic = ones(size(logic));

    %global U
    ux = reshape(logic.*mesh.UDIC(1,:), mesh.winDIC(2),mesh.winDIC(1));
    uy = reshape(logic.*mesh.UDIC(2,:),mesh.winDIC(2),mesh.winDIC(1));
    gl.Ux = ux(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));
    gl.Uy = uy(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));

    %global d
    dx = reshape(logic.*mesh.dDIC(1,:),mesh.winDIC(2),mesh.winDIC(1));
    dy = reshape(logic.*mesh.dDIC(2,:), mesh.winDIC(2),mesh.winDIC(1));
    gl.dx = dx(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));
    gl.dy = dy(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));

    %global deformed
    gl.Uxdef = gl.Ux;%+gl.dx*1;
    gl.Uydef = gl.Uy;%+gl.dy*1;

    %interpolate VM stress at nodes
    gl.smises = griddata(el.gpx,el.gpy,el.gpsmises,gl.Ux, gl.Uy,'linear');
%     gl.smises(gl.smises>3*nanmean(gl.smises(:))) = NaN;

elseif mesh.Operation=='Str'
%find non zero displacment nodes
logic = ((sqrt(mesh.dDIC(1,:).^2+mesh.dDIC(2,:).^2))~=0);

%global U
ux = reshape(logic.*mesh.UDIC(1,:), mesh.winDIC(2),mesh.winDIC(1));
uy = reshape(logic.*mesh.UDIC(2,:),mesh.winDIC(2),mesh.winDIC(1));
gl.Ux = ux(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));
gl.Uy = uy(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));

%global d
dx = reshape(logic.*mesh.dDIC(1,:),mesh.winDIC(2),mesh.winDIC(1));
dy = reshape(logic.*mesh.dDIC(2,:), mesh.winDIC(2),mesh.winDIC(1));
gl.dx = dx(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));
gl.dy = dy(mesh.winFE(1,2):mesh.winFE(2,2),mesh.winFE(1,1):mesh.winFE(2,1));

%global deformed
gl.Uxdef = gl.Ux;       %+gl.dx*1;
gl.Uydef = gl.Uy;       %+gl.dy*1;

%interpolate VM stress at nodes
gl.smises = griddata(el.gpx,el.gpy,el.gpsmises,gl.Ux, gl.Uy,'linear');
gl.SMISES = gl.smises;
gl.smises(gl.smises>3*nanmean(gl.smises(:))) = NaN;

gl.exx = griddata(el.gpx,el.gpy,el.gpexx,gl.Ux, gl.Uy,'linear');
gl.eyy = griddata(el.gpx,el.gpy,el.gpeyy,gl.Ux, gl.Uy,'linear');
gl.exy = griddata(el.gpx,el.gpy,el.gpexy,gl.Ux, gl.Uy,'linear');
end