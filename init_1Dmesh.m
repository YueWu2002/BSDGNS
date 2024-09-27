clear mesh;
mesh.a = -5.0;
mesh.b = 5.0;
mesh.Nx = 200;
mesh.nodes = linspace(mesh.a, mesh.b, mesh.Nx+1); % axis of global nodes
mesh.elems = zeros(2, mesh.Nx, 'int32'); % ids of global nodes of each element
mesh.elemnode = nan(1, 2, mesh.Nx); % axis of nodes of each element
mesh.xx = nan(p+1, mesh.Nx); % collocation points of each element (for plotting and initializing purposes only, not used in computing)
for k = 1: mesh.Nx
    mesh.elems(:,k) = [k, k+1];

    mesh.elemnode(:,1,k) = mesh.nodes(k);
    mesh.elemnode(:,2,k) = mesh.nodes(k+1);
    
    mesh.xx(:,k) = 0.5*(mesh.nodes(k) + mesh.nodes(k+1)) + 0.5*(mesh.nodes(k+1) - mesh.nodes(k)).*qx(:);
end