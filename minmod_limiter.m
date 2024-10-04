function [u] = minmod_limiter(mesh, u, p, M, Tl, Tr)
if p == 0
    return;
end

% compute the cell averages
u_avg = (0.5*ones(1,p+1)*M)*u(:,:);

% derive the difference between j and j-1
u_avg_diff = nan(1, mesh.Nx+1);
u_avg_diff(2:mesh.Nx) = u_avg(2:mesh.Nx) - u_avg(1:mesh.Nx-1);
% use the periodic boundary conditions
u_avg_diff(1) = u_avg(1) - u_avg(mesh.Nx);
u_avg_diff(mesh.Nx+1) = u_avg(1) - u_avg(mesh.Nx);

% compute the right lifting and left liftings (they should be the same when p==1)
u_rr = Tr*u(:,:) - u_avg;
u_ll = u_avg - Tl*u(:,:);

u_rr_new = minmod_3(u_rr(1:mesh.Nx), u_avg_diff(1:mesh.Nx), u_avg_diff(2:mesh.Nx+1));
u_ll_new = minmod_3(u_ll(1:mesh.Nx), u_avg_diff(1:mesh.Nx), u_avg_diff(2:mesh.Nx+1));

if p == 1
    % use uniqueness
    u = [Tl; Tr] \ [u_avg - u_ll_new(:)'; u_avg + u_rr_new(:)'];
    return;
end

% use L2 best approximation
B = [Tl; 0.5*ones(1,p+1)*M; Tr];
sol = full([M, B'; B, zeros(3,3)]) \ [sparse(p+1,mesh.Nx); u_ll(:)' - u_ll_new(:)'; sparse(1,mesh.Nx); u_rr_new(:)' - u_rr(:)'];
u(:,:) = u(:,:) + sol(1:p+1,:);

end