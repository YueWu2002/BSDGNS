function [ut, max_dt] = DG_1D_semidiscrete(mesh, u, p, qw, Minv, D, Tl, Tr, gamma)

ut = zeros(p+1, 5, mesh.Nx);
max_dtx = Inf(1, mesh.Nx); % maximum step size for 1st order FV scheme for each element

% intra-element update
for k = 1: mesh.Nx
    hx = mesh.nodes(mesh.elems(2,k)) - mesh.nodes(mesh.elems(1,k));
    Jx = [hx/2.0];
    
    [pressure, sonic] = Euler_EOS(u(:,1,k), u(:,2,k), u(:,3,k), u(:,4,k), u(:,5,k), gamma);
    if (any(pressure(:) < 0.0) || any(u(:,1,k) <= 0.0))
        error('Negative pressure or density! Computation aborted!');
    end
    S = sonic(:) + abs(u(:,2,k)./u(:,1,k)); % max eigen values in x-direction
    Fx = nan(p+1, 5); % advective flux in x-direction
    for j = 1: p+1
        [Fx(j,:), ~, ~] = Euler_advective_flux(u(j,:,k), pressure(j));
    end
    ut(:,:,k) = ut(:,:,k) + (1.0./det(Jx)) .* (Minv * (det(Jx) .* inv(Jx)' .* (D' * (qw(:) .* Fx(:,:)))));

    max_local_eigx = max(S(:));
    max_dtx(k) = min(max_dtx(k), hx/max_local_eigx);
end

% inter-element update
for k = 1: mesh.Nx-1
    idl = k; % left
    idr = k+1; % right

    hxl = mesh.nodes(mesh.elems(2,idl)) - mesh.nodes(mesh.elems(1,idl));
    Jxl = [hxl/2.0];

    hxr = mesh.nodes(mesh.elems(2,idr)) - mesh.nodes(mesh.elems(1,idr));
    Jxr = [hxr/2.0];

    LEFTc = Tr*u(:,:,idl);
    RIGHTc = Tl*u(:,:,idr);
    % [Fn, S] = Euler_HLLC_flux_primitive(Euler_con2pri(LEFTc, gamma), Euler_con2pri(RIGHTc, gamma), [1.0, 0.0, 0.0], gamma);
    [Fn, S] = Euler_LF_flux_primitive(Euler_con2pri(LEFTc, gamma), Euler_con2pri(RIGHTc, gamma), [1.0, 0.0, 0.0], gamma);

    max_dtx(idl) = min(max_dtx(idl), hxl/S);
    max_dtx(idr) = min(max_dtx(idr), hxr/S);

    ut(:,:,idl) = ut(:,:,idl) + (1.0./det(Jxl)) .* (Minv * (-Tr'*1*Fn(:)'));
    ut(:,:,idr) = ut(:,:,idr) + (1.0./det(Jxr)) .* (Minv * (+Tl'*1*Fn(:)'));
end

% boundary update

%{
% periodic BC

idl = mesh.Nx; % left
idr = 1; % right

hxl = mesh.nodes(mesh.elems(2,idl)) - mesh.nodes(mesh.elems(1,idl));
Jxl = [hxl/2.0];

hxr = mesh.nodes(mesh.elems(2,idr)) - mesh.nodes(mesh.elems(1,idr));
Jxr = [hxr/2.0];

LEFTc = Tr*u(:,:,idl);
RIGHTc = Tl*u(:,:,idr);
[Fn, S] = Euler_HLLC_flux_primitive(Euler_con2pri(LEFTc, gamma), Euler_con2pri(RIGHTc, gamma), [1.0, 0.0, 0.0], gamma);

max_dtx(idl) = min(max_dtx(idl), hxl/S);
max_dtx(idr) = min(max_dtx(idr), hxr/S);

ut(:,:,idl) = ut(:,:,idl) + (1.0./det(Jxl)) .* (Minv * (-Tr'*1*Fn(:)'));
ut(:,:,idr) = ut(:,:,idr) + (1.0./det(Jxr)) .* (Minv * (+Tl'*1*Fn(:)'));
%}


% shock tube problem (extrapolation boundary)
idr = 1;
hxr = mesh.nodes(mesh.elems(2,idr)) - mesh.nodes(mesh.elems(1,idr));
Jxr = [hxr/2.0];

LEFTc = Tl*u(:,:,idr); % extrapolation 
RIGHTc = Tl*u(:,:,idr);
[Fn, S] = Euler_HLLC_flux_primitive(Euler_con2pri(LEFTc, gamma), Euler_con2pri(RIGHTc, gamma), [1.0, 0.0, 0.0], gamma);
max_dtx(idr) = min(max_dtx(idr), hxr/S);
ut(:,:,idr) = ut(:,:,idr) + (1.0./det(Jxr)) .* (Minv * (+Tl'*1*Fn(:)'));

idl = mesh.Nx;
hxl = mesh.nodes(mesh.elems(2,idl)) - mesh.nodes(mesh.elems(1,idl));
Jxl = [hxl/2.0];

LEFTc = Tr*u(:,:,idl);
RIGHTc = Tr*u(:,:,idl); % extrapolation
[Fn, S] = Euler_HLLC_flux_primitive(Euler_con2pri(LEFTc, gamma), Euler_con2pri(RIGHTc, gamma), [1.0, 0.0, 0.0], gamma);
max_dtx(idl) = min(max_dtx(idl), hxl/S);
ut(:,:,idl) = ut(:,:,idl) + (1.0./det(Jxl)) .* (Minv * (-Tr'*1*Fn(:)'));



max_dt = min(max_dtx(:));

end