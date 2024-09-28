function [Fn, S] = Euler_HLLC_flux_primitive(U_l, U_r, normal, gamma)
% [Fn, S] = Euler_HLLC_flux_primitive(U_l, U_r, normal, gamma)
%   compute the HLLC numerical flux on [nx, ny, nz] direction (for ideal gas only) based on primitive variables
%   can adapt to problem dimension automatically
% 
% recommended: 
%   HLL-family fluxes (positivity-preserving), 
%   the classical HLL flux (entropy-stable if the signal speed is estimated properly and 1<gamma<=5/3), 
%   the HLLC flux (positivity-preserving) 
% 
% input:
%   primitive variables (U_l, U_r): [density, velocity on x,y,z-directions, pressure]
%   unit normal vector
%   specific heat ratio
% 
% output:
%   Fn: numerical flux 
%   S:  estimate of the maximal wave speed
% 
% Refs: 
%   [1] On the Choice of Wavespeeds for the HLLC Riemann Solver
%   [2] On positivity-preserving high order discontinuous Galerkin schemes for compressible Euler equations on rectangular meshes
%   [3] Entropy stable high order discontinuous Galerkin methods with suitable quarature rules for hyperbolic conservation laws

d = numel(U_l)-2; % problem dimension
U_l = U_l(:);
U_r = U_r(:);
normal = normal(:);
normal = normal(1:d);

rho_l = U_l(1);
v_l = U_l(2:(d+1));
p_l = U_l(end);
E_l = p_l/(gamma-1.0) + 0.5*rho_l*sum(v_l(:).^2);
c_l = sqrt(gamma*p_l/rho_l); % sonic speed
H_l = (E_l + p_l) / rho_l; % enthalpy
vn_l = dot(v_l(:), normal(:)); % normal speed

rho_r = U_r(1);
v_r = U_r(2:(d+1));
p_r = U_r(end);
E_r = p_r/(gamma-1.0) + 0.5*rho_r*sum(v_r(:).^2);
c_r = sqrt(gamma*p_r/rho_r); % sonic speed
H_r = (E_r + p_r) / rho_r; % enthalpy
vn_r = dot(v_r(:), normal(:)); % normal speed

% compute the Roe averaged state that is ONLY used for estimating the signal speeds
R_roe = sqrt(rho_r/rho_l);
deno = 1.0 + R_roe;
v_roe = (v_l(:) + R_roe*v_r(:)) / deno;
H_roe = (H_l + R_roe*H_r) / deno;
vn_roe = (vn_l + R_roe*vn_r) / deno;
c_roe = sqrt((gamma-1.0)*(H_roe - 0.5*sum(v_roe(:).^2)));

% Roe's estimated left- and right- signal speeds
S_l = min(vn_l - c_l, vn_roe - c_roe);
S_r = max(vn_r + c_r, vn_roe + c_roe);

% compute the HLLC numerical flux
Fn = nan(d+2,1);
S_m = 0.0;
if S_l > 0.0
    % supersonic flow toward right
    Fn(1) = vn_l * rho_l;
    Fn(2:(d+1)) = Fn(1) * v_l(:) + p_l * normal(:);
    Fn(end) = vn_l * (E_l + p_l);
elseif S_r < 0.0
    % supersonic flow toward left
    Fn(1) = vn_r * rho_r;
    Fn(2:(d+1)) = Fn(1) * v_r(:) + p_r * normal(:);
    Fn(end) = vn_r * (E_r + p_r);
else
    % subsonic, data from both left and right
    flow_l = rho_l * (vn_l - S_l);
    flow_r = rho_r * (vn_r - S_r);

    % compute the estimated contact wave speed
    S_m = ((flow_l*vn_l - flow_r*vn_r) + (p_l - p_r)) / (flow_l - flow_r);

    if S_m >= 0.0
        % to the left of the contact wave
        % Fn = F_l + S_l*(U_l^* - U_l)
        
        % checked by Mathematica
        ts_l = ((vn_l - S_m) / (S_l - S_m)) * S_l;
        tsflow_l = ts_l * flow_l;
        vv_l = vn_l - ts_l;
        Fn(1) = vv_l * rho_l;
        Fn(2:(d+1)) = Fn(1) * v_l(:) + (p_l + tsflow_l) * normal(:);
        Fn(end) = vv_l * (E_l + p_l) + tsflow_l * S_m;
    else
        % to the right of the contact wave
        % Fn = F_r + S_r*(U_r^* - U_r)

        % checked by Mathematica
        ts_r = ((vn_r - S_m) / (S_r - S_m)) * S_r;
        tsflow_r = ts_r * flow_r;
        vv_r = vn_r - ts_r;
        Fn(1) = vv_r * rho_r;
        Fn(2:(d+1)) = Fn(1) * v_r(:) + (p_r + tsflow_r) * normal(:);
        Fn(end) = vv_r * (E_r + p_r) + tsflow_r * S_m;
    end
end

S = max(abs([S_l, S_m, S_r]));

end