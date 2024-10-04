function [Fn, S] = Euler_HLLC_flux_primitive(U_l, U_r, normal, gamma, dim)
% [Fn, S] = Euler_HLLC_flux_primitive(U_l, U_r, normal, gamma, dim)
%   compute the HLLC numerical flux on normal direction (for ideal gas only) using primitive variables
%   (support any problem dimension)
% 
% recommended: 
%   HLL-family fluxes (positivity-preserving, resolve 1- and 3- discontinuity exactly), 
%   the classical HLL flux (entropy-stable if the signal speed is estimated properly and 1<gamma<=5/3), 
%   the HLLC flux (positivity-preserving, can resolve 1-, 2- and 3- discontinuity exactly) 
% 
% input:
%   U_l, U_r:   primitive variables, numel == dim+2
%   normal:     unit normal vector, numel == dim
%   gamma:      the specific heat ratio in the gamma law EOS, scalar
%   dim:        problem dimension, integer scalar
%   
% output:
%   Fn:         numerical flux, size = [dim+2, 1]
%   S:          estimate of the maximal wave speed, scalar
% 
% Refs: 
%   [1] On the Choice of Wavespeeds for the HLLC Riemann Solver
%   [2] On positivity-preserving high order discontinuous Galerkin schemes for compressible Euler equations on rectangular meshes
%   [3] Entropy stable high order discontinuous Galerkin methods with suitable quarature rules for hyperbolic conservation laws

normal = reshape(normal(1:dim),[dim,1]);

rho_l = U_l(1);
v_l = reshape(U_l(2:dim+1), [dim,1]);
p_l = U_l(dim+2);
c_l_sq = gamma*p_l/rho_l;
c_l = sqrt(c_l_sq); % sonic speed
E_l = p_l/(gamma-1.0) + 0.5*rho_l*sum(v_l(:).^2);
vn_l = dot(v_l(:), normal(1:dim)); % normal speed

rho_r = U_r(1);
v_r = reshape(U_r(2:dim+1), [dim,1]);
p_r = U_r(dim+2);
c_r_sq = gamma*p_r/rho_r;
c_r = sqrt(c_r_sq); % sonic speed
E_r = p_r/(gamma-1.0) + 0.5*rho_r*sum(v_r(:).^2);
vn_r = dot(v_r(:), normal(1:dim)); % normal speed

% compute the Roe averaged state that is ONLY used for estimating the signal speeds
R_roe = sqrt(rho_r/rho_l);
deno = 1.0 + R_roe;
vn_roe = (vn_l + R_roe*vn_r) / deno;
% verified by Mathematica
c_roe_sq = (0.5*(gamma-1.0)/(R_roe + 2.0 + 1.0/R_roe)) * sum((v_l(:) - v_r(:)).^2) + (c_l_sq + R_roe * c_r_sq) / deno;
c_roe = sqrt(c_roe_sq);

% Roe's estimated left- and right- signal speeds
S_l = min(vn_l - c_l, vn_roe - c_roe);
S_r = max(vn_r + c_r, vn_roe + c_roe);

% compute the HLLC numerical flux
Fn = nan(dim+2,1);
S_m = 0.0;
if S_l > 0.0
    % supersonic flow toward right
    Fn(1) = vn_l * rho_l;
    Fn(2:dim+1) = Fn(1) * v_l(:) + p_l * normal(1:dim);
    Fn(dim+2) = vn_l * (E_l + p_l);
elseif S_r < 0.0
    % supersonic flow toward left
    Fn(1) = vn_r * rho_r;
    Fn(2:dim+1) = Fn(1) * v_r(:) + p_r * normal(1:dim);
    Fn(dim+2) = vn_r * (E_r + p_r);
else
    % subsonic, data from both left and right
    flow_l = rho_l * (vn_l - S_l);
    flow_r = rho_r * (vn_r - S_r);

    % compute the estimated contact wave speed
    S_m = ((flow_l*vn_l - flow_r*vn_r) + (p_l - p_r)) / (flow_l - flow_r);

    if S_m >= 0.0
        % to the left of the contact wave
        % Fn = F_l + S_l*(U_cl^* - U_cl)
        
        % checked by Mathematica
        ts_l = ((vn_l - S_m) / (S_l - S_m)) * S_l;
        tsflow_l = ts_l * flow_l;
        vv_l = vn_l - ts_l;
        Fn(1) = vv_l * rho_l;
        Fn(2:dim+1) = Fn(1) * v_l(:) + (p_l + tsflow_l) * normal(1:dim);
        Fn(dim+2) = vv_l * (E_l + p_l) + tsflow_l * S_m;
    else
        % to the right of the contact wave
        % Fn = F_r + S_r*(U_cr^* - U_cr)

        % checked by Mathematica
        ts_r = ((vn_r - S_m) / (S_r - S_m)) * S_r;
        tsflow_r = ts_r * flow_r;
        vv_r = vn_r - ts_r;
        Fn(1) = vv_r * rho_r;
        Fn(2:dim+1) = Fn(1) * v_r(:) + (p_r + tsflow_r) * normal(1:dim);
        Fn(dim+2) = vv_r * (E_r + p_r) + tsflow_r * S_m;
    end
end

S = max(abs([S_l, S_m, S_r]));

end