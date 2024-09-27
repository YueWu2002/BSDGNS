function [Fn, S] = Euler_HLLC_flux(U_l, U_r, normal, gamma)
% [Fn, S] = Euler_HLLC_flux(U_l, U_r, normal, gamma)
%   compute the HLLC numerical flux on [nx, ny, nz] direction (for ideal gas only)
% 
% recommended: 
%   HLL-family fluxes (positivity-preserving), 
%   the classical HLL flux (entropy-stable if the signal speed is estimated properly and 1<gamma<=5/3), 
%   the HLLC flux (positivity-preserving) 
% 
% input:
%   conservative variables
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

% transform to primitive variables and compute related physical quantities
rho_l = U_l(1);
v_l = U_l(2:4)/U_l(1);
E_l = U_l(5);
p_l = (gamma-1.0)*(E_l - 0.5*rho_l*sum(v_l(:).^2)); % pressure
c_l = sqrt(gamma*p_l/rho_l); % sonic speed
H_l = (E_l + p_l) / rho_l; % enthalpy
vn_l = dot(v_l(:), normal(:)); % normal speed

% transform to primitive variables and compute related physical quantities
rho_r = U_r(1);
v_r = U_r(2:4)/U_r(1);
E_r = U_r(5);
p_r = (gamma-1.0)*(E_r - 0.5*rho_r*sum(v_r(:).^2)); % pressure
c_r = sqrt(gamma*p_r/rho_r); % sonic speed
H_r = (E_r + p_r) / rho_r;
vn_r = dot(v_r(:), normal(:)); % normal speed

% compute the Roe averaged state that is ONLY used for estimating the signal speeds
R_roe = sqrt(rho_r/rho_l);
v_roe = (v_l(:) + R_roe*v_r(:)) / (1.0 + R_roe);
H_roe = (H_l + R_roe*H_r) / (1.0 + R_roe);
vn_roe = (vn_l + R_roe*vn_r) / (1.0 + R_roe);
c_roe = sqrt((gamma-1.0)*(H_roe - 0.5*sum(v_roe(:).^2)));

% Roe's estimated left- and right- signal speeds
S_l = min(vn_l - c_l, vn_roe - c_roe);
S_r = max(vn_r + c_r, vn_roe + c_roe);

% compute the HLLC numerical flux
% Fn = nan(5,1);
if S_l > 0.0
    % supersonic flow toward right
    Fn = vn_l * U_l(:) + p_l * [0.0; normal(:); vn_l];
elseif S_r < 0.0
    % supersonic flow toward left
    Fn = vn_r * U_r(:) + p_r * [0.0; normal(:); vn_r];
else
    % subsonic, data from both left and right
    flow_l = rho_l * (vn_l - S_l);
    flow_r = rho_r * (vn_r - S_r);
    S_m = ((flow_l*vn_l - flow_r*vn_r) + (p_l - p_r)) / (flow_l - flow_r);

    if S_m >= 0.0
        % to the left of the contact wave
        % Fn = F_l + S_l*(U_l^* - U_l)

        % rho_lc = flow_l / (S_m - S_l);
        % p_c = p_l + flow_l * (vn_l - S_m); % same value for left and right
        % v_lc = v_l(:) - (vn_l - S_m)*normal(:); % flow speed in the middle-left region

        % Fn(1) = rho_lc*S_m;
        % Fn(1) = (S_m/(S_m - S_l)) * flow_l;
        % Fn(2:4) = Fn(1) * v_lc(:) + p_c * normal(:);
        % Fn(2:4) = Fn(1) * v_l(:) + (p_l - flow_l * (vn_l - S_m) * (S_l/(S_m - S_l))) * normal(:);
        % Fn(5) = ((S_l - vn_l)/(S_l - S_m)) * (E_l + p_l - rho_l*S_l*(vn_l - S_m)) * S_m;

        t_l = (vn_l - S_m) / (S_l - S_m); % t_l <= 1
        % Fn = (vn_l * U_l(:) + p_l * [0.0; normal(:); vn_l]) + (S_l * t_l) * ([0.0; flow_l*normal(:); flow_l*S_m - p_l] - U_l(:));
        % checked by Mathematica

        ts_l = t_l*S_l;
        vv_l = vn_l - ts_l;
        Fn = vv_l * U_l(:) + p_l * [0.0; normal(:); vv_l] + (ts_l * flow_l) * [0.0; normal(:); S_m];
    else
        % to the right of the contact wave
        % Fn = F_r + S_r*(U_r^* - U_r)

        % rho_rc = flow_r / (S_m - S_r);
        % p_c = p_r + flow_r * (vn_r - S_m); % same value for left and right
        % v_rc = v_r(:) - (vn_r - S_m) * normal(:); % flow speed

        % Fn(1) = rho_rc*S_m;
        % Fn(1) = (S_m/(S_m - S_r)) * flow_r;
        % Fn(2:4) = Fn(1) * v_rc(:) + p_c * normal(:);
        % Fn(2:4) = Fn(1) * v_r(:) + (p_r - flow_r * (vn_r - S_m) * (S_r/(S_m - S_r))) * normal(:);
        % Fn(5) = ((S_r - vn_r)/(S_r - S_m)) * (E_r + p_r - rho_r*S_r*(vn_r - S_m)) * S_m;

        t_r = (vn_r - S_m) / (S_r - S_m); % t_r <= 1
        % Fn = (vn_r * U_r(:) + p_r * [0.0; normal(:); vn_r]) + (S_r * t_r) * ([0.0; flow_r*normal(:); flow_r*S_m - p_r] - U_r(:));
        % checked by Mathematica

        ts_r = t_r*S_r;
        vv_r = vn_r - ts_r;
        Fn = vv_r * U_r(:) + p_r * [0.0; normal(:); vv_r] + (ts_r * flow_r) * [0.0; normal(:); S_m];
    end
end

S = max(abs([S_l, S_m, S_r]));

end