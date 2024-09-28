function [U_roe] = Euler_Roe_average(U_l, U_r, gamma)
% compute the Roe's average of primitive variables
% can adapt to dimension

d = numel(U_l) - 2;

rho_l = U_l(1);
v_l = U_l(2:(d+1));
p_l = U_l(end);
E_l = p_l/(gamma-1.0) + 0.5*rho_l*sum(v_l(:).^2);
H_l = (E_l + p_l) / rho_l; % enthalpy

rho_r = U_r(1);
v_r = U_r(2:(d+1));
p_r = U_r(end);
E_r = p_r/(gamma-1.0) + 0.5*rho_r*sum(v_r(:).^2);
H_r = (E_r + p_r) / rho_r; % enthalpy

R_roe = sqrt(rho_r/rho_l);
rho_roe = sqrt(rho_l*rho_r);
deno = 1.0 + R_roe;
v_roe = (v_l(:) + R_roe*v_r(:)) / deno;
H_roe = (H_l + R_roe*H_r) / deno;
p_roe = ((gamma-1.0)/gamma)*rho_roe*(H_roe - 0.5*sum(v_roe(:).^2));

U_roe = [rho_roe; v_roe; p_roe];

end