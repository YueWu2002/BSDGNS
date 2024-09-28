function [Fn, S] = Euler_LF_flux_primitive(U_l, U_r, normal, gamma)
% conservative Lax-Friedrichs flux from primitive variables

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
vn_l = dot(v_l(:), normal(:)); % normal speed

rho_r = U_r(1);
v_r = U_r(2:(d+1));
p_r = U_r(end);
E_r = p_r/(gamma-1.0) + 0.5*rho_r*sum(v_r(:).^2);
c_r = sqrt(gamma*p_r/rho_r); % sonic speed
vn_r = dot(v_r(:), normal(:)); % normal speed

S = max(abs(vn_l) + c_l, abs(vn_r) + c_r);

Fl = nan(d+2,1);
Fr = nan(d+2,1);

Fl(1) = vn_l * rho_l;
Fl(2:(d+1)) = Fl(1) * v_l(:) + p_l * normal(:);
Fl(end) = vn_l * (E_l + p_l);

Fr(1) = vn_r * rho_r;
Fr(2:(d+1)) = Fr(1) * v_r(:) + p_r * normal(:);
Fr(end) = vn_r * (E_r + p_r);

Fn = 0.5*(Fl(:) + Fr(:)) - 0.5*S*[rho_r - rho_l; rho_r*v_r(:) - rho_l*v_l(:); E_r - E_l];

end