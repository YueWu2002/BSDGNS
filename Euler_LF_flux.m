function [Fn, S] = Euler_LF_flux(U_l, U_r, normal, gamma)
% [Fn, S] = Euler_LF_flux(U_l, U_r, normal, gamma)
%   compute the local Lax-Friedrichs numerical flux, a.k.a. the Rusanov flux
%   (for ideal gas only)
% 
% input:
%   conservative variables
%   unit normal vector
% 
% output:
%   Fn: numerical flux 
%   S:  estimate of the maximal wave speed

Fnl = nan(5,1);
Fnr = nan(5,1);

U_l = U_l(:);
U_r = U_r(:);

p_l = (gamma-1.0)*(U_l(5) - 0.5*sum(U_l(2:4).^2)/U_l(1));
Fnl(1) = dot(normal(:), U_l(2:4));
vn_l = Fnl(1) / U_l(1);
Fnl(2:4) = vn_l*U_l(2:4) + p_l*normal(:);
Fnl(5) = (U_l(5) + p_l)*vn_l;
Sl = sqrt(gamma*p_l/U_l(1)) + abs(vn_l);

p_r = (gamma-1.0)*(U_r(5) - 0.5*sum(U_r(2:4).^2)/U_r(1));
Fnr(1) = dot(normal(:), U_r(2:4));
vn_r = Fnr(1) / U_r(1);
Fnr(2:4) = vn_r*U_r(2:4) + p_r*normal(:);
Fnr(5) = (U_r(5) + p_r)*vn_r;
Sr = sqrt(gamma*p_r/U_r(1)) + abs(vn_r);

S = max(Sl, Sr);
Fn = 0.5*(Fnl + Fnr - S*(U_r - U_l));

end