function [U_roe, c_roe] = Euler_Roe_average(U_l, U_r, gamma, dim)
% [U_roe, c_roe] = Euler_Roe_average(U_l, U_r, gamma, dim)
%   compute the Roe's average of primitive variables and the sonic speed
%   only support the gamma law EOS
%   (support any problem dimension)
%   (support vectorized input)
% 
% input:
%   U_l, U_r: primitive variables, size = [dim+2, ...]
%   gamma: the specific heat ratio in the gamma law EOS, scalar
%   dim: problem dimension, integer scalar
% 
% output:
%   U_roe:  the Roe average in primitive variables, size = [dim+2, ...] (same as U_l)
%   c_roe:  sonic speed, size = [1, ...]

input_shape = size(U_l); % row vector, input_shape(1) should== dim+2
U_l = reshape(U_l, dim+2,[]);
U_r = reshape(U_r, dim+2,[]);

rho_l = U_l(1,:);
v_l = U_l(2:dim+1,:);
p_l = U_l(dim+2,:);
c_l_sq = gamma*p_l./rho_l; % square of sonic speed

rho_r = U_r(1,:);
v_r = U_r(2:dim+1,:);
p_r = U_r(dim+2);
c_r_sq = gamma*p_r./rho_r; % square of sonic speed

R_roe = sqrt(rho_r./rho_l); % row vector
deno = 1.0 + R_roe;
rho_roe = sqrt(rho_l.*rho_r);
v_roe = (v_l + R_roe.*v_r) ./ deno;

% verified by Mathematica
c_roe_sq = (0.5*(gamma-1.0)./(R_roe + 2.0 + 1.0./R_roe)) .* sum((v_l - v_r).^2, 1) + (c_l_sq + R_roe .* c_r_sq) ./ deno;
c_roe = sqrt(c_roe_sq);
p_roe = (rho_roe .* c_roe_sq) / gamma;

U_roe = [rho_roe; v_roe; p_roe];

U_roe = reshape(U_roe, input_shape);
c_roe = reshape(c_roe, [1,input_shape(2:end)]);

end