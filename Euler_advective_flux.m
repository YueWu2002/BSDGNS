function [F] = Euler_advective_flux(U, p, normal, dim)
% [F] = Euler_advective_flux(U, p, normal, dim)
%   compute the advective flux of the Euler equation in conservative form in the normal direction
%   (support any equation of state)
%   (support any problem dimension)
%   (support vectorized input)
% 
% input: 
%   U:          conservative variables, size = [dim+2, ...]
%   p:          pressure,               size = [1, ...]
%   normal:     unit normal vector,     size = [dim, ...]
%   dim:        problem dimension,      integer scalar
% 
% output:
%   F:          flux in the normal direction, size = [dim+2, ...] (same as U)

input_shape = size(U); % row vector, input_shape(1) should== dim+2
U = reshape(U, dim+2,[]);
p = reshape(p, 1,[]); % row vector
normal = reshape(normal, dim,[]);
vn = dot(U(2:dim+1,:)./U(1,:), normal(1:dim,:), 1); % normal velocity (row vector)
F = reshape(vn .* U(:,:) + p .* [zeros(1,size(U,2)); normal(1:dim,:); vn], input_shape);
end