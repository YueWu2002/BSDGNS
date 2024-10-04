function [U_c] = Euler_pri2con(U_p, gamma, dim)
% [U_c] = Euler_pri2con(U_p, gamma, dim)
%   transform primitive variables to conservative variables
%   for compressible Euler equations of dim-dimension with gamma law EOS
%   (support any problem dimension)
%   (support vectorized input)
% 
% input:
%   U_p:    primitive variables,    size = [dim+2, ...]
%   gamma:  the specific heat ratio in the gamma law EOS, scalar
%   dim:    problem dimension,      integer scalar
% 
% output:
%   U_c:    conservative variables, size = [dim+2, ...] (same as U_p)

input_shape = size(U_p); % row vector, input_shape(1) should== dim+2
U_p = reshape(U_p, dim+2,[]);
U_c = reshape([U_p(1,:); U_p(1,:).*U_p(2:dim+1,:); U_p(dim+2,:)/(gamma-1.0) + 0.5*U_p(1,:).*sum(U_p(2:dim+1,:).^2, 1)], input_shape);
end
