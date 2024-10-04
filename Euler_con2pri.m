function [U_p, sonic] = Euler_con2pri(U_c, gamma, dim)
% [U_p, sonic] = Euler_con2pri(U_c, gamma, dim)
%   transform conservative variables to primitive variables and compute the speed of sound
%   for compressible Euler equations of dim-dimension with gamma law EOS
%   (support any problem dimension)
%   (support vectorized input)
% 
% input:
%   U_c:    conservative variables, size = [dim+2, ...]
%   gamma:  the specific heat ratio in the gamma law EOS, scalar
%   dim:    problem dimension,      integer scalar
% 
% output:
%   U_p:    primitive variables,    size = [dim+2, ...] (same as U_c)
%   sonic:  speed of sound,         size = [1, ...]

input_shape = size(U_c); % row vector, input_shape(1) should== dim+2
U_c = reshape(U_c, dim+2,[]);
U_p = [
    U_c(1,:);
    U_c(2:dim+1,:)./U_c(1,:);
    (gamma-1.0)*(U_c(dim+2,:) - 0.5*sum(U_c(2:dim+1,:).^2, 1)./U_c(1,:))
];
if nargout > 1
    sonic = sqrt(gamma*U_p(dim+2,:)./U_p(1,:));
    sonic = reshape(sonic, [1,input_shape(2:end)]);
end
U_p = reshape(U_p, input_shape);
end