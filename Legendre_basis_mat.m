function [M, Minv, D] = Legendre_basis_mat(n)
% [M, Minv, D] = Legendre_basis_mat(n)
%   compute the matrices for Legendre basis on [-1,1]
% 
% input: 
%   n:      max degree, n>=0
% 
% output: 
%   M:      diagonal mass matrix, size = [n+1, n+1]
%   Minv:   inverse of the diagonal mass matrix, size = [n+1, n+1]
%   D:      differentiation matrix, size = [n+1, n+1]

% checked. 

M = sparse(1:n+1, 1:n+1, (1.0)./((1:n+1) - 0.5), n+1, n+1, n+1);
Minv = sparse(1:n+1, 1:n+1, (1:n+1) - 0.5, n+1, n+1, n+1);

% p_{2n+1}'(x) = sum_{k=0}^{n-1} (4k+3)*p_{2k+2}(x)
% p_{2n+2}'(x) = sum_{k=0}^{n} (4k+1)*p_{2k+1}(x)
% (index starts from 1)

% asymptotic fill rate: 25%, small-size
D = zeros(n+1, n+1);
for j = 1: n+1
    for i = 1: j-1
        if (mod(i+j, 2) ~= 0)
            D(i,j) = 2*i - 1;
        end
    end
end
D = sparse(D);

end