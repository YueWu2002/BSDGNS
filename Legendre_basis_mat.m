function [M, Minv, D, Int] = Legendre_basis_mat(n)
% [M, Minv, D, Int] = Legendre_basis_mat(n)
%   compute the matrices for Legendre basis on [-1,1]
% 
% input: 
%   n:      max degree, n>=0
% 
% output: 
%   M:      diagonal mass matrix, size = [n+1, n+1]
%   Minv:   inverse of the diagonal mass matrix (elements are integer times halfs), size = [n+1, n+1]
%   D:      differentiation matrix (elements are integers), size = [n+1, n+1]
%   Int:    indefinite integration matrix, size = [n+2, n+1]

% Legendre polynomials: 
% (n+1)*p_{n+1}(x) = (2n+1)*p_{n}(x) - n*p_{n-1}(x)
% p_{2n+1}'(x) = sum_{k=0}^{n-1} (4k+3)*p_{2k+2}(x)
% p_{2n+2}'(x) = sum_{k=0}^{n} (4k+1)*p_{2k+1}(x)
% (2n+1)*p_{n}(x) = p_{n+1}'(x) - p_{n-1}'(x)
% 
% checked. 

M = sparse(1:n+1, 1:n+1, (1.0)./((1:n+1) - 0.5), n+1, n+1, n+1);
Minv = sparse(1:n+1, 1:n+1, (1:n+1) - 0.5, n+1, n+1, n+1);

% asymptotic fill rate: 25%
D = zeros(n+1, n+1);
for j = 1: n+1
    for i = 1: j-1
        if (mod(i+j, 2) ~= 0)
            D(i,j) = 2*i - 1;
        end
    end
end
D = sparse(D);


% indefinite integral matrix, the j-th column are the expansion coefficients of indefinite integrals of p_{j-1}
% asymptotic fill rate: 0%
% should have: D(1:n+1,1:n+1)*Int(1:n+1,1:n) == speye(n), (Int*D)(2:end,2:end) == speye(n)
Int = sparse([2:n+2, 1:n], [1:n+1, 2:n+1], [(1.0)./(2*(1:n+1)-1), (-1.0)./(2*(2:n+1)-1)], n+2, n+1, 2*n+1);

end