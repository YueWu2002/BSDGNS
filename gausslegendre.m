function [x, w] = gausslegendre(n)
% [x, w] = gausslegendre(n)
%   compute the n-point Gauss-Legendre quadrature points and weights on [-1,1]
%   algebraic exactness: 2n-1, convergence rate when using composite integration: 2n
% 
% input: 
%   n:  number of points, n>=1
% 
% output: (exact for 1<=n<=5)
%   x:  quadrature nodes, size = [n, 1]
%   w:  quadrature weights, size = [n, 1]
% 
% references: 
% [1] Shen, J., Tang, T., Wang, LL. (2011). Orthogonal Polynomials and Related Approximation Results. In: Spectral Methods. Springer Series in Computational Mathematics, vol 41. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-71041-7_3.
% [2] https://mathworld.wolfram.com/Legendre-GaussQuadrature.html

% checked. 

switch (n)
    case (1)
        x = [0.0];
        w = [2.0];
        return;
    case (2)
        x = [
            -sqrt(1.0/3);
            sqrt(1.0/3)
        ];
        w = [
            1.0;
            1.0
        ];
        return;
    case (3)
        x = [
            -sqrt(0.6);
            0.0;
            sqrt(0.6)
        ];
        w = [
            5.0/9;
            8.0/9;
            5.0/9
        ];
        return;
    case (4)
        x = [
            -sqrt((3.0 + 2*sqrt(1.2))/7);
            -sqrt((3.0 - 2*sqrt(1.2))/7);
            sqrt((3.0 - 2*sqrt(1.2))/7);
            sqrt((3.0 + 2*sqrt(1.2))/7)
        ];
        w = [
            0.5 - sqrt(30.0)/36;
            0.5 + sqrt(30.0)/36;
            0.5 + sqrt(30.0)/36;
            0.5 - sqrt(30.0)/36
        ];
        return;
    case (5)
        x = [
            -sqrt((35.0 + 2*sqrt(70.0))/63);
            -sqrt((35.0 - 2*sqrt(70.0))/63);
            0.0;
            sqrt((35.0 - 2*sqrt(70.0))/63);
            sqrt((35.0 + 2*sqrt(70.0))/63)
        ];
        w = [
            (322.0 - 13*sqrt(70.0))/900;
            (322.0 + 13*sqrt(70.0))/900;
            128.0/225;
            (322.0 + 13*sqrt(70.0))/900;
            (322.0 - 13*sqrt(70.0))/900
        ];
        return;
    otherwise
        % solve the nodes using the eigenvalue method (only for small-size problems)
        beta = zeros(1, n-1);
        for j = 1: n-1
            beta(j) = 1.0 / sqrt(4.0 - 1.0/j^2);
        end
        x = eig(sparse([(1:n-1), (2:n)], [(2:n), (1:n-1)], [beta, beta], n, n));
        x = 0.5*(x(1:end) - x(end:-1:1)); % anti-symmetrize
        
        % directly solve the weights using Taylor basis (deprecated)
        % bad conditioning and big error for big n: 
        %   cond(vander(x))=1.1243e+07 when n=20
        %   cond(vander(x))=4.9550e+14 when n=40
        %{
        b = zeros(n, 1);
        b(1:2:n) = (2.0)./(1:2:n);
        w = fliplr(vander(x))'\b;
        %}
        % remark: since the quadrature points have error, the integration accuracy will be the best if we use the above method to compute the weights, given n is not too big. 

        % directly solve the weights using Legendre basis (recommended)
        % very good conditioning: 
        %   cond(mat)=7.9122 when n=20
        %   cond(mat)=11.6409 when n=40
        %   cond(mat)=16.9256 when n=80
        b = sparse([1],[1],[2.0], n,1, 1);
        mat = nan(n,n);
        mat(1,:) = 1.0;
        if n >= 2
            mat(2,:) = x(:);
            for k = 3:n
                alpha = 1.0/(k-1);
                for j = 1: n
                    mat(k,j) = ((2.0-alpha)*x(j)).*mat(k-1,j) - (1.0-alpha)*mat(k-2,j);
                end
            end
        end
        w = mat\b;
        % remark: since the quadrature points have error, the integration accuracy will be the best if we use the above method to compute the weights, given n is not too big. 
end

end