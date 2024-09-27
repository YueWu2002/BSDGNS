function [x, w] = gausslobatto(n)
% [x, w] = gausslobatto(n)
%   compute the n-point Gauss-Lobatto quadrature points and weights on [-1,1]
%   algebraic exactness: 2n-3, convergence rate when using composite integration: 2n-2
% 
% input: 
%   n:  number of points, n>=2
% 
% output: (exact for 1<=n<=7)
%   x:  quadrature nodes, size = [n, 1]
%   w:  quadrature weights, size = [n, 1]
% 
% references: 
% [1] Shen, J., Tang, T., Wang, LL. (2011). Orthogonal Polynomials and Related Approximation Results. In: Spectral Methods. Springer Series in Computational Mathematics, vol 41. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-71041-7_3.
% [2] https://mathworld.wolfram.com/LobattoQuadrature.html

% checked.

switch (n)
    case(1)
        disp('gausslobatto: input n=1, use gauss-legendre node instead.');
        x = [0.0];
        w = [2.0];
        return;
    case (2)
        x = [
            -1.0;
            1.0
        ];
        w = [
            1.0;
            1.0
        ];
        return;
    case (3)
        x = [
            -1.0;
            0.0;
            1.0
        ];
        w = [
            1.0/3;
            4.0/3;
            1.0/3
        ];
        return;
    case (4)
        x = [
            -1.0;
            -sqrt(0.2);
            sqrt(0.2);
            1.0
        ];
        w = [
            1.0/6;
            5.0/6;
            5.0/6;
            1.0/6
        ];
        return;
    case (5)
        x = [
            -1.0;
            -sqrt(3.0/7);
            0.0;
            sqrt(3.0/7);
            1.0
        ];
        w = [
            0.1;
            49.0/90;
            32.0/45;
            49.0/90;
            0.1
        ];
        return;
    case (6)
        x = [
            -1.0;
            -sqrt((7.0 + 2*sqrt(7.0))/21);
            -sqrt((7.0 - 2*sqrt(7.0))/21);
            sqrt((7.0 - 2*sqrt(7.0))/21);
            sqrt((7.0 + 2*sqrt(7.0))/21);
            1.0
        ];
        w = [
            1.0/15;
            (14.0 - sqrt(7.0))/30;
            (14.0 + sqrt(7.0))/30;
            (14.0 + sqrt(7.0))/30;
            (14.0 - sqrt(7.0))/30;
            1.0/15
        ];
        return;
    case (7)
        x = [
            -1.0;
            -sqrt((15.0 + 2*sqrt(15.0))/33);
            -sqrt((15.0 - 2*sqrt(15.0))/33);
            0.0;
            sqrt((15.0 - 2*sqrt(15.0))/33);
            sqrt((15.0 + 2*sqrt(15.0))/33);
            1.0
        ];
        w = [
            1.0/21;
            (124.0 - 7*sqrt(15.0))/350;
            (124.0 + 7*sqrt(15.0))/350;
            256.0/525;
            (124.0 + 7*sqrt(15.0))/350;
            (124.0 - 7*sqrt(15.0))/350;
            1.0/21
        ];
        return;
    otherwise
        % solve the nodes using the eigenvalue method (only for small-size problems)
        beta = zeros(1, n-3);   % should require n>=4
        for j = 1: n-3
            beta(j) = sqrt((j*(j+2))/((2*j+1)*(2*j+3)));
        end
        x = eig(sparse([(1:n-3), (2:n-2)], [(2:n-2), (1:n-3)], [beta, beta], n-2, n-2));
        [x, ~] = sort(x);
        x = 0.5*(x(1:end) - x(end:-1:1)); % anti-symmetrize
        x = [-1.0; x(:); 1.0];
        
        % directly solve the weights using Taylor basis (deprecated)
        % bad conditioning and big error for big n: 
        %   cond(vander(x))=8.1200e+06 when n=20
        %   cond(vander(x))=3.5500e+14 when n=40
        %{
        b = zeros(n, 1);
        b(1:2:n) = (2.0)./(1:2:n);
        w = fliplr(vander(x))'\b;
        %}
        % remark: since the quadrature points have error, the integration accuracy will be the best if we use the above method to compute the weights, given n is not too big.

        % directly solve the weights using Legendre basis (recommended)
        % very good conditioning: 
        %   cond(mat)=8.1402 when n=20
        %   cond(mat)=12.1889 when n=40
        %   cond(mat)=17.8260 when n=80
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