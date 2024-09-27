function [x, w, M, Minv, D, Tl, Tr, Pmn, Pnm] = setup_LGLnodal(p)
% [x, w, M, Minv, D, Tl, Tr] = setup_LGLnodal(p)
%   compute the matrices for p-th order LGL nodal space on the 1D interval [-1,1]
%   using Lobatto nodal basis and Legendre modal basis
% 
% input: 
%   p:      polynommial order
% 
% output: 
%   x:      quadrature nodes, size = [p+1, 1]
%   w:      quadrature weights, size = [p+1, 1]
%   M:      mass matrix with respect to nodal DOFs, size = [p+1, p+1]
%   Minv:   inverse of the mass matrix, size = [p+1, p+1]
%   D:      differentiation matrix, size = [p+1, p+1]
%   Tl:     left trace matrix, size = [1, p+1]
%   Tr:     right trace matrix, size = [1, p+1]
%   Pmn:    maps modal coefs to nodal values, size = [p+1, p+1]
%   Pnm:    maps nodal values to modal coefs, size = [p+1, p+1]
% 
% Remark: 
%   SBP property: (D'*diag(w) + diag(w)*D) - (Tr'*Tr - Tl'*Tl) = 0 for 1D and 2D 

% checked. 

[x, w] = gausslobatto(p+1);
Tl = sparse([1],[1],[1.0] ,1,p+1, 1);
Tr = sparse([1],[p+1],[1.0], 1,p+1, 1);

switch (p)
    case(0)
        M = [2.0];
        Minv = [0.5];
        D = [0.0];
        Pmn = [1.0];
        Pnm = [1.0];
        return;
    case (1)
        M = [
            2.0/3, 1.0/3;
            1.0/3, 2.0/3
        ];
        Minv = [
            2, -1;
            -1, 2
        ];
        D = [
            -0.5, 0.5;
            -0.5, 0.5
        ];
        Pmn = [
            1.0, -1.0;
            1.0, 1.0
        ];
        Pnm = [
            0.5, 0.5;
            -0.5, 0.5
        ];
        return;
    case (2)
        M = [
            4.0/15, 2.0/15, -1.0/15;
            2.0/15, 16.0/15, 2.0/15;
            -1.0/15, 2.0/15, 4.0/15
        ];
        Minv = [
            4.5, -0.75, 1.5;
            -0.75, 1.125, -0.75;
            1.5, -0.75, 4.5
        ];
        D = [
            -1.5, 2, -0.5;
            -0.5, 0, 0.5;
            0.5, -2, 1.5
        ];
        Pmn = [
            1.0, -1.0, 1.0;
            1.0, 0.0, -0.5;
            1.0, 1.0, 1.0
        ];
        Pnm = [
            1.0/6, 2.0/3, 1.0/6;
            -0.5, 0.0, 0.5;
            1.0/3, -2.0/3, 1.0/3
        ];
        return;
    case (3)
        M = [
            1.0/7, sqrt(5.0)/42, -sqrt(5.0)/42, 1.0/42;
            sqrt(5.0)/42, 5.0/7, 5.0/42, -sqrt(5.0)/42;
            -sqrt(5.0)/42, 5.0/42, 5.0/7, sqrt(5.0)/42;
            1.0/42, -sqrt(5.0)/42, sqrt(5.0)/42, 1.0/7
        ];
        Minv = [
            8.0, -2.0/sqrt(5.0), 2.0/sqrt(5.0), -2.0;
            -2.0/sqrt(5.0), 1.6, -0.4, 2.0/sqrt(5.0);
            2.0/sqrt(5.0), -0.4, 1.6, -2.0/sqrt(5.0);
            -2.0, 2.0/sqrt(5.0), -2.0/sqrt(5.0), 8.0
        ];
        D = [
            -3.0, 1.25*(1.0 + sqrt(5.0)), -1.25*(-1.0 + sqrt(5.0)), 0.5;
            0.25*(-1.0 - sqrt(5.0)), 0.0, 0.5*sqrt(5.0), 0.25*(1.0 - sqrt(5.0));
            0.25*(-1.0 + sqrt(5.0)), -0.5*sqrt(5.0), 0.0, 0.25*(1.0 + sqrt(5.0));
            -0.5, 1.25*(-1.0 + sqrt(5.0)), -1.25*(1.0 + sqrt(5.0)), 3.0
        ];
        Pmn = [
            1.0, -1.0, 1.0, -1.0;
            1.0, -1.0/sqrt(5.0), -0.2, 1.0/sqrt(5.0);
            1.0, 1.0/sqrt(5.0), -0.2, -1.0/sqrt(5.0);
            1.0, 1.0, 1.0, 1.0
        ];
        Pnm = [
            1.0/12, 5.0/12, 5.0/12, 1.0/12;
            -0.25, -0.25*sqrt(5.0), 0.25*sqrt(5.0), 0.25;
            5.0/12, -5.0/12, -5.0/12, 5.0/12;
            -0.25, 0.25*sqrt(5.0), -0.25*sqrt(5.0), 0.25
        ];
        return;
    case (4)
        M = [
            4.0/45, 7.0/270, -4.0/135, 7.0/270, -1.0/90;
            7.0/270, 196.0/405, 28.0/405, -49.0/810, 7.0/270;
            -4.0/135, 28.0/405, 256.0/405, 28.0/405, -4.0/135;
            7.0/270, -49.0/810, 28.0/405, 196.0/405, 7.0/270;
            -1.0/90, 7.0/270, -4.0/135, 7.0/270, 4.0/45
        ];
        Minv = [
            12.5, -15.0/14, 0.9375, -15.0/14, 2.5;
            -15.0/14, 225.0/98, -45.0/112, 45.0/98, -15.0/14;
            0.9375, -45.0/112, 225.0/128, -45.0/112, 0.9375;
            -15.0/14, 45.0/98, -45.0/112, 225.0/98, -15.0/14;
            2.5, -15.0/14, 0.9375, -15.0/14, 12.5
        ];
        D = [
            -5.0, (7.0/12)*(7.0 + sqrt(21.0)), -8.0/3, -(7.0/12)*(-7.0 + sqrt(21.0)), -0.5;
            -(3.0/28)*(7.0 + sqrt(21.0)), 0.0, 8.0/sqrt(21.0), -0.5*sqrt(7.0/3), -(3.0/28)*(-7.0 + sqrt(21.0));
            0.375, -0.875*sqrt(7.0/3), 0.0, 0.875*sqrt(7.0/3), -0.375;
            (3.0/28)*(-7.0 + sqrt(21.0)), 0.5*sqrt(7.0/3), -8.0/sqrt(21.0), 0.0, (3.0/28)*(7.0 + sqrt(21.0));
            0.5, (7.0/12)*(-7.0 + sqrt(21.0)), 8.0/3, -(7.0/12)*(7.0 + sqrt(21.0)), 5.0
        ];
        Pmn = [
            1.0, -1.0, 1.0, -1.0, 1.0;
            1.0, -sqrt(3.0/7), 1.0/7, (3.0/7)*sqrt(3.0/7), -3.0/7;
            1.0, 0.0, -0.5, 0.0, 0.375;
            1.0, sqrt(3.0/7), 1.0/7, -(3.0/7)*sqrt(3.0/7), -3.0/7;
            1.0, 1.0, 1.0, 1.0, 1.0
        ];
        Pnm = [
            0.05, 49.0/180, 16.0/45, 49.0/180, 0.05;
            -0.15, -0.35*sqrt(7.0/3), 0.0, 0.35*sqrt(7.0/3), 0.15;
            0.25, 7.0/36, -8.0/9, 7.0/36, 0.25;
            -0.35, 0.35*sqrt(7.0/3), 0.0, -0.35*sqrt(7.0/3), 0.35;
            0.2, -7.0/15, 8.0/15, -7.0/15, 0.2
        ];
        return;
    otherwise % p>=5
        % Pmn maps modal coefs to nodal values (exact)
        Pmn = nan(p+1, p+1);
        for i = 1: p+1
            Pmn(:, i) = legendrep(x, i-1);
        end

        % Pnm maps nodal values to modal coefs, size = [p+1, p+1]
        Pnm = inv(Pmn);

        % compute modal matrices (exact)
        [M_m, Minv_m, D_m] = Legendre_basis_mat(p);

        % compute nodal matrices
        D = Pmn * D_m / Pmn; % D = Pmn * D_m * Pnm;

        Minv = Pmn * Minv_m * Pmn';
        Minv = 0.5*(Minv + Minv'); % symmetrize

        M = inv(Minv);  % here, we use inv() since Pmn is exact while Pnm is not
        % faster implementation (not too much difference regarding error)
        % R = chol(Minv); % Minv = R'*R
        % M = R \ inv(R');
        
        % M = Pnm' * M_m * Pnm; % Tests show big error if this is used. 
        M = 0.5*(M + M'); % symmetrize
end

return;
end