function [Lmat, Rmat] = Euler_eigenmat_con(U_p, gamma, dim)
% [Lmat, Rmat] = Euler_eigenmat_con(U_p, gamma)
% compute the left and right eigen matrices in x-direction for compressible Euler equations
% input: primitive variables
% output: matrix for the conservative flux
d = numel(U_p) - 2;
U_p = U_p(:);
rho = U_p(1);
v = U_p(2:(d+1));
vsq = sum(v(:).^2);
p = U_p(end);
T = p/rho;
c = sqrt(gamma*T);
H = (gamma/(gamma-1.0)) * T + 0.5 * vsq;

Lmat = zeros(d+2,d+2);
Rmat = zeros(d+2,d+2);

Lmat(1,1:2) = [-v(1), 1.0] / rho;
Lmat(3,1:2) = Lmat(1,1:2);
coef = (gamma-1.0)/(rho*c);
vec = [0.5*vsq, -v(:)', 1.0];
Lmat(1,:) = Lmat(1,:) - coef*vec;
Lmat(3,:) = Lmat(3,:) + coef*vec;
Lmat(2,:) = (gamma/c)*coef*vec;
Lmat(2,1) = Lmat(2,1) - gamma/rho;

Rmat(:,2) = -(rho/gamma) * [1; v(:); 0.5*vsq];
Rmat([2,end],1) = (0.5*rho) * [1; v(1)];
Rmat([2,end],3) = Rmat([2,end],1);
coef = 0.5*rho/c;
vec = [1.0; v(:); H];
Rmat(:,1) = Rmat(:,1) - coef * vec;
Rmat(:,3) = Rmat(:,3) + coef * vec;

if d > 1
    Lmat(4,[1,3]) = [-v(2), 1] / rho;
    Rmat([3,end],4) = rho * [1; v(2)];
    if d > 2
        Lmat(5,[1,4]) = [-v(3), 1] / rho;
        Rmat([4,end],5) = rho * [1; v(3)];
    end
end

end