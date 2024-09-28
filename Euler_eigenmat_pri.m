function [Lmat, Rmat] = Euler_eigenmat_pri(U_p, gamma)
% [Lmat, Rmat] = Euler_eigenmat_pri(U_p, gamma)
% compute the left and right eigen matrices in x-direction for compressible Euler equations
% input: primitive variables
% output: matrix for the primitive flux
d = numel(U_p) - 2;
U_p = U_p(:);
rho = U_p(1);
v = U_p(2:(d+1));
p = U_p(end);
T = p/rho;
c = sqrt(gamma*T);

Lmat = zeros(d+2,d+2);
Rmat = zeros(d+2,d+2);

temp = 1.0/(rho*c);
Lmat([1,3],[2,end]) = [
    1.0, -temp;
    1.0, temp
];
Lmat(2,[1,end]) = [-gamma/rho, 1.0/p];

temp1 = 0.5*rho/c;
temp2 = 0.5*rho*c;
Rmat([1,2,end],1) = [-temp1; 0.5; -temp2];
Rmat(1,2) = -rho/gamma;
Rmat([1,2,end],3) = [temp1; 0.5; temp2];

if d > 3
    Lmat(4,3) = 1.0;
    Rmat(3,4) = 1.0;
    if d > 4
        Lmat(5,4) = 1.0;
        Rmat(4,5) = 1.0;
    end
end

end