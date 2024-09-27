function [U_c] = Euler_pri2con(U_p, gamma)
% from primitive to conservative
U_c = [U_p(1); U_p(1)*U_p(2); U_p(1)*U_p(3); U_p(1)*U_p(4); U_p(5)/(gamma-1.0) + 0.5*U_p(1)*(U_p(2)^2 + U_p(3)^2 + U_p(4)^2)];
end