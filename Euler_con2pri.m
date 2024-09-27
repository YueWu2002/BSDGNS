function [U_p] = Euler_con2pri(U_c, gamma)
% from conservative to primitive
U_p = [U_c(1); U_c(2)/U_c(1); U_c(3)/U_c(1); U_c(4)/U_c(1); (gamma-1.0)*(U_c(5) - 0.5*(U_c(2)^2 + U_c(3)^2 + U_c(4)^2)/U_c(1))];
end