function [U_p] = Euler_con2pri(U_c, gamma)
% from conservative to primitive
% can adapt to dimension automatically
d = numel(U_c) - 2;
U_c = U_c(:);
U_p = [U_c(1); U_c(2:(d+1))/U_c(1); (gamma-1.0)*(U_c(end) - 0.5*sum(U_c(2:(d+1)).^2)/U_c(1))];
end