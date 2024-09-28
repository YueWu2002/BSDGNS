function [U_c] = Euler_pri2con(U_p, gamma)
% from primitive to conservative
% can adapt to dimension automatically
d = numel(U_p) - 2;
U_p = U_p(:);
U_c = [U_p(1); U_p(1)*U_p(2:(d+1)); U_p(end)/(gamma-1.0) + 0.5*U_p(1)*sum(U_p(2:(d+1)).^2)];
end