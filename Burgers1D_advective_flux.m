function [F, S] = Burgers1D_advective_flux(u)
F = 0.5*u.^2;
S = abs(u);
end