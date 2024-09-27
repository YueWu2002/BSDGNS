function [p, c] = Euler_EOS(rho, mx, my, mz, E, gamma)
% [p, c] = Euler_EOS(rho, mx, my, mz, E, gamma)
%   compute the pressure and sonic speed for ideal gas

p = (gamma - 1.0) * (E(:) - 0.5*(mx(:).^2 + my(:).^2 + mz(:).^2)./rho(:));
c = sqrt(gamma * p(:)./rho(:));

end