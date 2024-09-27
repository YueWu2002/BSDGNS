function [Fx, Fy, Fz] = Euler_advective_flux(U, p)
% [Fx, Fy, Fz] = Euler_advective_flux(U, p)
%   compute the advective flux of the Euler equation in conservative form
%   (for general equation of state)
% 
% input: (scalar)
%   U:          conservative variables
%   p:          pressure
% 
% output:
%   Fx, Fy, Fz: flux in x,y,z directions

U = U(:);

Fx = nan(5,1);
vx = U(2)/U(1);
Fx(1) = U(2);
Fx(2) = vx*U(2) + p;
Fx(3) = vx*U(3);
Fx(4) = vx*U(4);
Fx(5) = (U(5) + p)*vx;

if nargout > 1
    Fy = nan(5,1);
    vy = U(3)/U(1);
    Fy(1) = U(3);
    Fy(2) = vy*U(2);
    Fy(3) = vy*U(3) + p;
    Fy(4) = vy*U(4);
    Fy(5) = (U(5) + p)*vy;

    if nargout > 2
        Fz = nan(5,1);
        vz = U(4)/U(1);
        Fz(1) = U(4);
        Fz(2) = vz*U(2);
        Fz(3) = vz*U(3);
        Fz(4) = vz*U(4) + p;
        Fz(5) = (U(5) + p)*vz;
    end
end

end