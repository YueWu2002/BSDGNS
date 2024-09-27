function [Fn, S] = Burgers1D_Riemann_flux(u_l, u_r)
% [Fn, S] = Burgers1D_Riemann_flux(u_l, u_r)
% 

if u_l >= u_r
    % shock or contact discontinuiy
    u_m = 0.5*(u_l + u_r); % shock speed
    if u_m < 0.0
        % use the right state
        Fn = 0.5*u_r^2;
    else
        % use the left state
        Fn = 0.5*u_l^2;
    end
else
    % rarefaction
    u_m = 0.0;
    if u_r <= 0.0
        % use the right state
        Fn = 0.5*u_r^2;
    elseif u_l >= 0.0
        % use the left state
        Fn = 0.5*u_l^2;
    else
        % in the rarefaction
        Fn = 0.0;
    end
end

S = max(abs([u_l, u_m, u_l]));
end