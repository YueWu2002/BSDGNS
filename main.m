[qx, qw, M, Minv, D, Tl, Tr, Pmn, Pnm] = setup_LGLnodal(p);

init_1Dmesh

gamma = 1.4;

u = nan(p+1, 5, mesh.Nx);

%{
% smooth transport problem
for k = 1: mesh.Nx
    for i = 1: p+1
        u(i,1,k) = 1.0 + 0.2*sin((2*pi) * mesh.xx(i,k));
        u(i,2,k) = 1.0 * u(i,1,k); % velocity_x = 1.0
        u(i,3,k) = 0.0;
        u(i,4,k) = 0.0;
        u(i,5,k) = 1.0/(gamma-1.0) + 0.5*(u(i,2,k)^2 + u(i,3,k)^2 + u(i,4,k)^2) / u(i,1,k); % pressure = 1.0
    end
end
%}

% shock tube problem, check point T=1.3
for k = 1: mesh.Nx
    for i = 1: p+1
        if mesh.xx(i,k) < 0.0
            u(i,1,k) = 0.445;
            u(i,2,k) = 0.698 * u(i,1,k); % velocity_x = 0.698
            u(i,3,k) = 0.0;
            u(i,4,k) = 0.0;
            u(i,5,k) = 3.528/(gamma-1.0) + 0.5*(u(i,2,k)^2 + u(i,3,k)^2 + u(i,4,k)^2) / u(i,1,k); % pressure = 3.528
        else
            u(i,1,k) = 0.5;
            u(i,2,k) = 0.0 * u(i,1,k); % velocity_x = 0.0
            u(i,3,k) = 0.0;
            u(i,4,k) = 0.0;
            u(i,5,k) = 0.571/(gamma-1.0) + 0.5*(u(i,2,k)^2 + u(i,3,k)^2 + u(i,4,k)^2) / u(i,1,k); % pressure = 0.571
        end
    end
end


% LF和HLLC在简单算例收敛阶测试成功

T = 1.3; % final time
need_video = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if need_video
    plot_fps = 100;
    frame_dt = 1.0 / plot_fps;
    ready2plot = true;
    next_frametime = 0.0;

    my_video = VideoWriter('rho.avi');
    my_video.FrameRate = 30; % default
    open(my_video);

    % plot initial condition
    if ready2plot
        frameplot_hand = figure("Visible","off");
        plot(mesh.xx(:,:), squeeze(u(:,1,:)), 'x-');
        grid on;
        axis([mesh.a, mesh.b, 0.0, 1.8]);
        title(['t=',num2str(t)]);
        frame = getframe(frameplot_hand);
        close(frameplot_hand);
        writeVideo(my_video, frame);
        next_frametime = min(next_frametime + frame_dt, T);
        ready2plot = false;
    end
end

tic;
flag = true;
t = 0.0;
while flag
    [ut, max_dt] = DG_1D_semidiscrete(mesh, u, p, qw, Minv, D, Tl, Tr, gamma);
    TD_CFL = 1.0; % CFL from the time discretization
    dt = TD_CFL * max_dt / (2*p+1); % a robust choice
    if need_video
        if t + dt >= next_frametime
            dt = next_frametime - t;
            ready2plot = true;
        end
    end
    if t + dt >= T
        % ready to stop
        dt = T - t;
        flag = false;
    end

    % SSP-RK(3,3) (low-storage implementation)
    u1 = u + dt*ut;
    [ut, ~] = DG_1D_semidiscrete(mesh, u1, p, qw, Minv, D, Tl, Tr, gamma);
    u1 = 0.75*u + 0.25*(u1 + dt*ut);
    [ut, ~] = DG_1D_semidiscrete(mesh, u1, p, qw, Minv, D, Tl, Tr, gamma);
    u = (1.0/3)*u + (2.0/3)*(u1 + dt*ut);
    t = t + dt;

    if need_video
        if ready2plot
            frameplot_hand = figure("Visible", "off");
            plot(mesh.xx(:,:), squeeze(u(:,1,:)), 'x-');
            grid on;
            axis([mesh.a, mesh.b, 0.0, 1.8]);
            title(['t=', num2str(t)]);
            frame = getframe(frameplot_hand);
            close(frameplot_hand);
            writeVideo(my_video, frame);
            next_frametime = min(next_frametime + frame_dt, T);
            ready2plot = false;
        end
    end

    disp(['t=', num2str(t)]);
end
toc % output the time consumption
if need_video
    close(my_video);
end


%{
% compute error for the smooth transport problem
ue = nan(p+1, 5, mesh.Nx);
for k = 1: mesh.Nx
    for i = 1: p+1
        ue(i,1,k) = 1.0 + 0.2*sin((2*pi) * (mesh.xx(i,k) - T));
        ue(i,2,k) = 1.0 * ue(i,1,k);
        ue(i,3,k) = 0.0;
        ue(i,4,k) = 0.0;
        ue(i,5,k) = 1.0;
    end
end

err = sqrt(sum( squeeze(abs(u(:,1,:) - ue(:,1,:)).^2) .* ((0.5*(mesh.b-mesh.a)/mesh.Nx) * qw(:)), 'all'))
%}

figure(1);
plot(mesh.xx(:,:), squeeze(u(:,1,:)), 'x-');
grid on;
title('density');

for i = 1:3
    figure(i+1);
    plot(mesh.xx(:,:), squeeze(u(:,i+1,:)./u(:,1,:)), 'x-');
    grid on;
    title(['v_',num2str(i)]);
end

figure(5);
plot(mesh.xx(:,:), squeeze( (gamma-1.0)*(u(:,5,:) - 0.5*(u(:,2,:).^2 + u(:,3,:).^2 + u(:,4,:).^2)./u(:,1,:)) ), 'x-');
grid on;
title('pressure');