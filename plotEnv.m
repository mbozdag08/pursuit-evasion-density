function [] = plotEnv(idx, Init, Rx1, Rx2, Xe, Xp, xrArc, f_fun, ge, u_fun, rho_fun, quiver_density, dirname, vis, save)
% Plot the static environtment
%   idx:                The index of the timepoint in the solution array
%   Init:               Initial parameters
%   Rx1 & Rx2:          2D-grid values (x and y axis)
%   Xe and Xp:          Evader and pursuer trajectories
%   xrArc:              The reach (target) arc
%   f_fun, g, u_fun:    The system x' := f(x) + g(x)u(x) (g is constant)
%   rho_fun:            Density function rho
%   quiver_density:     The density of arrows for the phase portrait
%   dirname:            Directory for the save
%   vis:                Visibility flag for the figures (on/off)
%   save:               Save flag for the figure
    
    fprintf("\nPlotting the environment...\n");
    tic

    % Get initial conditions
    xr0 = Init.xr0;
    Rie = Init.Rie;
    Rip = Init.Rip;
    Rp = Init.Rp;
    Rr = Init.Rr;
    R = Init.R;

    % Extract the intersection points between Xt and Xe
    ptA = double([xrArc.x1(1); xrArc.x2(1)]);
    ptB = double([xrArc.x1(2); xrArc.x2(2)]);
    % Sort by angle to keep consistent order
    angles = [atan2(ptA(2), ptA(1)), atan2(ptB(2), ptB(1))];
    [~, sort_idx] = sort(mod(angles, 2*pi));  % modulo for wraparound
    pts = [ptA, ptB];
    xt1 = pts(:, sort_idx(1));
    xt2 = pts(:, sort_idx(2));
    % Calculate the angles for the big circle (centered at 0,0)
    thetaR1 = atan2(xt1(2), xt1(1));
    thetaR2 = atan2(xt2(2), xt2(1));
    % Calculate the angles for the small circle (centered at xt0)
    thetaRt1 = atan2(xt1(2) - xr0(2), xt1(1) - xr0(1));
    thetaRt2 = atan2(xt2(2) - xr0(2), xt2(1) - xr0(1));

    fig1 = figure('visible', vis);
    hold on;
    axis("equal");
    xticks(-R-1:1:R+1);
    yticks(-R-1:1:R+1);
    xlim([-R-1, R+1]);
    ylim([-R-1, R+1]);

    ax = gca;  % Get current axes
    ax.FontSize = 13;  % Set the axis tick label font size
    ax.LabelFontSizeMultiplier = 1;  % Scale label font size
    ax.TitleFontSizeMultiplier = 1;  % Scale title font size (if you add one)
    
    % Plot the initial and final positions
    plot(Xe(1,1),Xe(1,2), 'b.');                            % Initial position of evader
    plot(Xe(idx,1),Xe(idx,2), 'b.', 'MarkerSize', 15);      % Current position of evader
    plot(Xp(1,1),Xp(1,2), 'r.');                            % Initial position of pursuer
    plot(Xp(idx,1),Xp(idx,2), 'r.', 'MarkerSize', 15);      % Current position of pursuer
    
    % Plot the trajectories
    plot(Xe(1:idx,1), Xe(1:idx,2), 'b-', 'LineWidth', 1.6);  % Evader trajectory
    plot(Xp(1:idx,1), Xp(1:idx,2), 'r-', 'LineWidth', 1.6);  % Pursuer trajectory
    
    % Initial regions for evader and pursuer
    rectangle('Position', [Xp(1,1)-Rip, Xp(1,2)-Rip, 2*Rip, 2*Rip], 'Curvature', [1, 1], 'EdgeColor', 'r', 'LineStyle', '-.', 'LineWidth', 1);
    rectangle('Position', [Xe(1,1)-Rie, Xe(1,2)-Rie, 2*Rie, 2*Rie], 'Curvature', [1, 1], 'EdgeColor', 'b', 'LineStyle', '-.', 'LineWidth', 1);
    % Plot a red circle around the pursuer
    rectangle('Position', [Xp(idx,1)-Rp, Xp(idx,2)-Rp, 2*Rp, 2*Rp], 'Curvature', [1, 1], 'EdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 1.2);
    % Plot the target arc
    plotArc(0, 0, R, thetaR1, thetaR2, 'g-', 1.8);
    plotArc(xr0(1), xr0(2), Rr, thetaRt1, thetaRt2, 'g-', 1.8);
    % Plot the bounded region
    plotArc(0, 0, R, thetaR2, thetaR1 + 2*pi, 'r-', 1.8);

    % Plot rho
    % Grid X for Xp(idx)
    Xgrid = [Rx1(:), Rx2(:), Xp(idx,1)*ones(numel(Rx1), 1), Xp(idx,2)*ones(numel(Rx2), 1)];
    % Calculate rho to plot the rho=0
    rho_arr = rho_fun(Xgrid');
    % Reshape rho_eval to fit the grid dimensions
    rho_grid = reshape(rho_arr, size(Rx1));
    % Add contour lines for rho=0 with a dashed line
    hold on;
    contour(Rx1, Rx2, rho_grid, [0,0], 'k--', 'LineWidth', 1.3);
    
    % Plot the phase portrait for xp(t_curr)
    x_dot = @(t,x) f_fun([x(1);x(2);Xp(idx,1);Xp(idx,2)]) + [ge * u_fun([x(1);x(2);Xp(idx,1);Xp(idx,2)]); 0; 0];
    plotpp(x_dot, 'quiverDensity', quiver_density);
    
    if save == 1
        set(fig1, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
        savefig(fig1,sprintf("%s/PE_idx_%d_qDen_%d.fig",dirname, idx, quiver_density));
        print(fig1, sprintf("%s/PE_idx_%d_qDen_%d.png", dirname, idx, quiver_density), '-dpng', '-r300'); % Save with 300 DPI resolution
    end
    hold off;
   
    toc

end