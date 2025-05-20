function [] = animateEnv(t_sol, Init, Rx1, Rx2, Xe, Xp, xrArc, f_fun, ge, u_fun, rho_fun, quiver_density, frameRate, quality, dirname)
% Animate the dynamical system
%   t_sol:              Time array from the ODE solver
%   Init:               Initial parameters
%   Rx1 & Rx2:          2D-grid values (x and y axis)
%   Xe and Xp:          Evader and pursuer trajectories
%   xrArc:              The safe arc (reach set)
%   f_fun, g, u_fun:    The system x' := f(x) + g(x)u(x) (g is constant)
%   rho_fun:            Density function rho
%   quiver_density:     The density of arrows for the phase portrait
%   frameRate:          Frame rate for the video
%   quality:            Quality of the video
%   dirname:            Directory for the save
        
    fprintf("\nAnimation...\n");
    close all;
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

    % Index the time array
    indexes = 1:length(t_sol); % Or specify a specific range of t_curr values
    
    % Create a video writer object
    v = VideoWriter(sprintf("%s/Animation_frate_%d_quality_%d.mp4", dirname, frameRate, quality), 'MPEG-4');
    v.FrameRate = frameRate; % Set desired frame rate
    v.Quality = quality; % Set video quality (0 to 100)
    open(v);
    
    % Prepare figure and static elements
    fig_frame = figure('Position', [100, 100, 800, 600], 'Visible','off');
    
    % Create main axes for trajectory and other plots
    main_axes = axes('Parent', fig_frame);
    hold(main_axes, 'on');
    axis(main_axes, 'equal');
    axis(main_axes, 1.1*[-R R -R R]);
    grid(main_axes, 'off');
    phase_axes = axes('Parent', fig_frame, 'color', 'none');
    hold(phase_axes, 'on');
    axis(phase_axes, 'equal');
    axis(phase_axes, 1.1*[-R R -R R]);
    grid(phase_axes, 'off');
    
    % To plot dynamic stuff
    set(fig_frame, 'CurrentAxes', main_axes);
    
    % Plot static elements
    plot(xr0(1), xr0(2), 'g.'); % Target zone center
    plot(Xe(1, 1), Xe(1, 2), 'b.'); % Initial position of evader
    plot(Xp(1, 1), Xp(1, 2), 'r.'); % Initial position of pursuer
    % Initial region of evader
    rectangle('Position', [Xp(1, 1) - Rip, Xp(1, 2) - Rip, 2 * Rip, 2 * Rip], 'Curvature', [1, 1], 'EdgeColor', 'r', 'LineStyle', '-.', 'LineWidth', .7); 
    % Initial region of pursuer
    rectangle('Position', [Xe(1, 1) - Rie, Xe(1, 2) - Rie, 2 * Rie, 2 * Rie], 'Curvature', [1, 1], 'EdgeColor', 'b', 'LineStyle', '-.', 'LineWidth', .7); 
    % Target region
    plotArc(0, 0, R, thetaR1, thetaR2, 'g-', 1.5);
    plotArc(xr0(1), xr0(2), Rr, thetaRt1, thetaRt2, 'g-', 1.5);
    % Bounded region
    plotArc(0, 0, R, thetaR2, thetaR1 + 2*pi, 'r-', 1.5);

    
    % Plot and create and object for the dynamic unsafe zone
    unsafe_zone = rectangle('Position', [Xp(1, 1) - Rp, Xp(1, 2) - Rp, 2 * Rp, 2 * Rp], 'Curvature', [1, 1], 'EdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 1); 

    % Plot the initial rho contour
    Xgrid = [Rx1(:), Rx2(:), Xp(1, 1) * ones(numel(Rx1), 1), Xp(1, 2) * ones(numel(Rx2), 1)];
    rho_arr = rho_fun(Xgrid');
    rho_grid = reshape(rho_arr, size(Rx1));
    [~, rho_contour_object] = contourf(Rx1, Rx2, rho_grid, [0, 0], 'k--', 'LineWidth', 1, 'FaceAlpha', 0);
    
    % Plot the initial Phase Plot
    set(fig_frame, 'CurrentAxes', phase_axes);
    plotpp(@(t, x) f_fun([x(1); x(2); Xp(1, 1); Xp(1, 2)]) + ...
               [ge * u_fun([x(1); x(2); Xp(1, 1); Xp(1, 2)]); 0; 0],  ...
               'xlim', 1.1*[-R, R], 'ylim', 1.1*[-R, R], 'quiverDensity', quiver_density);
    set(fig_frame, 'CurrentAxes', main_axes);
    
    fprintf("________________ \n");
    fprintf("frame: 1 \n");
    
    % Save the figure as an image and capture the frame
    img_file = sprintf('frame_1.png');
    print(fig_frame, img_file, '-dpng', '-r300'); % Save at 300 DPI for high quality
    img = imread(img_file);
    
    % Resize the image to even dimensions if necessary
    [rows, cols, ~] = size(img);
    if mod(rows, 2) ~= 0
        rows = rows + 1; % Make the number of rows even
    end
    if mod(cols, 2) ~= 0
        cols = cols + 1; % Make the number of columns even
    end
    img_resized = imresize(img, [rows, cols]);
    
    % Convert the image to a frame and write to the video
    writeVideo(v, im2frame(img_resized));
    
    % Delete the temporary image file
    delete(img_file);
    
    % Loop over each t_curr to create frames
    for k = 2:length(indexes)
    
        fprintf("frame: %d \n", k);
    
        idx = indexes(k);
        t_prev = indexes(k-1);
    
        % Plot the trajectories
        plot(Xe(t_prev:idx, 1), Xe(t_prev:idx, 2), 'b-', 'LineWidth', 1.2); % Evader trajectory
        plot(Xp(t_prev:idx, 1), Xp(t_prev:idx, 2), 'r-', 'LineWidth', 1.2); % Pursuer trajectory
    
        % Update the pursuer's unsafe zone
        set(unsafe_zone, 'Position', [Xp(idx, 1) - Rp, Xp(idx, 2) - Rp, 2 * Rp, 2 * Rp]);
    
        % Update rho Contour
        Xgrid = [Rx1(:), Rx2(:), Xp(idx, 1) * ones(numel(Rx1), 1), Xp(idx, 2) * ones(numel(Rx2), 1)];
        rho_arr = rho_fun(Xgrid');
        rho_grid = reshape(rho_arr, size(Rx1));
        set(rho_contour_object, 'ZData', rho_grid);
    
        % Update the phase plot
        set(fig_frame, 'CurrentAxes', phase_axes);
        cla(phase_axes);
        plotpp(@(t, x) f_fun([x(1); x(2); Xp(idx, 1); Xp(idx, 2)]) + ...
               [ge * u_fun([x(1); x(2); Xp(idx, 1); Xp(idx, 2)]); 0; 0], ...
               'xlim', 1.1*[-R, R], 'ylim', 1.1*[-R, R], 'quiverDensity', quiver_density);
        set(fig_frame, 'CurrentAxes', main_axes);
    
        % Save the figure as an image and capture the frame
        img_file = sprintf('frame_%d.png', k);
        print(fig_frame, img_file, '-dpng', '-r300'); % Save at 300 DPI for high quality
        img = imread(img_file);
    
        % Resize the image to even dimensions if necessary
        [rows, cols, ~] = size(img);
        if mod(rows, 2) ~= 0
            rows = rows + 1; % Make the number of rows even
        end
        if mod(cols, 2) ~= 0
            cols = cols + 1; % Make the number of columns even
        end
        img_resized = imresize(img, [rows, cols]);
    
        % Convert the image to a frame and write to the video
        writeVideo(v, im2frame(img_resized));
    
        % Delete the temporary image file
        delete(img_file);
    
    end
    fprintf("________________ \n");
    
    % Close the video writer
    close(fig_frame);
    close(v);
    toc

end