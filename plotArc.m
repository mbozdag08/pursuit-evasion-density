function plotArc(xc, yc, r, theta_start, theta_end, color, linewidth)
    % Ensure angles are between 0 and 2*pi
    theta_start = mod(theta_start, 2*pi);
    theta_end = mod(theta_end, 2*pi);
    
    % If the end angle is less than the start angle, wrap around the circle
    if theta_end < theta_start
        theta_end = theta_end + 2*pi;
    end
    
    % Define the arc range with 100 points for smoothness
    theta = linspace(theta_start, theta_end, 100);
    
    % Compute x and y coordinates of the arc
    x_arc = xc + r * cos(theta);
    y_arc = yc + r * sin(theta);
    
    % Plot the arc
    plot(x_arc, y_arc, color, 'LineWidth', linewidth);
end
