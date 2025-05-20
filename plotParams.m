function [] = plotParams(Init, Xe, Xp, dirname, vis, save)
% Plot some system parameters
%   Init:               Initial parameters
%   Xe and Xp:          Evader and pursuer trajectories
%   dirname:            Directory for the save
%   vis:                Visibility flag for the figures (on/off)
%   save:               Save flag for the figure (1/0)
    
    fprintf("\nPlotting parameters...\n");
    tic

    % Get initial conditions
    xr0 = Init.xr0;
    Rp = Init.Rp;
    Rr = Init.Rr;
    R = Init.R;
    
    fontSize = 20;
    
    % Distance to X
    figX = figure('visible',vis);
    hold on;
    plot(vecnorm(Xe,2,2), 'c', 'LineWidth', 1.6);
    plot(vecnorm(Xp,2,2), 'r', 'LineWidth', 1.6);
    yline(R,'r--', 'LineWidth', 1.6);
    xlim([0 size(Xe,1)]);
    ylim([0 inf]);
    legend("x_e", "x_p");
    xlabel("Time (t)");
    ylabel("Distance");
    title("Distance to Boundary");
    ax = gca;
    ax.FontSize = fontSize;           % Tick label font size
    ax.LabelFontSizeMultiplier = 1;
    ax.TitleFontSizeMultiplier = 1;
    if save == 1
        set(figX, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
        savefig(figX,sprintf("%s/d_X.fig",dirname));
        print(figX, sprintf("%s/d_X.png", dirname), '-dpng', '-r300'); % Save with 300 DPI resolution
    end
    hold off;
    
    %Distance to unsafe zone
    figXu = figure('visible',vis);
    plot(vecnorm(Xe-Xp,2,2), 'b', 'LineWidth', 1.6);
    yline(Rp,'r--', 'LineWidth', 1.6);
    xlim([0 size(Xe,1)]);
    ylim([0 inf]);
    xlabel("Time (t)");
    ylabel("Distance");
    title("Distance Between x_e and x_p");
    ax = gca;
    ax.FontSize = fontSize;           % Tick label font size
    ax.LabelFontSizeMultiplier = 1;
    ax.TitleFontSizeMultiplier = 1;
    if save == 1
        set(figXu, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
        savefig(figXu,sprintf("%s/d_Xa.fig",dirname));
        print(figXu, sprintf("%s/d_Xa.png", dirname), '-dpng', '-r300'); % Save with 300 DPI resolution
    end
    
    % Distance to the target center
    figXt = figure('visible',vis);
    plot(vecnorm(Xe-xr0',2,2), 'b', 'LineWidth', 1.2);
    yline(Rr,'g--', 'LineWidth', 1.2);
    xlim([0 size(Xe,1)]);
    ylim([0 inf]);
    xlabel("Time (t)");
    ylabel("Distance");
    title("Distance Between x_e and x_r");
    ax = gca;
    ax.FontSize = fontSize;            % Tick label font size
    ax.LabelFontSizeMultiplier = 1;
    ax.TitleFontSizeMultiplier = 1;
    if save == 1
        set(figXt, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
        savefig(figXt,sprintf("%s/d_Xr.fig",dirname));
        print(figXt, sprintf("%s/d_Xr.png", dirname), '-dpng', '-r300'); % Save with 300 DPI resolution
    end

    toc

end