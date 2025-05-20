function [value, isterminal, direction] = stopEvent(x, hXe_fun, hr_fun, ha_fun, stopXt, stopXu, Par)
% ODE stop events
%   x:          States
%   hXe_fun:    Big region X for evader
%   hr_fun:     Reach (target) region Xr for evader
%   ha_fun:     Avoid (unsafe) region Xa 
%   stopXt:     Stop if we reach the target
%   stopXu:     Stop if we reach the unsafe set

    ne = Par.ne;

    % Compute the current values of hXe, ht, and hu
    hXe_value = hXe_fun(x(1:ne));              % Evaluates hXe based on current state x
    ht_value = hr_fun(x(1:ne));                % Evaluates ht based on x(1:2)
    hu_value = ha_fun(x);                 % Evaluates hu based on current state x
    % Event: hXe = 0 and ht < 0
    % - Considering the signs, hXe*ht crosses 0 in -1 direction (decreasing)
    htComb = hXe_value*ht_value;
    % Create the value array for both events
    value = [htComb; hu_value];               % Combine the two event conditions

    % Set terminal conditions:
    isterminal = [stopXt; stopXu];            % Stop the integration for both events
    direction = [-1; -1];                     % No direction for the first event, decreasing for hu
end
