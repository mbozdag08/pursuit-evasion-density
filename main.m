% Code for "Safe Control for Pursuit-Evasion with Density Functions"
% Mustafa Bozdag, Arya Honarpisheh
% 2025.05.20

close all;  clear;  clc;  yalmip('clear');  rng(8);

%% Initial Conditions 
% - Radii
Init.R = 4;                 % Radius of the entire set we should be staying in (around 0,0)
Init.Rie = 0.3;             % Radius of the initial set for the evader
Init.Rip = 0.3;             % Radius of the initial set for the pursuer
Init.Rp = 0.5;              % Radius of the unsafe zone around the pursuer
Init.Rr = 0.5;              % Radius of the circle that includes the target arc
% - Coordiantes
Init.xe0 = [-2;0];          % Initial coordinates of the evader
Init.xp0 = [-2;2];          % Initial coordinates of the pursuer
theta_r = 45;               % Initial coordinates of the target arc center
Init.xr0 = [Init.R*cos(theta_r/180*pi); Init.R*sin(theta_r/180*pi)];

%% Parameters

% - System parameters
Par.type = 1;       % System type (for optimization, only evader matters)
Par.ne   = 2;       % Num of states for the evador
Par.np   = 2;       % Num of states for the pursuer
Par.n 	 = Par.ne + Par.np;	    % Num of states
Par.n_u  = 2;       % Num of inputs (dimension of u)
Par.uMax = 0.015;    % Maximum control action for the evader
Par.wMax = 0.01;   % Maximum control action for the pursuer
% - Polynomial degrees
Par.drho = 10;	    % Degree of rho
Par.dpsi = 10;      % Degree of psi
Par.alphaV = 18;    % The exponent of Lyapunov function
Par.du   = 10;      % Degree of u
Par.ds   = 6;       % Degree of sos polynomial coeffs
Par.dl   = 6;       % Degree of non-sos polynomial coeffs
% - Tolerences for the optimization
Par.tolXe = 5e-2;
Par.tolXp = 5e-2;
Par.tolu = 8e-2;
Par.inp_norm = 'inf';

%% Optimization

% Define variables for the SDP
sdp_vars = sdpvar(Par.n, 1);
% Get the system
[f, ge] = getSystem(Par, sdp_vars);
% Relative gap termination tolerance used by the interior-point optimizer for conic problems
relgaptol = 1e-3;
% Solve the optimization problem
options = sdpsettings('solver', 'mosek', 'verbose', 1, 'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', relgaptol);

fprintf("\nOptimization...\n");
tic
opt_out = robustController(Init, Par, f, ge, sdp_vars, options);
toc

if opt_out.sol.problem ~= 0
    proceed = input('\nInfeasible problem. Do you want to proceed? (y/n): \n', 's');
    if lower(proceed) ~= 'y'
        fprintf('Stopping the execution.\n');
        return
    end
end

%% Simulation Setup

% Get system parameters
n = Par.n;
ne = Par.ne;
np = Par.np;
n_u = Par.n_u;
uMax = Par.uMax;

% Change the system type (save the type used for optimization)
typeOpt = Par.type;
% For simulation, both evader and pursuer matters
Par.type = 2;               
% Generate symbolic variables
vars = sym('x', [n, 1]);

% Get system
[f,ge] = getSystem(Par, Init, vars);

% Get initial conditions
xe0 = Init.xe0;
xp0 = Init.xp0;
xt0 = Init.xr0;
rie = Init.Rie;
rip = Init.Rip;
rp = Init.Rp;
rt = Init.Rr;
R = Init.R;

if length(vars) == 4
    states = vars;
elseif length(vars) == 5
    states = vars(1:4);
    states(3) = vars(4);
    states(4) = vars(5);
    theta = vars(3);
end

% Define the safety conditions
hie = (states(1) - xe0(1))^2 + (states(2) - xe0(2))^2 - rie^2;          % Initial condition for evader
hip = (states(3) - xp0(1))^2 + (states(4) - xp0(2))^2 - rip^2;          % Initial condition for evader
hXe = states(1)^2 + states(2)^2 - R^2;                                  % Restriction over the whole set X
hXp = states(3)^2 + states(4)^2 - R^2;                                  % Initial condition for pursuer
ha = (states(1) - states(3))^2 + (states(2) - states(4))^2 - rp^2;      % Distance constraint (unsafe zone) between evader and pursuer
hre = (states(1)-xt0(1))^2 + (states(2)-xt0(2))^2 - rt^2;               % Target zone condition for the evader
V = (states(1) - 1*xt0(1))^2 + (states(2) - 1*xt0(2))^2;                % Lyapunov function for the denominator of rho and psi

%% u and rho
fprintf("\nCalculating u and rho...\n");
tic
[u, rho] = getU(Par, vars, V, opt_out);
toc

%% ODE

% Initial state
x0 = [Init.xe0; Init.xp0];
% Convert u and f to functions for the ODE
% - To avoid matlabFunction problems when u(1) or u(2) is a constant or 0
u(1) = u(1)+vars(1)*1e-16; %   u(2) = u(2)+vars(1)*1e-16;
u_fun = matlabFunction(u, 'vars', {vars});
% Convert f to a symbolic expression (for system type 1 where f=[0;0;0;0])
f = sym(f); 
f_fun = matlabFunction(f, 'vars', {vars});

% Create functions for the stop event
hXe_fun = matlabFunction(hXe, 'Vars', {vars(1:ne)});
ht_fun = matlabFunction(hre, 'Vars', {vars(1:ne)});
hu_fun = matlabFunction(ha, 'Vars', {vars});

% Set the ODE options
reltol = 1e-3;
abstol = 1e-3;
ref = 8;
stopXt = 1;  % Stop the ODE when we reach the target
stopXu = 1;  % Stop the ODE when we go to the unsafe set
ode_opt = odeset('RelTol', reltol, 'AbsTol', abstol, 'Refine', ref, ...
                    'Events', @(t, x) stopEvent(x, hXe_fun, ht_fun, hu_fun, stopXt, stopXu, Par));

xdot = @(t,x)  f_fun(x) + [ge * u_fun(x); 0; 0];    % xdot for the ODE solver
T = 500;                                            % Final time
N = 50;                                             % Number of points
t_span = linspace(0, T, N);

fprintf("\nSolving the ODE...\n")
tic
% Call ode45 with event handling, specifying the time points to return
[t_sol, x_sol, te, xe, ie] = ode45(xdot, t_span, x0, ode_opt);
toc
% Check if the stop event was triggered
if ~isempty(te)
    % Event occurred, stop the loop
    fprintf('Event occurred at t = %.2f \n', te(1));
end

% Get seperate arrays for the coordinates of the evader and pursuer
Xe = x_sol(:,1:n/2);
Xp = x_sol(:,n/2+1:n);

% Get the u values
U = u_fun(x_sol(1:end-1,:)');

fprintf("Is Safe?:  \t%d\n", all(vecnorm(Xe-Xp,2,2)>=rp));
fprintf("Is in X?:  \t%d\n", all(vecnorm(Xe,2,2)<=R) || ~isempty(te));
fprintf("Is in Xt?: \t%d\n", ~isempty(te));

%% Plots

savefigs = 1;   % Save figs
vis = 'on';    % Show figs

% Define the folder name based on the current date and time
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
dirname = sprintf('./Figs/%s_SOS', timestamp);

% Create the folder
if ~exist(dirname, 'dir') & (savefigs==1)
    mkdir(dirname);
end

% Define the parameter information
parameters = {
    'method', opt_out.method;
    'xe0_1', xe0(1);
    'xe0_2', xe0(2);
    'xp0_1', xp0(1);
    'xp0_2', xp0(2);
    'theta_t', theta_r;
    'typeOpt', typeOpt;
    'typeSim', Par.type;
    'n', Par.n;
    'ne', Par.ne;
    'np', Par.np;
    'n_u', Par.n_u;
    'uMax', Par.uMax;
    'wMax', Par.wMax;
    'drho', Par.drho;
    'dpsi', Par.dpsi;
    'alphaV', Par.alphaV;
    'tolXe', Par.tolXe;
    'tolXp', Par.tolXp;
    'tolu', Par.tolu;
    'T', T;
    'N', N;
    'Velocity Norm', Par.inp_norm
};

% Save parameters in a text file within the folder
% Open the text file for writing
if savefigs == 1
    paramFile = fullfile(dirname, 'parameters.txt');
    fid = fopen(paramFile, 'w');
    % Loop through the parameters
    for i = 1:size(parameters, 1)
        paramName = parameters{i, 1};
        paramValue = parameters{i, 2};
        % Check if the parameter is a string or a number
        if ischar(paramValue) || isstring(paramValue)
            fprintf(fid, '%s: %s\n', paramName, paramValue);  % Use %s for strings
        else
            fprintf(fid, '%s: %g\n', paramName, paramValue);  % Use %g for numbers
        end
    end
    % Close the file
    fclose(fid);
end

% Solve the system of equations hXe = 0 and ht = 0 to get the target arc
xtArc = solve([hXe == 0, hre == 0], vars(1:ne));
% Create a function for the density function rho
rho_fun = matlabFunction(rho, 'vars',{vars});
% Create a grid of points for the square encapsulating the circle
grid_spacing = 0.1; % Distance between grid points
[Rx1, Rx2] = meshgrid(-2*R:grid_spacing:2*R, -2*R:grid_spacing:2*R);
% Choose the index of the timepoint for the static plot
idxs = round(linspace(1, length(t_sol), 4));
% Choose quiver density for the phase portrait
quiver_density = 20;

% Plot the environment for idxs
for idx = idxs
    plotEnv(idx, Init, Rx1, Rx2, Xe, Xp, xtArc, f_fun, ge, u_fun, rho_fun, quiver_density, dirname, vis, savefigs);
end

% Plot parameters
plotParams(Init, Xe, Xp, dirname, vis, savefigs);

% Animation
if savefigs==1
    frameRate = 5; % Set desired frame rate
    quality = 100; % Set video quality (0 to 100)
    animateEnv(t_sol, Init, Rx1, Rx2, Xe, Xp, xtArc, f_fun, ge, u_fun, rho_fun, quiver_density, frameRate, quality, dirname);
end
