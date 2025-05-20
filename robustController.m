function out = robustController(Init, Par, f, ge, vars, options)
% Solve the optimization problem to get the full-information controller
% Density with Lyapunov denominator, 1-step: Find rho and u simultaneously
%   Init:       System initial conditions
%   Par:        System parameters 
%   f,g:        System dynamics
%   vars:       SDP Variables for f and g
%   options:    Solver options

% Get initial conditions
xe0 = Init.xe0;
xp0 = Init.xp0;
xr0 = Init.xr0;
Rie = Init.Rie;
Rip = Init.Rip;
Rp = Init.Rp;
Rr = Init.Rr;
R = Init.R;

% Get system parameters
n   = Par.n;
ne = Par.ne;
np = Par.np;
n_u = Par.n_u;
uMax = Par.uMax;
wMax = Par.wMax;
drho = Par.drho;
dpsi = Par.dpsi;
alpha = Par.alphaV;
ds = Par.ds;
dl = Par.dl;
du = Par.du;
tolXe = Par.tolXe;
tolXp = Par.tolXp;
tolu = Par.tolu;

% Get the polynomial dynamics (only fe matters)
fe = f(1:ne);

% Create a Polynomial for rho = rho_bar / V^alpha
[rho_bar, c_rho_bar] = polynomial(vars, drho, 0);

% Create Polynomials for psi =psi_bar / V^alpha
psi_bar = [];
c_psi_bar = [];
% Loop through and generate corresponding polynomials
for i = 1:n_u
    eval(sprintf("[psi_bar%d, c_psi_bar%d] = polynomial(vars, dpsi, 0)", i, i));
    eval(sprintf("psi_bar = [psi_bar; psi_bar%d]",i));
    eval(sprintf("c_psi_bar = [c_psi_bar, c_psi_bar%d]",i));
end

% Define N, e here
N = [eye(np);-eye(np)];
e = wMax*ones(2*np,1);

% Get max degree
dmax = max(dpsi, drho) - 1;
% Define y here
y = [];
cy = [];
% Loop through and generate corresponding polynomials
for i = 1:size(N,1)
    eval(sprintf("[y%d, cy%d] = polynomial(vars, dmax, 0)", i, i));
    eval(sprintf("y = [y; y%d]",i));
    eval(sprintf("cy = [cy, cy%d]",i));
end

% SOS polynomials (sigma)
[sie, cie] = polynomial(vars, ds, 0);
[sip, cip] = polynomial(vars, ds, 0);
[su, cu] = polynomial(vars, ds, 0);
[ste, cte] = polynomial(vars, ds, 0);
[ste1, cte1] = polynomial(vars, ds, 0);
[lXe, dXe] = polynomial(vars, dl, 0);
[lXe1, dXe1] = polynomial(vars, dl, 0);
[lXp, dXp] = polynomial(vars, dl, 0);
[sXe, cXe] = polynomial(vars, ds, 0);
[sXp, cXp] = polynomial(vars, ds, 0);

[sNXe1, cNXe1] = polynomial(vars, ds, 0);
[sNXe2, cNXe2] = polynomial(vars, ds, 0);
[sNXp1, cNXp1] = polynomial(vars, ds, 0);
[sNXp2, cNXp2] = polynomial(vars, ds, 0);

[sdXe1, cdXe1] = polynomial(vars, ds, 0);
[sdXp1, cdXp1] = polynomial(vars, ds, 0);

[suXe1a, cuXe1a] = polynomial(vars, du, 0);
[suXe2a, cuXe2a] = polynomial(vars, du, 0);
[suXp1a, cuXp1a] = polynomial(vars, du, 0);
[suXp2a, cuXp2a] = polynomial(vars, du, 0);
[suu1a, cuu1a] = polynomial(vars, du, 0);
[suu2a, cuu2a] = polynomial(vars, du, 0);

[suXe1b, cuXe1b] = polynomial(vars, du, 0);
[suXe2b, cuXe2b] = polynomial(vars, du, 0);
[suXp1b, cuXp1b] = polynomial(vars, du, 0);
[suXp2b, cuXp2b] = polynomial(vars, du, 0);
[suu1b, cuu1b] = polynomial(vars, du, 0);
[suu2b, cuu2b] = polynomial(vars, du, 0);


if length(vars) == 4
    states = vars;
elseif length(vars) == 5
    states = vars(1:4);
    states(3) = vars(4);
    states(4) = vars(5);
    theta = vars(3);
end

% Define the safety conditions
hie = (states(1) - xe0(1))^2 + (states(2) - xe0(2))^2 - Rie^2;          % Initial condition for evader
hip = (states(3) - xp0(1))^2 + (states(4) - xp0(2))^2 - Rip^2;          % Initial condition for evader
hXe = states(1)^2 + states(2)^2 - R^2;                                  % Restriction over the whole set X
hXp = states(3)^2 + states(4)^2 - R^2;                                  % Initial condition for pursuer
hu = (states(1) - states(3))^2 + (states(2) - states(4))^2 - Rp^2;      % Distance constraint (unsafe zone) between evader and pursuer
hte = (states(1)-xr0(1))^2 + (states(2)-xr0(2))^2 - Rr^2;               % Target zone condition for the evader
V = (states(1) - 1*xr0(1))^2 + (states(2) - 1*xr0(2))^2;                % Lyapunov function for the denominator of rho and psi

% SOS constraints for safety and control
F = [];

% 1. Initial set conditions:    rho >= 0 in Xi
F = [F, sos(rho_bar + sie*hie + sip*hip)];

% 2. Unsafe zone conditions:    rho < 0 in Xu and rho = 0 in boundary(X)
F = [F, sos(-rho_bar - lXe*hXe - ste*hte)];
F = [F, sos(-rho_bar - lXp*hXp)];
F = [F, sos(-rho_bar + su*hu + sXe*hXe + sXp*hXp)];

% 3. Make y sos
F = [F, sos(y)];

% 4. Get r(x) = -grad_p(rho) and set d1 = y*N - r = 0
% -> div_p(V) = 0 since V is not a function of xp
r = -1*jacobian(rho_bar, vars(ne+1:n));
d1 = y'*N-r;
F = [F, sos(d1 + sNXe1*hXe + sNXp1*hXp)];
F = [F, sos(-d1 + sNXe2*hXe + sNXp2*hXp)];

% 5. Calculate the divergence of psi wrt xe
fprintf("\nCalculating divergence of psi...\n")
div_psi_bar = 0; % Initialize the divergence
for i = 1:ne
    div_psi_bar = div_psi_bar + jacobian(rho_bar*fe + ge*psi_bar, vars(i));
end
d2 = V*y'*e - V*div_psi_bar + alpha*dot(jacobian(V,vars(1:ne)), rho_bar*fe + ge*psi_bar);
F = [F, sos(-d2 + sdXe1*hXe + sdXp1*hXp)];

% 6. Target conditions:         rho > 0 in Xt
F = [F, sos(rho_bar + ste1*hte - lXe1*hXe)];

% 7. Bound u or theta
if Par.inp_norm == 'inf'
    disp('The inf norm of input is bounded.')
    F = [F, sos(uMax*rho_bar*ones(n_u,1) - psi_bar ...
            + suXe1a*(hXe + tolXe) + suXp1a*(hXp + tolXp) - suu1a*(hu - tolu))]; ...
    F = [F, sos(uMax*rho_bar*ones(n_u,1) + psi_bar ...
            + suXe2a*(hXe + tolXe) + suXp2a*(hXp + tolXp) - suu2a*(hu - tolu))]; ...
elseif Par.inp_norm == 1
    disp('The l1 norm of input is bounded.')
    F = [F, sos(uMax*rho_bar + [-1, -1]*psi_bar ...
            + suXe1a*(hXe + tolXe) + suXp1a*(hXp + tolXp) - suu1a*(hu - tolu))];
    F = [F, sos(uMax*rho_bar + [1, 1]*psi_bar ...
            + suXe2a*(hXe + tolXe) + suXp2a*(hXp + tolXp) - suu2a*(hu - tolu))];
    F = [F, sos(uMax*rho_bar + [1, -1]*psi_bar ...
            + suXe1b*(hXe + tolXe) + suXp1b*(hXp + tolXp) - suu1b*(hu - tolu))];
    F = [F, sos(uMax*rho_bar + [-1, 1]*psi_bar ...
            + suXe2b*(hXe + tolXe) + suXp2b*(hXp + tolXp) - suu2b*(hu - tolu))];
else
    disp('The input norm is invalid.')
end

% 8. Ensure slack variables are SOS
F = [F, sos(sie), sos(sip), sos(su), sos(ste), sos(ste1), sos(sXe), sos(sXp),  ...
        sos(sNXe1), sos(sNXe2), ...
        sos(sNXp1), sos(sNXp2), ...
        sos(sdXe1), sos(sdXp1), ...
        sos(suXe1a), sos(suXp1a), ...
        sos(suXe2a), sos(suXp2a), ...
        sos(suu1a), sos(suu2a)];

% Set up the objective function (empty since we are solving a feasibility problem)
objective = [];

% Solve the SOS problem using solvesos to find rho
fprintf("\nSolving for rho...\n")

[sol, v, Q] = solvesos(F, objective, options, [c_rho_bar; c_psi_bar(:); ...
                            cy(:); ...
                            cie; cip; cu; cte; cte1; dXe; dXe1; dXp; cXe; cXp;  ...
                            cNXe1; cNXe2; ...
                            cNXp1; cNXp2; ...
                            cdXe1; cdXp1; ...
                            cuXe1a; cuXp1a; ...
                            cuXe2a; cuXp2a; ...
                            cuu1a; cuu2a; ...
                            cuXe1b; cuXp1b; ...
                            cuXe2b; cuXp2b; ...
                            cuu1b; cuu2b]);

% Check if the problem was solved successfully
if sol.problem == 0
    disp('Feasible solution found for rho.');
    disp(sol.info);
else
    fprintf('\n!!! NO FEASIBLE SOLUTION FOUND FOR RHO!!!\n');
    disp(sol.info);
end

out.method = "rcDV";            % 'r'obust 'c'ontroller with 'D'ensity and 'V' (Lyapunov)
out.c_rho = value(c_rho_bar);
out.c_psi = value(c_psi_bar);
out.rho = rho_bar;
out.psi = psi_bar;
out.sol = sol;
out.F = F;
out.Q = Q;
out.v = v;

end
