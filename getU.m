function [u, rho] = getU(Par, vars, V, opt_out)
% Create the control action u depending on the output type
%   Par:    Parameters of the system
%   vars:   Symbolic variables of the system
%   V:      Lyapunov function for the denominator
%   out:    Output of the optimization problem

% Get the method used to get the control action:
% - fic: Full Information Controller
% -     ficD:  Density, 1-step
% -     ficD2: Density, 2-step
% -     ficV:  Density with Lyapunov denominator, 1-step
% -     ficB:  Barrier, 1-step
% -     ficB2: Barrier, 2-step
% - rc:  Robust Controller (to predator action)
% -     rcD:   Density, 2-step, robust
% -     rcDV:  Density with Lyapunov denominator, 1-step, robust
method = opt_out.method;

% Get parameters
drho = Par.drho;
dpsi = Par.dpsi;
du = Par.du;
n_u = Par.n_u;
alpha = Par.alphaV;

% Get rho 
c_rho = opt_out.c_rho;
rho = c_rho'*monomials(vars, 0:drho);

if strcmp(method,"ficDV") || strcmp(method,"rcDV")
    % --- u = psi_n/rho_n ---
    % Get psi
    c_psi = opt_out.c_psi;
    c_psi(isnan(c_psi)) = 0;
    psi = c_psi'*monomials(vars, 0:dpsi);
    % Get u = psi_n/rho_n
    u = psi/rho;
    % Divide by Lyapunov functions
    rho = rho / V^alpha;
    psi = psi / V^alpha;
    % --- ---
    
elseif strcmp(method,"ficD")
    % --- u = psi/rho ---
    % Get psi
    c_psi = opt_out.c_psi;
    c_psi(isnan(c_psi)) = 0;
    psi = c_psi'*monomials(vars, 0:dpsi);
    % Get u = psi/rho
    u = psi./rho;

elseif strcmp(method,"ficB") || strcmp(method,"ficB2") || strcmp(method,"ficD2") || strcmp(method,"rcD")
    % --- u polynomial ---
    % Get the coefficients of u and create it as a symbolic vector
    c_u = opt_out.c_u;
    c_u(isnan(c_u)) = 0;
    u = c_u'*monomials(vars, 0:du);

else
    error("\nThe method to find the controller is not defined. Please check.\n")
    
end

end