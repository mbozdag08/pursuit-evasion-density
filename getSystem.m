function [f,ge] = getSystem(Par, Init, vars)
% Choose between different types of dynamical systems
%   Par:    System parameters
%   vars:   State variables of the system

% Get the system dimension (num of states)
n = Par.n;
ne = Par.ne;
np = Par.np;
% Get the system type
type = Par.type;

switch type

    case 1
        fprintf("\nSystem Type 1: \n" + ...
                "-------------- \n" + ...
                "ne: %d,    np: %d \n" + ...
                "Evader:  Holonomic \n" + ...
                "Pursuer: Stationary \n\n", ne, np);

        f = [0; 0; 0; 0];
        ge = [1 0; 0 1];

    case 2
        fprintf("\nSystem Type 2: \n" + ...
                "-------------- \n" + ...
                "ne: %d,    np: %d \n" + ...
                "Evader:  Holonomic \n" + ...
                "Pursuer: Line-of-sight (constant w) \n\n", ne, np);

        % The distance between the pursuer and evader
        d = norm(vars(1:n/2) - vars(n/2+1:n), 2);
        % Maximum control action for the pursuer
        wMax = Par.wMax;

        f = [0; 0; 
            wMax*(vars(1)-vars(3)) / (d + 1e-8); 
            wMax*(vars(2)-vars(4)) / (d + 1e-8)];
        ge = [1 0; 0 1];

    case 3
        fprintf("\nSystem Type 3: \n" + ...
                "-------------- \n" + ...
                "ne: %d,    np: %d \n" + ...
                "Evader:  Holonomic \n" + ...
                "Pursuer: Get the Middle \n\n", ne, np);

        % The distance between the pursuer and evader
        xc = (Init.xr0 + vars(1:ne)) / 2;
        d = norm(vars(ne+1:end) - xc, 2);
        % Maximum control action for the pursuer
        wMax = Par.wMax;
        f = [0; 0; 
            wMax*(xc(1) - vars(3)) / (d + 1e-8); 
            wMax*(xc(2) - vars(4)) / (d + 1e-8)];
        ge = [1 0; 0 1];    

    case 4
        fprintf("\nSystem Type 4: \n" + ...
                "-------------- \n" + ...
                "ne: %d,    np: %d \n" + ...
                "Evader:  Van der Pol oscillator \n" + ...
                "Pursuer: Stationary \n\n", ne, np);
        
        mu = 0.5;
        f = [0.005*vars(2);
            0.005*(mu*(1-vars(1)^2)*vars(2) - vars(1)); 
            0; 
            0];
        ge = [1 0; 0 1];

    case 5
        fprintf("\nSystem Type 5: \n" + ...
                "-------------- \n" + ...
                "ne: %d,    np: %d \n" + ...
                "Evader:  Van der Pol oscillator \n" + ...
                "Pursuer: Line-of-sight (constant w) \n\n", ne, np);

        % The distance between the pursuer and evader
        d = norm(vars(1:n/2) - vars(n/2+1:n), 2);
        % Maximum control action for the pursuer
        wMax = Par.wMax;
    
        mu = 0.5;
        f = [0.005*vars(2); 
            0.005*(mu*(1-vars(1)^2)*vars(2) - vars(1)); 
            wMax*(vars(1)-vars(3)) / (d + 1e-8); 
            wMax*(vars(2)-vars(4)) / (d + 1e-8)];
        ge = [1 0; 0 1];

end
 
end