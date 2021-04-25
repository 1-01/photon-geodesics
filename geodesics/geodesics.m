function geodesic = geodesics(X, Y, dx_htt, dy_htt, ...
                              photonTheta, photonOrigin)
% 
% Calculate geodesics from the perturbed Einstein equation (in 2 spatial
% dimensions) in which the speed of light is 1.
% 

% Allocate memory for varargout
geodesic = cell(size(photonTheta));

% Specify ODE options
options = odeset('RelTol', 1e-6, 'MaxStep', 0.01, 'Refine', 8, 'Events', @odevents);

% Plot the geodesics
k = 1;
for theta = photonTheta
    % Determine the direction in which this photon is initially travelling.
    % The angular coordinate theta is measured relative to the positive
    % x-axis.
    v0 = [cosd(theta), sind(theta)];
    
    % Specify initial conditions
    x_initial = [0, photonOrigin, 0, 1, v0, 0]';
    
    % Solve the geodesic equation
    sol = ode45(@(l, u) odefun(l, u, X, Y, dx_htt, dy_htt), [0, 2], x_initial, options);
    
    % Unpack the solution as (x0, x1, x2, x3) = (t, x, y, z)
    xmu = sol.y(1:4,:)';
    
    % Questionable lines
    if (any(isnan(xmu)))
        xmu = xmu(1:end-1,:);
    end
    
    % Assign output
    geodesic{k} = xmu;
    
    % Increment counter
    k = k + 1;
end