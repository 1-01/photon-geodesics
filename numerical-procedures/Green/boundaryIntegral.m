function I = boundaryIntegral(r, t, tPartitions)
% 
% Evaluate the exterior integral of the Green's function solution using
% 2-point Gauss quadrature
% 
% Mp is the requested resolution of the mesh for the dummy variable tp.
% 

% Compute the mass of the perturber according to the Gaussian distribution
% of matter density,
%                    2
%                -s r
%       rho = a e     ,
%  for which the mass is
%             __  a
%         M = || ---.
%                 s
[~, M] = matterDensity(nan,nan);

% Define the external partition width and begin at the absolute lower bound
thetaWidth = 2*pi/tPartitions;
lowerBound = 0;

% Begin the integral that is to be accumlated
I = 0;

% Evaluate the double integral at this position (r, t, r') on the boundary
% of the domain where r' = 1.
for i = 1:tPartitions
    % Adjust the upper bound
    upperBound = lowerBound + thetaWidth;
    
    % Evaluate the integral at this position (r, t, r') with t' mapped
    % to [-1, 1] using 2-point Gauss quadrature
    td = gqp(lowerBound, upperBound);
    I = I + 0.5*thetaWidth*(phi(r, t, td(1)) + phi(r, t, td(2)));
    
    % Move to the next interval
    lowerBound = upperBound;
end

%% Nested function

    % Provide the expression for the normal derivative of the Green's
    % function to the unit disk evaluated on the boundary (r' = 1).
    function f = phi(r, t, td)
        f = (M/2/pi)*(1 - r^2)/(1 - 2*r*cos(t - td) + r^2);
    end
end