function [rho, M] = matterDensity(x, y)
% 
% Define the matter density of an object at the origin as a Gaussian
% distribution with "amplitude" a and extent s
% 
a = 4.6;
s = 500;
rho = a*exp(-(x.^2 + y.^2)*s);

% Provide the exact form of the mass
M = pi*a/s;