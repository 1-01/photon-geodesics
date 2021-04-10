function xd = gqp(a, b)
% 
% Define the Gauss quadrature points at which each integral in [-1, 1]
% should be evaluated
% 

xd = [1, -1]*(b - a)/2/sqrt(3) + (b + a)/2;
