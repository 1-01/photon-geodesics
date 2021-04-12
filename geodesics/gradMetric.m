function [dx_htt, dy_htt] = gradMetric(x, y, Phi)
% 
% Evaluate the gradient of the gravitational potential obtained as a
% solution to Poisson's equation
% 

% Take gradients of the gravitational potential and convert into the form
% of the perturbed metric
[dx_htt, dy_htt] = gradient(Phi, x, y);

dx_htt = -2*dx_htt;
dy_htt = -2*dy_htt;