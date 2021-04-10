function dGdn = poissonGDiskNormalD(r, t, tp)
% 
% Evaluate the normal derivative of the Green's function for Poisson's
% equation on the unit disk in polar coordinates
% 

dGdn = 1/(2*pi)*(1 - r.^2)./(1 - 2*r.*cos(t - tp) + r.^2);