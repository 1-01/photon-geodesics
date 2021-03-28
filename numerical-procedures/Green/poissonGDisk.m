function G = poissonGDisk(r, t, rp, tp)
% 
% Evaluate the Green's function on the unit disk for Poisson's equation in
% polar coordinates
% 

G = 1/(4*pi)*log((rp.^2 + r.^2 - 2*r.*rp.*cos(t - tp))./(1 + r.^2.*rp.^2 - 2*r.*rp.*cos(t - tp)));