% 
% Matt Werner
% Chase Leuchtmann
% 
% Solve Poisson's equation using a Green's function representation of the
% solution which is obtained via numerical integration.
% 
% April 10, 2021
% 

% Obtain the perturber mass using numerical integration
M = integral2(@matterDensity, ...
    -1, 1, ...
    -1, 1, ...
    'AbsTol', 1e-10, ...
    'RelTol', 1e-10, ...
    'Method', 'iterated');

% Define the number of points along each dimension at which the solution is
% to be evaluated
N = 50;

% Define the coordinates at which the solution is to be evaluated
rvec = linspace(0,1,N);
tvec = linspace(0,2*pi,N);

% Allocate memory for the solution
PhiMatlabG = nan(N);

% Begin evaluating the solution via Green's function
i = 1;
for r = rvec
    % Set the t coordinate index
    j = 1;
    for t = tvec
        
        % Define integrands
        G4pirho = @(rp, tp) poissonGDisk(r, t, rp, tp).*4.*pi.*matterDensity(rp.*cos(tp), rp.*sin(tp)).*rp;
        dGdn_NewtonianPotential = @(tp) M*poissonGDiskNormalD(r, t, tp);

        % Evaluate the solution at this position (r,t) corresponding with
        % indices (i,j)
        PhiMatlabG(i,j) = integral2(G4pirho, 0, 1-1e7*eps, 0, 2*pi, ...
                            'AbsTol', 1e-10, ...
                            'RelTol', 1e-10, ...
                            'Method', 'iterated') ...
                        - integral(dGdn_NewtonianPotential, 0, 2*pi, ...
                            'AbsTol', 1e-10, ...
                            'RelTol', 1e-10);
                        
        % Enforce that the potential can't actually reach 0 in finite space
        if (PhiMatlabG(i,j) == 0)
            PhiMatlabG(i,j) = nan;
        end
                        
        % Increment the t coordinate index
        j = j + 1;
    end
    % Increment the r coordinate index
    i = i + 1;
end

% Plot solution obtained from Green's function
[T, R] = meshgrid(tvec, rvec);
mesh(R.*cos(T), R.*sin(T), PhiMatlabG, 'EdgeColor', 'interp');
colormap copper
title("Gravitational Potential from Green's Function", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$y$", 'interpreter', 'latex')
zlabel("Gravitational Potential $\Phi$", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')

%% Calculate geodesics
departureAngles = [262, 285];
Geodesics = geodesics(X, Y, dx_htt, dy_htt, departureAngles);

%% Plot
% Plot domain boundary
plotGeodesics(Geodesics, departureAngles)