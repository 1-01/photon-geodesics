% 
% Matt Werner
% Chase Leuchtmann
% 
% Numerical integration implementation of the solution of Poisson's
% equation using the Green's function for the unit disk
% 
% April 10, 2021
% 

% Specify the domain
rvec = linspace(0, 1, 50);
tvec = linspace(0, 2*pi, 50);

% Allocate memory for the solution
Phi = nan(numel(rvec), numel(tvec));

% Begin procedure for obtaining the solution using Green's function on the
% unit disk
i = 1; % Start radius counter
for r = rvec
    j = 1;
    for t = tvec
        % Fill in the value of the gravitational potential at this (r, t)
        % located at index (i, j)
        Phi(i,j) = interiorIntegral(r, t, 100, 100) - boundaryIntegral(r, t, 100);
        
        % Enforce that the potential can't actually reach 0 in finite space
        if (Phi(i,j) == 0)
            Phi(i,j) = nan;
        end
        
        % Increment theta counter
        j = j + 1;
    end
    % Increment radius counter
    i = i + 1;
end

%% Calculate geodesics
% Interpolate the solution of Poisson's equation since the original grid
% was made in polar coordinates (meaning the Cartesian grid is nonuniform)

% Meshing the solution coordinates (polar)
[T, R] = meshgrid(tvec, rvec);

% Create a uniform mesh (Cartesian)
x = linspace(-1,1,500);
y = x;
[X, Y] = meshgrid(x, y);

% Calculate the polar coordinates in the uniform Cartesian mesh
cartR = sqrt(X.^2 + Y.^2);
cartT = atan2(Y, X);
cartT(cartT < 0) = cartT(cartT < 0) + 2*pi;

% Interpolate the obtained solution (polar coordinates) for a uniform
% mesh in Cartesian coordinates
cartPhi = interp2(T, R, Phi, cartT, cartR);

% Compute the gradient of the perturbing metric
[dx_htt, dy_htt] = gradMetric(x, y, cartPhi);

% Plot the interpolated solution obtained from using a Green's function
figure
mesh(X, Y, cartPhi, 'EdgeColor', 'interp')
colormap copper
title("Gravitational Potential from Green's Function", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$y$", 'interpreter', 'latex')
zlabel("Gravitational Potential $\Phi$", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')


% Calculate geodesics
initialPosition = [0, 0.95];
departureAngles = [262, 285];
Geodesics = geodesics(X, Y, dx_htt, dy_htt, departureAngles, initialPosition);

%% Plot Geodesics
plotGeodesics(Geodesics, departureAngles)