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
dr = 0.01;
dt = pi/10;

rvec = 0:dr:1;
tvec = 0:dt:2*pi;

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
        Phi(i,j) = MMSinteriorIntegral(r, t, 250, 250) + MMSboundaryIntegral(r, t, 250);
        
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

% Compute the gradient of the perturbing metric
[dx_htt, dy_htt] = gradMetric(x, y, cartPhi);

% Plot the solution obtained from using a Green's function on the mesh
figure
mesh(R.*cos(T), R.*sin(T), Phi)
% mesh(X, Y, cartPhi, 'EdgeColor', 'interp') % Interpolated solution
colormap copper
title("Gravitational Potential from Green's Function", 'interpreter', 'latex')
xlabel("$x^1$", 'interpreter', 'latex')
ylabel("$x^2$", 'interpreter', 'latex')
zlabel("Gravitational Potential $\Phi$", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')
view(315,15)%140,29


PhiSensor = interp2(T, R, Phi, deg2rad(45), 0.5*sqrt(2))

Boundary = @(x,y) -(x.^4 + y.^4);
PhiExact = Boundary(0.5, 0.5)

PhiNormalizedError = abs(PhiSensor-PhiExact)/abs(PhiExact)