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
R = linspace(0, 1, 50);
T = linspace(0, 2*pi, 50);

% Allocate memory for the solution
Phi = nan(numel(R), numel(T));

% Begin procedure for obtaining the solution using Green's function on the
% unit disk
i = 1; % Start radius counter
for r = R
    j = 1;
    for t = T
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

% Plot solution obtained from Green's function
[T, R] = meshgrid(T, R);
mesh(R.*cos(T), R.*sin(T), Phi, 'EdgeColor', 'interp')
colormap copper
title("Gravitational Potential from Green's Function", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$y$", 'interpreter', 'latex')
zlabel("Gravitational Potential $\Phi$", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')