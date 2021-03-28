function bcMatrix = rcoeffunction(location, state)


G = 1;

% Mass
global circleRadius
M = integral2(@matterDensity, ...
    -circleRadius, circleRadius, ...
    -circleRadius, circleRadius, ...
    'AbsTol', 1e-10, 'RelTol', 1e-10, 'Method', 'iterated');

% Radius from the source
r = sqrt(location.x.^2 + location.y.^2);

% Use Newtonian potential on the boundary
bcMatrix = -G*M./r;
end