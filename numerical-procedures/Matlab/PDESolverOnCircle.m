% Applied Numerical Methods
% GR Photon Trajectories (Poisson Eq.)
% March 25, 2021
% Chase Leuchtmann, Matt Werner
clearvars, clc, close all

PDE = createpde();

% Create the geometry
global circleRadius
circleRadius = 1; % Circle radius
C1 = [1, 0, 0, circleRadius]';
gm = C1;
sf = 'C1';
ns = char('C1');
ns = ns';
g = decsg(gm, sf, ns);

% Include the geometry in the model and plot it
geometryFromEdges(PDE, g);
pdegplot(PDE, 'EdgeLabels', 'on')
axis equal

% Create the mesh
PDEmesh = generateMesh(PDE, 'Hmax', 0.04);
% Plot the mesh on the previous plot
hold on, pdemesh(PDEmesh), hold off;

% Equation coefficients
m = 0; d = 0; c = -1; a = 0;
% Specify functional form of the PDE coefficients
specifyCoefficients(PDE, 'm',0, 'd',0, 'c',c, 'a',a, 'f', @fcoeffunction);

% BCs
applyBoundaryCondition(PDE, 'dirichlet', 'Edge', 1:4, 'h', 1, 'r', @rcoeffunction, 'Vectorized', 'on');

%% Solve the PDE
Phi = solvepde(PDE); % for stationary problems

plots

%% Interpolate for gradient of Grav. Potential
x = linspace(-circleRadius, circleRadius, 2000);
y = linspace(-circleRadius, circleRadius, 2000)';
[X, Y] = meshgrid(x, y);
[dPhidx, dPhidy] = evaluateGradient(Phi, X, Y);

% Reshape
dPhidx = reshape(dPhidx, size(X));
dPhidy = reshape(dPhidy, size(Y));

% Convert gradient to that of htt
dx_htt = -2*dPhidx;
dy_htt = -2*dPhidy;

%% Calculate geodesics
options = odeset('AbsTol', 1e-6, 'MaxStep', 0.05, 'Events', @odevents);

c = 1;

figure(2)
% hold on
% plot3(0,0,0,'b.','MarkerSize',100)

figure(4)
% plot([0.1, 0.1, -0.1, -0.1, 0.1], [-0.1, 0.1, 0.1, -0.1, -0.1])
plot(0,0,'.k')
hold on
% Plot the domain's boundary
theta = linspace(0, 2*pi, 10000);
xcirc = circleRadius*cos(theta);
ycirc = circleRadius*sin(theta);
plot(xcirc, ycirc, '-k')
% Plot the geodesics (262, 285)
for theta = [262, 285]
    v0 = [c*cosd(theta), c*sind(theta)];
    x_initial = [0, 0, 0.95*circleRadius, 0, 1, v0, 0]';
    sol = ode45(@(t, u) odefun(t, u, X, Y, dx_htt, dy_htt), [0, 2*circleRadius], x_initial, options);
    xmu = sol.y(1:4,:)';
    
    % Questionable lines
    if (any(isnan(xmu)))
        xmu = xmu(1:end-1,:);
    end
    
    % Determine line color for current geodesic
    geodcolor = rand(1,3);
    
    % Plots
    figure(4)
    plot(xmu(:, 2), xmu(:, 3), 'Color', geodcolor)
    axis equal, grid on
    xlabel("$x^1$", 'interpreter', 'latex')
    ylabel("$x^2$", 'interpreter', 'latex')
    title("Photon Geodesics", 'interpreter', 'latex')
    
    % Plot on gravitational potential plot (figure 2)
    figure (2)
    hold on
    xmu23OnPotential = interpolateSolution(Phi, xmu(:, 2), xmu(:, 3));
    plot3(xmu(:, 2), xmu(:, 3), xmu23OnPotential, 'Color', geodcolor)
    xlabel("$x^1$", 'interpreter', 'latex')
    ylabel("$x^2$", 'interpreter', 'latex')
    zlabel("$\Phi$", 'interpreter', 'latex')
    title("Photon Geodesics on the Gravitational Potential", 'interpreter', 'latex')
    colormap pink
    view(-135, 37)
    
    pause(0.001)
end



%% Functions
function f = odefun(t, u, X, Y, dx_htt, dy_htt)
    disp("lambda = " + t)
    %
    % du/dt = f(t, u)
    %
    
    x = u(2);
    y = u(3);
    % Evaluate the perturbed metric gradient at this spacetime (t, x, y, z)
    dx_httval = interp2(X, Y, dx_htt, x, y);
    dy_httval = interp2(X, Y, dy_htt, x, y);
    
    % Commit derivatives to more common names
    dx1 = u(5);
    dx2 = u(6);
    dx3 = u(7);
    dx4 = u(8);
    
    
    % Evaluate the geodesic equation at this spacetime (t, x, y, z)
    f = zeros(8, 1);
    f(1:4) = u(5:8);
    f(5) = dx_httval*dx1*dx2 + dy_httval*dx1*dx3;
    f(6) = 0.5*dx_httval*dx1^2 + 0.5*dx_httval*dx2^2 + 0.5*dx_httval*dx3^2 + 0.5*dx_httval*dx4^2 + dy_httval*dx2*dx3;
    f(7) = 0.5*dy_httval*dx1^2 + 0.5*dy_httval*dx2^2 + 0.5*dy_httval*dx3^2 + 0.5*dy_httval*dx4^2 + dx_httval*dx2*dx3;
    f(8) = dx_httval*dx2*dx4 + dy_httval*dx3*dx4;
end

function [value,isterminal,direction] = odevents(t, u)
    global circleRadius
    % Add an event when we reach the edge of the domain
    value = u(2)^2 + u(3)^2 - circleRadius^2;
    isterminal = 1;
    direction = 0;
    
    if (isnan(value))
        value = 0;
    end
end