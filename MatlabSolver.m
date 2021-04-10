% Applied Numerical Methods
% GR Photon Trajectories (Poisson Eq.)
% March 25, 2021
% Chase Leuchtmann, Matt Werner
clearvars, clc, close all

PDE = createpde();

% Create the geometry
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

