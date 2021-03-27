% Applied Numerical Methods
% GR Photon Trajectories (Poisson Eq.)
% March 25, 2021
% Chase Leuchtmann, Matt Werner
clearvars, clc, close all
global p

PDE = createpde();
R1 = [3,4,-1,1,1,-1,1,1,-1,-1]';
gm = R1;
sf = 'R1';

% Create the geometry
ns = char('R1');
ns = ns';
g = decsg(gm, sf, ns);

% Include the geometry in the model and plot it
geometryFromEdges(PDE, g);
pdegplot(PDE, 'EdgeLabels', 'on')
axis equal
xlim([-1.1, 1.1])

% Create the mesh
p = 0.01;
PDEmesh = generateMesh(PDE, 'Hmax', p/2, 'GeometricOrder', 'quadratic');

% BCs
applyBoundaryCondition(PDE, 'dirichlet', 'Edge', 1:4, 'u', 0);

% Equation coefficients
m = 0;
d = 0;
c = -1;
a = 0;
% Specify functional form of nonhomogeneous side f
specifyCoefficients(PDE, 'm',0, 'd',0, 'c',c, 'a',a, 'f', @fcoeffunction);
% Solve the PDE
sol = solvepde(PDE); % for stationary problems

% Plot the mesh
figure
pdemesh(PDEmesh);
% Plot the result
figure
pdeplot(PDE, 'XYData', sol.NodalSolution)

close all
figure
pdeplot(PDE, 'XYData', sol.NodalSolution, 'ZData', sol.NodalSolution)
grid on


% Nonhomogeneous function
function f = fcoeffunction(location, state)
global p

N = 1; % Number of equations
nr = length(location.x); % Number of columns
f = zeros(N, nr); % Allocate f

% Newton's gravitational constant in natural units
G = 1;

% Location of the perturber
x0 = 0.3;
y0 = -0.6;
% Matter density function
% rho = exp((-(location.x).^2 - (location.y).^2)/(1/100)) + ...
%       3*exp((-(location.x - x0).^2 - (location.y - y0).^2)/(1/100)) + ...
%       3*exp((-(location.x - 2*x0).^2 - (location.y + y0/2).^2)/(1/100));

deltaxy = rectangularPulse(-p,p,location.x).*rectangularPulse(-p,p,location.y)/p;

rho = deltaxy;



% Now the particular functional form of f
% f(1,:) = location.x - location.y + state.u(1,:);
f(1,:) = 4*pi*G*rho;
end