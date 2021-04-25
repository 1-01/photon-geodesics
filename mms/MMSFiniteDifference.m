%% Method of Manufactured Solutions
% Matt Werner
% Chase Leuchtmann
% 
% 2nd-order central difference implementation of the solution of Poisson's
% equation
% 
% April 22, 2021
% 4.054688347183699e-04

% Clearing workspace, closing windows, and clearing command window
clear;
close all;
clc;
format long

% Mesh resolution, set h = dx = dy
h = 0.01/2;

% Geometry & boundary conditions, must contain unit disk
xmin = -1.1;
xmax = -1*xmin; % Symmetrical
ymin = xmin;
ymax = -1*ymin; % Symmetrical

% Boundary conditions
Boundary = @(x, y) -(x.^4+y.^4); % For boundary and beyond unit disk

% Mesh (h is given as input)
nx = (xmax-xmin)/h - 1; % Number of "interior nodes" in x-direction
ny = (ymax-ymin)/h - 1; % Number of "interior nodes" in y-direction
N = nx*ny; % Total number of interior nodes

% Allocate memory for A and b (sparse matrix & vector)
% A = sparse(N,N); 
% b = sparse(N,1);
A = spalloc(N,N, 5*N); 
b = spalloc(N,1, N);

reverseStr = '';
% Loop through all the interior nodes, and fill A & b.
for k=1:N
    % In the following, we look at the four neighbours of (i,j) that are
    % involved in the 5-point formula
    i = mod(k-1,nx) + 1; % get (i,j) from k.
    j = (k-i)/nx + 1; % get (i,j) from k.
    
    % Get current (x, y) coordinates from (i, j) indices
    currentX = xmin + i*h;
    currentY = ymin + j*h;
    
    % Fill the coefficients of central differencing scheme
    % - Phi(i-1,j) - Phi(i+1,j) - Phi(i,j-1) - Phi(i,j+1) + 4Phi(i,j) = -h^2*f(i,j)
    % - Phi(k-1)   - Phi(k+1)   - Phi(k-nx)  - Phi(k+nx)  + 4Phi(k)   = -h^2*f(k)
    
    % Check if the currentX and currentY coordinates are either an exterior
    % node or boundary node
    if currentX^2+currentY^2 >= 1 % Current node is exterior or on boundary
        A(k, :) = 0; % Make all other columns in row k equal to zero
        A(k,k) = 1; % Always 1 if exterior or on boundary
        b(k) = Boundary(currentX, currentY);
    else % Current node must be INSIDE of the unit disk
        % Fill the k-th row of A and the k-th element of b
        A(k,k) = 4; % The diagonal element, corresponding to T(i,j), is always 4
        b(k) = -h^2*(-12*currentX^2-12*currentY^2);
        
        % Check left neighbor
        leftX = xmin + (i-1)*h;
        leftY = ymin + j*h;
        if leftX^2 + leftY^2 >= 1
            b(k) = b(k) + Boundary(leftX, leftY);
        else
            A(k, k-1) = -1;
        end
        
        % Check right neighbor
        rightX = xmin + (i+1)*h;
        rightY = ymin + j*h;
        if rightX^2 + rightY^2 >= 1
            b(k) = b(k) + Boundary(rightX, rightY);
        else
            A(k, k+1) = -1;
        end
        
        % Check bottom neighbor
        bottomX = xmin + i*h;
        bottomY = ymin + (j-1)*h;
        if bottomX^2 + bottomY^2 >= 1
            b(k) = b(k) + Boundary(bottomX, bottomY);
        else
            A(k, k-nx) = -1;
        end
        
        % Check top neighbor
        topX = xmin + i*h;
        topY = ymin + (j+1)*h;
        if topX^2 + topY^2 >= 1
            b(k) = b(k) + Boundary(topX, topY);
        else
            A(k, k+nx) = -1;
        end
    end
   % Display the progress
   percentDone = 100*k/N;
   msg = sprintf('Percent done: %3.1f', percentDone);
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
% Solve AT = b for nodal temperature values
Phi = full(A\b);

x = xmin:h:xmax;
y = ymin:h:ymax;
[X, Y] = meshgrid(x, y);
Phivis = zeros(size(X)); % just for visualization
 
% Left boundary
PhiLeft = Boundary(xmin, ymin:h:ymax);
Phivis(1,:) = PhiLeft; 

% Right boundary
PhiRight = Boundary(xmax, ymin:h:ymax);
Phivis(nx+2,:) = PhiRight;

% Bottom boundary
PhiBottom = Boundary(xmin:h:xmax, ymin);
Phivis(:,1) = PhiBottom; 

% Top boundary
PhiTop = Boundary(xmin:h:xmax, ymax);
Phivis(:,ny+2) = PhiTop;

% Creating visualization matrix
for i = 2:nx+1
  for j = 2:ny+1
      Phivis(i,j) = Phi((j-2)*nx + i-1);
  end
end

% Plot gravitational potential Phi
figure
mesh(X, Y, Phivis', 'EdgeColor', 'interp')
colormap copper
title("Gravitational Potential from Central Finite Difference", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$y$", 'interpreter', 'latex')
zlabel("Gravitational Potential $\Phi$", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')

% Plot known manufactured solution
figure
mesh(X, Y, Boundary(X, Y), 'EdgeColor', 'interp')
colormap copper
title("Manufactured Solution", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$y$", 'interpreter', 'latex')
zlabel("Manufactured Solution ($\Phi$)", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')

% Comparing results from finite difference and manufactured solution
% Interpolate the obtained solution
PhiSensor = interp2(X, Y, Phivis', 0.5, 0.5)

PhiExact = Boundary(0.5, 0.5)

PhiNormalizedError = abs(PhiSensor-PhiExact)/abs(PhiExact)

% Plot both on subplot
figure
subplot(1, 2, 1)
mesh(X, Y, Phivis', 'EdgeColor', 'interp')
colormap copper
title("Gravitational Potential from Central Finite Difference", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$y$", 'interpreter', 'latex')
zlabel("Gravitational Potential $\Phi$", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')
subplot(1, 2, 2)
mesh(X, Y, Boundary(X, Y), 'EdgeColor', 'interp')
colormap copper
title("Manufactured Solution", 'interpreter', 'latex')
xlabel("$x$", 'interpreter', 'latex')
ylabel("$y$", 'interpreter', 'latex')
zlabel("Manufactured Solution ($\Phi$)", 'interpreter', 'latex')
h = colorbar;