% 
% Matt Werner
% Chase Leuchtmann
% 
% 2nd-order central difference implementation of the solution of Poisson's
% equation
% 
% April 16, 2021
% 

% Mesh resolution, set h = dx = dy
h = 0.025;

% Boolean for plotting 2D mesh (not useful)
plot2DMesh = false;

% Geometry & boundary conditions, must contain unit disk
xmin = -1.015;
xmax = -1*xmin; % Symmetrical
ymin = xmin;
ymax = -1*ymin; % Symmetrical

% Boundary conditions
G = 1; % G attains unit value
[~, M] = matterDensity(nan,nan); % Mass of perturbing body
NewtPotential = @(x, y) -G.*M./(sqrt(x.^2+y.^2)); % For boundary and beyond unit disk

% Mesh (h is given as input)
nx = floor((xmax-xmin)/h - 1); % Number of "interior nodes" in x-direction
ny = floor((ymax-ymin)/h - 1); % Number of "interior nodes" in y-direction
N = nx*ny; % Total number of interior nodes

% Allocate memory for A and b (sparse matrix & vector)
% A = sparse(N,N); 
% b = sparse(N,1);
A = spalloc(N,N, 5*N); 
b = spalloc(N,1, N);

reverseStr = '';
% Loop through all the interior nodes, and fill A & b.
for k=1:N
%     disp("k = " + num2str(k))
    % In the following, we look at the four neighbours of (i,j) that are
    % involved in the 5-point formula
    i = mod(k-1,nx) + 1; % get (i,j) from k.
    j = (k-i)/nx + 1; % get (i,j) from k.
    
    % Get current (x, y) coordinates from (i, j) indices
    currentX = xmin + i*h;
    currentY = ymin + j*h;
    
    % Fill the coefficients of central differencing scheme
    % - Phi(i-1,j) - Phi(i+1,j) - Phi(i,j-1) - Phi(i,j+1) + 4Phi(i,j) = h^2*f(i,j)
    % - Phi(k-1)   - Phi(k+1)   - Phi(k-nx)  - Phi(k+nx)  + 4Phi(k)   = h^2*f(k)
    
    % Check if the currentX and currentY coordinates are either an exterior
    % node or boundary node
    if currentX^2+currentY^2 >= 1 % Current node is exterior or on boundary
        A(k, :) = 0; % Make all other columns in row k equal to zero
        A(k,k) = 1; % Always 1 if exterior or on boundary
        b(k) = NewtPotential(currentX, currentY);
    else % Current node must be INSIDE of the unit disk
        % Fill the k-th row of A and the k-th element of b
        A(k,k) = 4; % The diagonal element, corresponding to T(i,j), is always 4
        b(k) = -4*pi*h^2*matterDensity(currentX, currentY);
        
        % Check left neighbor
        leftX = xmin + (i-1)*h;
        leftY = ymin + j*h;
        if leftX^2 + leftY^2 >= 1
            b(k) = b(k) + NewtPotential(leftX, leftY);
        else
            A(k, k-1) = -1;
        end
        
        % Check right neighbor
        rightX = xmin + (i+1)*h;
        rightY = ymin + j*h;
        if rightX^2 + rightY^2 >= 1
            b(k) = b(k) + NewtPotential(rightX, rightY);
        else
            A(k, k+1) = -1;
        end
        
        % Check bottom neighbor
        bottomX = xmin + i*h;
        bottomY = ymin + (j-1)*h;
        if bottomX^2 + bottomY^2 >= 1
            b(k) = b(k) + NewtPotential(bottomX, bottomY);
        else
            A(k, k-nx) = -1;
        end
        
        % Check top neighbor
        topX = xmin + i*h;
        topY = ymin + (j+1)*h;
        if topX^2 + topY^2 >= 1
            b(k) = b(k) + NewtPotential(topX, topY);
        else
            A(k, k+nx) = -1;
        end
    end
   % Display the progress
   percentDone = 100 * k / N;
   msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

% Solve AT = b for nodal temperature values
Phi = full(A\b);

% % Output Temperature at a sensor location (10 cm, 10 cm)
% i = (0.25*(xmax-xmin))/h;
% j = (0.25*(ymax-ymin))/h;
% Tsensor = T(nx*(j-1)+i);
% fprintf('h = %e, T(%d,%d) = %e.\n', h, i, j, Tsensor);

x = xmin:h:xmax;
y = ymin:h:ymax;
[X, Y] = meshgrid(x, y);
Phivis = zeros(size(X)); % just for visualization

% Left boundary
PhiLeft = NewtPotential(xmin, ymin:h:ymax);
Phivis(1,:) = PhiLeft; 

% Right boundary
PhiRight = NewtPotential(xmax, ymin:h:ymax);
Phivis(nx+2,:) = PhiRight;

% Bottom boundary
PhiBottom = NewtPotential(xmin:h:xmax, ymin);
Phivis(:,1) = PhiBottom; 

% Top boundary
PhiTop = NewtPotential(xmin:h:xmax, ymax);
Phivis(:,ny+2) = PhiTop;

% Creating visualization matrix
for i = 2:nx+1
  for j = 2:ny+1
      Phivis(i,j) = Phi((j-2)*nx + i-1);
  end
end

% Visualization
if plot2DMesh == true
    figure
    h = surf(X,Y,Phivis');
    colormap('copper');
    view(0,90);
    shading interp
    set(h,'edgecolor','k');
    cb = colorbar; 
    cb.Label.String = 'Gravitational Potential Phi';
    xlabel('$x$', 'interpreter', 'latex')
    ylabel('$y$', 'interpreter', 'latex')
    daspect([1 1 1]);
    % hold on
    % theta = 0:0.001:360;
    % plot(cosd(theta), sind(theta), 'k-');
end

% Plot gravitational potential Phi
figure
Phivis(X.^2 + Y.^2 > 1) = NaN;
mesh(X, Y, Phivis', 'EdgeColor', 'interp')
colormap copper
title("Gravitational Potential from Finite Difference", 'interpreter', 'latex')
xlabel("$x^1$", 'interpreter', 'latex')
ylabel("$x^2$", 'interpreter', 'latex')
zlabel("Gravitational Potential $\Phi$", 'interpreter', 'latex')
h = colorbar;
title(h, "$\Phi$", 'interpreter', 'latex')
view(315,15)
xticks(-1:0.5:1)
yticks(-1:0.5:1)
xlim([-1, 1])
xlim([-1, 1])
ylim([-1, 1])
zticks(-0.25:0.05:-0.05)

% Compute the gradient of the perturbing metric
[dx_htt, dy_htt] = gradMetric(x, y, Phivis);

% Calculate geodesics
initialPosition = [0, 0.95];
departureAngles = [262, 285];
Geodesics = geodesics(X, Y, dx_htt, dy_htt, departureAngles, initialPosition);

% Plot geodesics
plotGeodesics(Geodesics, departureAngles)

