function f = fcoeffunction(location, state)
% 
% Nonhomogeneous function
% 

N = 1; % Number of equations
nr = length(location.x); % Number of columns
f = zeros(N, nr); % Allocate f

% Newton's gravitational constant in natural units
G = 1;

% Matter density function
rho = matterDensity(location.x, location.y);
% rho = 1;

% Now the particular functional form of f
% f(1,:) = location.x - location.y + state.u(1,:);
f(1,:) = 4*pi*G*rho;
end