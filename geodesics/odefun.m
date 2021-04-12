function f = odefun(l, u, X, Y, dx_htt, dy_htt)
    disp("lambda = " + l)
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