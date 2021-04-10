function I = interiorIntegral(r, t, rPartitions, tPartitions)
% 
% Evaluate the interior integral of the Green's function solution using
% 2-point Gauss quadrature
% 
% Np and Mp are the requested resolutions of the meshes for the dummy
% variables rp and tp.
% 

% Define the external partition width and begin at the absolute lower bound
tWidth = 2*pi/tPartitions;
tlowerBound = 0;

% Begin the integral that is to be accumlated
I = 0;

% Evaluate the double integral at this position (r, t)
for i = 1:tPartitions
    % Adjust the upper bound
    tupperBound = tlowerBound + tWidth;
    
    % Evaluate the outer integral(s) at this position (r, t) with t' mapped
    % to [-1, 1] using 2-point Gauss quadrature
    td = gqp(tlowerBound, tupperBound);
    I = I + 0.5*tWidth*(innerIntegral(r, t, rPartitions, td(1)) + innerIntegral(r, t, rPartitions, td(2)));
    
    % Move to the next interval
    tlowerBound = tupperBound;
end

%% Nested functions

    % Define the inner (/internal/iterated) integral which is to be fully
    % evaluated throughout every partition segment used to evaluate the
    % outer (/external/main) integral.
    function f = innerIntegral(r, t, rPartitions, td)
        % Define the internal partition width and begin at the absolute
        % lower bound
        rWidth = 1/rPartitions;
        rlowerBound = 0;
        
        % Begin the integral that is to be accumulated
        f = 0;
        for j = 1:rPartitions
            % Adjust the upper bound
            rupperBound = rlowerBound + rWidth;
            
            % Evaluate the integral at this position (r, t, td) with r'
            % mapped to [-1, 1] using 2-point Gauss quadrature
            rd = gqp(rlowerBound, rupperBound);
            f = f + 0.5*rWidth*(phi(r, t, td, rd(1)) + phi(r, t, td, rd(2)));
            
            % Move to the next interval
            rlowerBound = rupperBound;
        end
    end


    % Provide the expression for the Green's function multiplied by the
    % nonhomogeneous term making up Poisson's equation evaluated at the
    % points (r, t) = (rd, td).
    function f = phi(r, t, td, rd)
        f = log((r^2 - 2*r*rd*cos(t - td) + rd^2)/(1 - 2*r*rd*cos(t - td) + r^2*rd^2))*matterDensity(rd*cos(td), rd*sin(td))*rd;
    end
end