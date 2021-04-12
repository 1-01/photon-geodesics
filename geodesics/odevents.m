function [value,isterminal,direction] = odevents(t, u)
    global circleRadius
    % Add an event when we reach the edge of the domain
    value = u(2)^2 + u(3)^2 - 1;
    isterminal = 1;
    direction = 0;
    
    if (isnan(value))
        value = 0;
    end
end