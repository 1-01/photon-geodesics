function plotGeodesics(Geodesics, angles)
% 
% Plot the geodesics indicating the source, lens, and observer (where they
% intersect) using InterX from the Mathworks File Exchange
% 

% Find where curves intersect
P = InterX(Geodesics{1}(:,2:3)',Geodesics{2}(:,2:3)');
% Sort the intersections based on descending y so that the starting
% position is always first
[B, I] = sort(P(2,:), 'descend');
P(1,:) = P(1,I);
P(2,:) = B;


% Try to assign the last output
try
    % Try to index P for the last intersection
    [xO, yO] = deal(P(1,end), P(2,end));
catch e
    error("No intersection.")
end


% Plot domain boundary (circle)
fig = figure;
fig.Position = [2234, 575, 557, 529];
plot(0, 0, '.k', 'MarkerSize', 15, 'HandleVisibility', 'off')
hold on


% Add text
text(0.025,-0.05, '$L$', 'Interpreter', 'latex')
text(Geodesics{1}(1,2)-0.075, Geodesics{1}(1,3)-0.05, '$S$', 'Interpreter', 'latex')
text(xO-0.09, yO-0.05, '$O$', 'Interpreter', 'latex')
% Tangent lines
for k = 1:numel(Geodesics)
	[~, index] = min(abs(Geodesics{k}(:,3) - yO));
    m = (Geodesics{k}(index-1,3) - Geodesics{k}(index,3))/(Geodesics{k}(index-1,2) - Geodesics{k}(index,2));
    
    xtmp = linspace(-1,1,10000);
    ytmp = yO + m*(xtmp - xO);
    
    xtmp(xtmp.^2 + ytmp.^2 > 1 | ytmp < yO) = nan;
    ytmp(xtmp.^2 + ytmp.^2 > 1 | ytmp < yO) = nan;
    
    % Plot observer lines
    a = plot(xtmp, ytmp, 'k--', 'HandleVisibility', 'off');
    uistack(a, 'down')
    
    % Plot observer dots on the boundary
    [~, idx] = min(abs(xtmp.^2 + ytmp.^2 - 1));
    plot(xtmp(idx), ytmp(idx), '.k', 'MarkerSize', 15, 'HandleVisibility', 'off')
    if (k == 1)
        text(xtmp(idx)-0.09, ytmp(idx)+0.03, sprintf("$S_{%1.0f}$", k), 'Interpreter', 'latex')
    else
        text(xtmp(idx)+0.0, ytmp(idx)+0.05, sprintf("$S_{%1.0f}$", k), 'Interpreter', 'latex')
    end
    
end


t = linspace(0, 2, 10000);
plot(cospi(t), sinpi(t), '-k', 'HandleVisibility', 'off')
% Determine line color for current geodesic
geodcolor = [10, 28, 55; ...    blue
             80, 14, 42; ...    red
             17, 71, 03; ...    green
             85, 33, 10]/100; % orange
% Plot the geodesics
for k = 1:numel(Geodesics)
    plot(Geodesics{k}(:,2), Geodesics{k}(:,3), 'Color', geodcolor(k,:))
end
axis equal, grid on
xlabel("$x^1$", 'interpreter', 'latex')
ylabel("$x^2$", 'interpreter', 'latex')
title("Photon Geodesics", 'interpreter', 'latex')
xticks(-1:0.25:1);
yticks(-1:0.25:1);
xlim([-1.125, 1.125])
ylim([-1.125, 1.125])
hleg = legend(num2str(angles(1))+"$^\circ$", num2str(angles(2))+"$^\circ$", ...
                'interpreter', 'latex', ...
                'Position', [0.221110749611811,0.479914933026175,0.195579258501,0.096880908994909]);
title(hleg, "Departure Angle $\theta_0$", 'interpreter', 'latex')


% Add dots
plot(Geodesics{1}(1,2), Geodesics{1}(1,3), '.k', 'MarkerSize', 15, 'HandleVisibility', 'off')
plot(xO, yO, '.k', 'MarkerSize', 15, 'HandleVisibility', 'off')