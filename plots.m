%% Plot the result
figure
pdeplot(PDE, 'XYData', Phi.NodalSolution, 'ZData', Phi.NodalSolution)
grid on

figure
rho = matterDensity(Phi.Mesh.Nodes(1,:), Phi.Mesh.Nodes(2,:));
plot3(Phi.Mesh.Nodes(1,:),  Phi.Mesh.Nodes(2,:), rho, '.')