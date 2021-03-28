function rho = matterDensity(x, y)
% Location of the perturber
x0 = 0;
y0 = 0;
rho = 7*exp(-((x - 2*x0).^2 + (y + y0/2).^2)*500);
  
  
% p = 0.1;
% deltaxy = rectangularPulse(x0-p,x0+p,x).*rectangularPulse(y0-p,y0+p,y);
% 
% % Calculate the mass density in SI units (kg/m3)
% rho = perturberMatterDensity*deltaxy;
% rho = 2*deltaxy;
end