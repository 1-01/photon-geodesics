
M = integral2(@matterDensity, ...
    -1, 1, ...
    -1, 1, ...
    'AbsTol', 1e-10, 'RelTol', 1e-10, 'Method', 'iterated');


N = 50;
PhiGF = nan(N);
R = linspace(0,1,N);
T = linspace(0,2*pi,N);
tic
i = 1;
for r = R
    j = 1;
    for t = T
        Gintfun = @(rp, tp) poissonGDisk(r, t, rp, tp).*4.*pi.*matterDensity(rp.*cos(tp), rp.*sin(tp)).*rp;
        dGdnintfun = @(tp) poissonGDiskNormalD(r, t, tp)*M;

        PhiGF(i,j) = integral2(Gintfun, 0, 1, 0, 2*pi, ...
                            'AbsTol', 1e-10, 'RelTol', 1e-10, 'Method', 'iterated') ...
                    - integral(dGdnintfun, 0, 2*pi, ...
                            'AbsTol', 1e-10, 'RelTol', 1e-10);
        j = j + 1;
    end
    i = i + 1;
end
toc

% Plot solution obtained from Green's function
[T, R] = meshgrid(T, R);
surf(R.*cos(T), R.*sin(T), PhiGF);

% interpolateSolution(Phi, r*cos(t), r*sin(t))