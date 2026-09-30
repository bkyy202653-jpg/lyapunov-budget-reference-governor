% Are the kappa = 500 violations an integration artefact? Re-run failing runs with smaller steps.
P = params_lbg(); P.c = [1.2; 0.8]; P.umax = 2.5; P.k1 = 1; P.k2 = 1; P.g = 20;
C = load('results/Gamma_cert.mat'); Gt = C.Gt; E = load('results/envelope_cert_settle.mat'); X0set = E.X0set;
a = linspace(0, 2*pi, 9); a(end) = []; THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
for r = E.Rset
    Pr = P; Pr.r = r; o = simulate(Pr, THb, X0b, 'LBG', 500, Gt);
    bad = find(o.maxviol > 1e-6 | ~isfinite(o.ts));
    if isempty(bad), fprintf('r=%5.2f: no failure at dt=2e-3\n', r); continue; end
    for dt = [2e-4 5e-5]
        Q = Pr; Q.dt = dt; Q.T = 15; o2 = simulate(Q, THb(:,bad), X0b(:,bad), 'LBG', 500, Gt);
        fprintf('r=%5.2f: %2d failing runs at dt=2e-3 -> dt=%.0e: maxviol %.1e, unsettled %d, min Delta %.2e\n', ...
            r, numel(bad), dt, max(o2.maxviol), sum(~isfinite(o2.ts)), min(o2.minDelta));
    end
end
