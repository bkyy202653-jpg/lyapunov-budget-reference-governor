% Why does LBG stall on the arm at r = -0.9 although Gamma_c(r) = 2.9 F? Check s_inf vs Gamma_c(v_inf) and theta error.
P = params_arm(); F = P.R^2/(2*P.g); L = load('results/Gamma_cert_arm.mat'); Gt = L.Gt; A = load('results/arm.mat');
X0set = A.X0set; a = linspace(0, 2*pi, 9); a(end) = []; THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
for r = [-0.9 1.0 0.6]
    Pr = P; Pr.r = r; Pr.T = 60; o = simulate(Pr, THb, X0b, 'LBG', 20, Gt);
    ve = o.rec.v(end,:); se = o.rec.s(end,:); Ve = o.rec.V(end,:); Vz = 0; %#ok<NASGU>
    st = abs(ve - r) > 0.02;
    fprintf('r=%5.2f (T=60): stalled %d/%d; median v_end %.3f; s_end/F %.2f; Gam(v_end)/F %.2f; V_theta_end/F %.2f (V_end ~ V_theta since z->0)\n', ...
        r, sum(st), numel(st), median(ve(st)), median(se(st))/F, median(interp1(Gt.v, Gt.G, ve(st)))/F, median(Ve(st))/F);
    % which true parameters stall? friction component of theta - theta0
    d = THb(:, st) - P.th0;
    fprintf('        stalled runs: friction offset th1-th10 in [%.2f, %.2f], gravity offset in [%.2f, %.2f]\n', min(d(1,:)), max(d(1,:)), min(d(2,:)), max(d(2,:)));
end
