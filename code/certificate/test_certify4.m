P = params_lbg(); P.c = [1.2; 0.8]; P.umax = 2.5; P.k1 = 1; P.k2 = 1; P.g = 20; P.margin = 0.02;
% (a) infeasibility witnesses: at these levels a feasible point violates a constraint, the certifier must NOT return 1
[ok1, i1] = certify_level(-0.9, -0.9, 0.001481, P, 3e6);
[ok2, i2] = certify_level(0.85, 0.85, 0.001508, P, 3e6);
fprintf('witness v=-0.90 gam=0.001481: ok=%d (boxes %d)\n', ok1, i1.nproc);
fprintf('witness v=+0.85 gam=0.001508: ok=%d (boxes %d)\n', ok2, i2.nproc);
% (b) bisection of the certified level vs a level where fmincon finds a violation
opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp'); rng(2);
for v = [-0.9 -0.3 0 0.3 0.7 0.85]
    lo = 0; hi = 0.2;
    for it = 1:22
        m = (lo + hi)/2; ok = certify_level(v, v, m, P, 1e6);
        if ok == 1, lo = m; else, hi = m; end
    end
    % falsify slightly above the certified level
    worst = 0;
    for trial = 1:40
        y0 = [0.05*randn(2,1); 0.2*randn(2,1); randn(2,1)]; y0(5:6) = y0(5:6)/max(1, norm(y0(5:6)));
        for obj = 1:3
            [y, fv] = fmincon(@(y) -cval(y, v, P, obj), y0, [], [], [], [], [], [], @(y) nlc(y, 1.02*lo, P), opts);
            if all(nlc(y, 1.02*lo, P) <= 1e-9), worst = max(worst, -fv); end
        end
    end
    fprintf('v=%5.2f certified level %.6f ; fmincon at 1.02x: worst/lim = %.4f\n', v, lo, worst);
end
function val = cval(y, v, P, which)
th = P.th0 + P.R*y(5:6); h = th - y(3:4); x1 = y(1) + v;
x2 = y(2) - P.k1*y(1) - h(1)*x1^2; u = prestab(x1, x2, h(1), h(2), v, P);
lim = (1 - P.margin)*[P.c(1); P.c(2); P.umax]; vals = [abs(x1); abs(x2); abs(u)]./lim; val = vals(which);
end
function [c, ceq] = nlc(y, gam, P)
c = [0.5*(y(1)^2 + y(2)^2) + 0.5*(y(3)^2 + y(4)^2)/P.g - gam; y(5)^2 + y(6)^2 - 1]; ceq = [];
end
