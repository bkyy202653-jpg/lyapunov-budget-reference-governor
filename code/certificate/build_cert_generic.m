function Gt = build_cert_generic(P, vn, cap, iters, fname)
% Certified piecewise-linear Gamma_c (Lemma 3 of the paper) for parameters P on nodes vn; log-bisection per cell.
% Uses parfor over cells when a pool is available.
N = numel(vn) - 1; b = zeros(1, N); und = zeros(1, N);
parfor j = 1:N
    lo = log(1e-5); hi = log(0.3); ok0 = certify_level(vn(j), vn(j+1), 1e-5, P, cap);
    if ok0 ~= 1, b(j) = 0; continue; end
    u = 0;
    for it = 1:iters
        m = (lo + hi)/2; ok = certify_level(vn(j), vn(j+1), exp(m), P, cap);
        u = u + (ok == -1);
        if ok == 1, lo = m; else, hi = m; end
    end
    b(j) = exp(lo); und(j) = u;
end
g = [b(1), min(b(1:end-1), b(2:end)), b(end)];
Gt.v = vn; Gt.G = g; Gt.Dmax = 1; Gt.b = b; Gt.LG = max(abs(diff(g))./diff(vn)); Gt.und = sum(und);
Gt.note = 'certify_level.m, log-bisection, undecided = not certified, relative inflation 1e-12';
save(fname, 'Gt', 'P');
end
