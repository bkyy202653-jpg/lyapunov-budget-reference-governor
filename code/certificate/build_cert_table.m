% Build the certified piecewise-linear threshold Gamma_c on I = [-0.90, 0.85] (Lemma 3 of the paper):
% certify one level b_j per full cell [v_j, v_j+1] by branch-and-bound interval arithmetic (certify_level.m),
% then g_j = min(b_{j-1}, b_j) at interior nodes and linear interpolation.
P = params_lbg(); P.c = [1.2; 0.8]; P.umax = 2.5; P.k1 = 1; P.k2 = 1; P.g = 20; P.margin = 0.02;
T = load('results/tight_lbgtune.mat'); Gs = T.bestGt;          % sampled table, only used as a bisection bracket
vn = round((-0.90:0.01:0.85)*100)/100; N = numel(vn) - 1;
b = zeros(1, N); nb = zeros(1, N); und = zeros(1, N);
tic;
for j = 1:N
    hi = 1.3*max(interp1(Gs.v, Gs.G, [vn(j) vn(j+1)])) + 1e-4; lo = 0;
    for it = 1:16
        m = (lo + hi)/2; [ok, info] = certify_level(vn(j), vn(j+1), m, P, 1e6);
        nb(j) = nb(j) + info.nproc; und(j) = und(j) + (ok == -1);
        if ok == 1, lo = m; else, hi = m; end
    end
    b(j) = lo;
    if mod(j, 25) == 0, fprintf('cell %d/%d  v=[%.2f %.2f]  b=%.6f  (%.0fs)\n', j, N, vn(j), vn(j+1), b(j), toc); end
end
g = [b(1), min(b(1:end-1), b(2:end)), b(end)];
Gt.v = vn; Gt.G = g; Gt.Dmax = 1; Gt.b = b; Gt.LG = max(abs(diff(g))./diff(vn));
Gt.note = 'certified by certify_level.m (interval B&B, relative inflation 1e-12, no directed rounding)';
save results/Gamma_cert.mat Gt nb und
fprintf('done in %.0fs: min b %.6f, max b %.6f, L_Gamma %.5f, undecided tests %d\n', toc, min(b), max(b), Gt.LG, sum(und));
gs = interp1(Gs.v, Gs.G, vn);
for v = [-0.9 -0.6 -0.3 0 0.3 0.5 0.7 0.85]
    fprintf('v=%5.2f  certified %.5f  sampled %.5f  ratio %.3f\n', v, interp1(vn, g, v), interp1(Gs.v, Gs.G, v), interp1(vn, g, v)/interp1(Gs.v, Gs.G, v));
end
