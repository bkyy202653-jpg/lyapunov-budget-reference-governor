% Quick reproduction of the manipulator numbers: re-evaluates the DEPLOYED designs of run_arm_p1.m (friction known,
% p = 1) and the LBG of run_arm.m (friction also unknown, p = 2) without the tuning searches, and compares with the
% reference results. Needs in results/: Gamma_cert_arm_p1.mat, Gamma_cert_arm_p1_nom.mat, Gamma_cert_arm.mat,
% arm_p1.mat, arm.mat.
ok = true;
% ---------------- p = 1 (paper, Section "Single-link manipulator")
P = params_arm(); P.fixed1 = true; F = P.R^2/(2*P.g);
L = load('results/Gamma_cert_arm_p1.mat'); Gt = L.Gt; L = load('results/Gamma_cert_arm_p1_nom.mat'); Gn = L.Gt;
A = load('results/arm_p1.mat');
THt = [P.th0(1)*ones(1,5); P.th0(2) + P.R*[-1 -0.5 0 0.5 1]];
[I, J] = ndgrid(1:size(THt,2), 1:size(A.X0set,2)); THb = THt(:, I(:)); X0b = A.X0set(:, J(:));
Pn = P; Pn.g = 0; Pn.R = 0; Pn.k1 = A.bestN(1); Pn.k2 = A.bestN(1);
Pf = P; c = A.pc(A.ipf,:); Pf.k1 = c(1); Pf.k2 = c(1); Pf.g = c(2);
M = {'LBG',   @(Pr) simulate(Pr, THb, X0b, 'LBG', 20, Gt),                    P,  A.oL;
     'CE',    @(Pr) simulate(Pr, THb, X0b, 'RERG', 20, Gt),                   P,  A.oC;
     'ERG-N', @(Pr) simulate(Pr, THb, X0b, 'RERG', A.bestN(2), Gn),           Pn, A.oN;
     'PF+',   @(Pr) simulate(Pr, THb, X0b, 'PF', c(3:4), Gt),                 Pf, A.pf;
     'BLF',   @(Pr) simulate_blf(Pr, THb, X0b, A.bc(A.ibl,:)),                P,  A.bl};
fprintf('\n=== manipulator, p = 1: envelope %s rad, %d runs per method ===\n', mat2str(A.Rset), numel(A.Rset)*size(THb,2));
fprintf('%-6s | %6s %6s %6s %6s %6s | reference\n', 'method', 'med', 'p90', 'worst', 'unset', 'viol');
for m = 1:size(M, 1)
    ts = []; vx = [];
    for r = A.Rset, Pr = M{m,3}; Pr.r = r; q = M{m,2}(Pr); ts = [ts q.ts]; vx = [vx q.maxviolx]; end %#ok<AGROW>
    row = [median(ts) prctile(ts, 90) max(ts) sum(~isfinite(ts)) sum(vx > 1e-6)];
    o = M{m,4}; rf = [o.med o.p90 o.worst o.nst o.nvx];
    same = isequal(round(row, 2), round(rf, 2)); ok = ok && same;
    fprintf('%-6s | %6.2f %6.2f %6.2f %6d %6d | %6.2f %6.2f %6.2f %6d %6d  %s\n', M{m,1}, row, rf, ternary(same, 'MATCH', 'DIFFERS'));
end
fprintf('\ntransfer (state-violation %% | median ts | unsettled runs of %d):\n', size(THb,2));
for r = A.Rtr
    s = sprintf('r=%5.2f Gam/F=%5.2f', r, interp1(Gt.v, Gt.G, r)/F);
    for m = [1 4 5 2]
        Pr = M{m,3}; Pr.r = r; q = M{m,2}(Pr);
        s = [s sprintf(' | %s %5.1f%% %6.2f %3d', M{m,1}, 100*mean(q.maxviolx > 1e-6), median(q.ts), sum(~isfinite(q.ts)))]; %#ok<AGROW>
        if m == 1, ok = ok && all(q.maxviolx <= 1e-6); end
    end
    fprintf('%s\n', s);
end
% ---------------- p = 2 (paper, "Illustration of Theorem 2(c)"): LBG only
P2 = params_arm(); L = load('results/Gamma_cert_arm.mat'); G2 = L.Gt; A2 = load('results/arm.mat');
a = linspace(0, 2*pi, 9); a(end) = []; TH2 = [P2.th0 + P2.R*[cos(a); sin(a)], P2.th0];
[I, J] = ndgrid(1:size(TH2,2), 1:size(A2.X0set,2)); TH2b = TH2(:, I(:)); X02b = A2.X0set(:, J(:));
ts = []; vx = [];
for r = A2.Rset, Pr = P2; Pr.r = r; q = simulate(Pr, TH2b, X02b, 'LBG', 20, G2); ts = [ts q.ts]; vx = [vx q.maxviolx]; end %#ok<AGROW>
nt = 0; nvt = 0;
for r = A2.Rtr, Pr = P2; Pr.r = r; q = simulate(Pr, TH2b, X02b, 'LBG', 20, G2); nt = nt + numel(q.ts); nvt = nvt + sum(q.maxviolx > 1e-6); end
same = sum(~isfinite(ts)) == A2.oL.nst && sum(vx > 1e-6) == A2.oL.nvx; ok = ok && same && nvt == 0;
fprintf('\n=== manipulator, p = 2: envelope %d runs: unsettled %d (%.1f%%), violating %d (reference %d, %d) %s; transfer %d runs, violating %d\n', ...
    numel(ts), sum(~isfinite(ts)), 100*mean(~isfinite(ts)), sum(vx > 1e-6), A2.oL.nst, A2.oL.nvx, ternary(same, 'MATCH', 'DIFFERS'), nt, nvt);
fprintf('\nManipulator quick check: %s (compare the transfer rows with results/logs/arm_p1.log)\n', ternary(ok, 'ALL MATCH', 'SOME VALUES DIFFER'));

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
