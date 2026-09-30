% Quick reproduction of the academic-example numbers: re-evaluates the DEPLOYED designs (the tuning searches of
% run_envelope.m, run_baselines2.m and run_stress_scen.m are skipped) and compares with the reference results.
% Table 1 (envelope + transfer), Section "Tighter actuator and larger uncertainty" (S1: umax = 1.8, S2: R = 0.4).
% Needs in results/: Gamma_cert*.mat, envelope_cert_settle.mat, baselines2.mat, stress_S1.mat, stress_S2.mat.
P = params_lbg(); P.c = [1.2; 0.8]; P.umax = 2.5;
C = load('results/Gamma_cert.mat'); Gc = C.Gt;
E = load('results/envelope_cert_settle.mat'); bp = E.bp; Rset = E.Rset; X0set = E.X0set;
B = load('results/baselines2.mat'); S = B.S;
% the envelope initial states are exactly the 5 x 5 grid points that satisfy Assumption 2 with Gamma_cert.mat
Pa2 = P; Pa2.k1 = 1; Pa2.k2 = 1; Pa2.g = 20; [g1, g2] = ndgrid(-0.5:0.25:0.5, -0.6:0.3:0.6); X0c = [g1(:) g2(:)]'; adm = false(1, 25);
for j = 1:25
    [~,~,z1,z2] = prestab(X0c(1,j), X0c(2,j), P.th0(1), P.th0(2), X0c(1,j), Pa2);
    adm(j) = 0.5*(z1^2 + z2^2) + P.R^2/(2*Pa2.g) <= interp1(Gc.v, Gc.G, X0c(1,j));
end
fprintf('envelope initial states = certified admissible grid points (%d of 25): %d\n', sum(adm), isequal(sortrows(X0c(:,adm)'), sortrows(X0set')));
a = linspace(0, 2*pi, 9); a(end) = []; THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
Rtr = [-0.9 -0.8 -0.45 0.4 0.6 0.75 0.8 0.85];
% deployed designs (ERG-N: best design of run_baselines2.m, k = 3, kappa = 500, see results/logs/baselines2.log)
Pl = P; Pl.k1 = 1; Pl.k2 = 1; Pl.g = 20;
Pp = P; Pp.k1 = bp(1); Pp.k2 = bp(1); Pp.g = bp(2);
Pn = P; Pn.k1 = 3; Pn.k2 = 3; Pn.g = 0; Pq = Pn; Pq.R = 0; Gn = make_table(Pq, -1.0:0.02:1.2, 3000);
M = {'LBG',  @(Pr) simulate(Pr, THb, X0b, 'LBG', 20, Gc),              Pl;
     'PF+',  @(Pr) simulate(Pr, THb, X0b, 'PF', bp(3:4), Gc),           Pp;
     'BLF',  @(Pr) simulate_blf(Pr, THb, X0b, S.BLF.opt),               P;
     'ERG-N', @(Pr) simulate(Pr, THb, X0b, 'RERG', 500, Gn),            Pn;
     'CE',   @(Pr) simulate(Pr, THb, X0b, 'RERG', 20, Gc),              Pl};
% reference: envelope from baselines2.mat (ERG-N from the log), transfer violation % from baselines2.log
ref = [S.LBG.med S.LBG.worst S.LBG.nst S.LBG.pvx; S.PF.med S.PF.worst S.PF.nst S.PF.pvx; ...
       S.BLF.med S.BLF.worst S.BLF.nst S.BLF.pvx; 4.16 Inf 204 0; S.CE.med S.CE.worst S.CE.nst S.CE.pvx];
reftr = [0 0 0 0 0 0 0 0; 20.4 1.9 0 0 0 1.9 10.2 18.5; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 11.1; 5.6 5.6 0 0 25.0 32.4 32.4 32.4];
fprintf('\n=== Table 1: academic example, envelope (%d runs per method) and %d transfer setpoints ===\n', numel(Rset)*size(THb,2), numel(Rtr));
fprintf('%-6s | %7s %7s %6s %7s | %-22s | %s\n', 'method', 'med ts', 'worst', 'unset', 'viol%', 'reference (same order)', 'transfer viol% per setpoint (max)');
ok = true;
for m = 1:size(M, 1)
    ts = []; vx = [];
    for r = Rset, Pr = M{m,3}; Pr.r = r; q = M{m,2}(Pr); ts = [ts q.ts]; vx = [vx q.maxviolx]; end %#ok<AGROW>
    row = [median(ts) max(ts) sum(~isfinite(ts)) 100*mean(vx > 1e-6)];
    tr = zeros(size(Rtr));
    for i = 1:numel(Rtr), Pr = M{m,3}; Pr.r = Rtr(i); q = M{m,2}(Pr); tr(i) = 100*mean(q.maxviolx > 1e-6); end
    same = isequal(round(row, 2), round(ref(m,:), 2)) && isequal(round(tr, 1), reftr(m,:));
    ok = ok && same;
    fprintf('%-6s | %7.2f %7.2f %6d %6.1f%% | %5.2f %6.2f %4d %5.1f%% | %s (max %.1f%%)  %s\n', M{m,1}, row, ref(m,:), ...
        mat2str(round(tr, 1)), max(tr), ternary(same, 'MATCH', 'DIFFERS'));
end
% ---- S1 (umax = 1.8) and S2 (R = 0.4): LBG with the recomputed table, BLF re-tuned for the scenario (bB)
for sc = {'S1', 1.8, 0.25; 'S2', 2.5, 0.40}'
    Q = params_lbg(); Q.c = [1.2; 0.8]; Q.umax = sc{2}; Q.R = sc{3}; Q.k1 = 1; Q.k2 = 1; Q.g = 20; Q.margin = 0.02;
    T = load(['results/stress_' sc{1} '.mat']); G = load(['results/Gamma_cert_' sc{1} '.mat']); Gt = G.Gt;
    th = [Q.th0 + Q.R*[cos(a); sin(a)], Q.th0];
    [I, J] = ndgrid(1:size(th,2), 1:size(T.X0set,2)); THs = th(:, I(:)); X0s = T.X0set(:, J(:));
    ts = []; vx = []; vb = [];
    for r = T.Rset
        Pr = Q; Pr.r = r; q = simulate(Pr, THs, X0s, 'LBG', 20, Gt); ts = [ts q.ts]; vx = [vx q.maxviolx]; %#ok<AGROW>
        q = simulate_blf(Pr, THs, X0s, T.bB); vb = [vb q.maxviolx]; %#ok<AGROW>
    end
    trL = zeros(size(T.Rtr)); trB = trL;
    for i = 1:numel(T.Rtr)
        Pr = Q; Pr.r = T.Rtr(i);
        q = simulate(Pr, THs, X0s, 'LBG', 20, Gt); trL(i) = 100*mean(q.maxviolx > 1e-6);
        q = simulate_blf(Pr, THs, X0s, T.bB); trB(i) = 100*mean(q.maxviolx > 1e-6);
    end
    same = abs(median(ts) - T.oL.med) < 5e-3 && abs(max(ts) - T.oL.worst) < 5e-3 && sum(vx > 1e-6) == T.oL.nvx;
    ok = ok && same && all(trL == 0);
    fprintf('\n[%s] umax=%.1f R=%.2f, %d runs per setpoint. LBG envelope: med %.2f worst %.2f unsettled %d viol %d (reference %.2f %.2f %d %d) %s\n', ...
        sc{1}, Q.umax, Q.R, size(THs,2), median(ts), max(ts), sum(~isfinite(ts)), sum(vx > 1e-6), T.oL.med, T.oL.worst, T.oL.nst, T.oL.nvx, ternary(same, 'MATCH', 'DIFFERS'));
    fprintf('[%s] BLF %s envelope violating runs: %d\n', sc{1}, mat2str(T.bB), sum(vb > 1e-6));
    fprintf('[%s] transfer r = %s\n       LBG viol%% %s\n       BLF viol%% %s\n', sc{1}, mat2str(T.Rtr), mat2str(round(trL, 1)), mat2str(round(trB, 1)));
end
fprintf('\nAcademic example quick check: %s (compare the transfer rows with results/logs/stress.log)\n', ternary(ok, 'ALL MATCH', 'SOME VALUES DIFFER'));

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
