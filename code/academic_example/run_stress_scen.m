function run_stress_scen(tag, umax, R)
% Recertification scenario (paper: tighter actuator S1, larger uncertainty S2): recertify Gamma_c, build envelope, run LBG, re-tune BLF, transfer test.
P = params_lbg(); P.c = [1.2; 0.8]; P.umax = umax; P.R = R; P.k1 = 1; P.k2 = 1; P.g = 20; P.margin = 0.02;
F = P.R^2/(2*P.g);
start_pool();   % 4 workers were used for the reported timings
fn = ['results/Gamma_cert_' tag '.mat'];
tic;
if exist(fn, 'file'), L = load(fn); Gt = L.Gt;
else, Gt = build_cert_generic(P, round((-0.90:0.02:0.86)*100)/100, 3e5, 14, fn); end
fprintf('[%s] umax=%.2f R=%.2f F=%.5f: table done (%.0fs), undecided %d, Gamma_c/F at {-0.9 -0.6 -0.3 0 0.3 0.5 0.7 0.85} = %s\n', ...
    tag, umax, R, F, toc, Gt.und, mat2str(interp1(Gt.v, Gt.G, [-0.9 -0.6 -0.3 0 0.3 0.5 0.7 0.85])/F, 3));
% envelope
[g1, g2] = ndgrid(-0.5:0.25:0.5, -0.6:0.3:0.6); X0c = [g1(:) g2(:)]'; adm = false(1, 25);
for j = 1:25
    [~,~,z1,z2] = prestab(X0c(1,j), X0c(2,j), P.th0(1), P.th0(2), X0c(1,j), P);
    adm(j) = 0.5*(z1^2 + z2^2) + F <= interp1(Gt.v, Gt.G, X0c(1,j));
end
X0set = X0c(:, adm);
Rall = [-0.6 -0.3 0.3 0.5 0.7]; Rset = Rall(interp1(Gt.v, Gt.G, Rall) >= 1.5*F);
Rtr = setdiff([-0.9 -0.8 -0.45 0.4 0.6 0.75 0.8 0.85 Rall], Rset);
a = linspace(0, 2*pi, 9); a(end) = []; THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
fprintf('[%s] admissible x0: %d/25, envelope setpoints %s, transfer %s, runs per setpoint %d\n', tag, sum(adm), mat2str(Rset), mat2str(Rtr), size(THb,2));
if isempty(X0set) || isempty(Rset), fprintf('[%s] empty envelope, stop\n', tag); return; end
% LBG
oL = evalset(@(Pr) simulate(Pr, THb, X0b, 'LBG', 20, Gt), P, Rset);
fprintf('[%s] LBG: med %.2f worst %.2f unsettled %d viol(state) %d viol(any) %d\n', tag, oL.med, oL.worst, oL.nst, oL.nvx, oL.nv);
% BLF grid
cfgs = [];
for k = [1 2 3], for g = [0 5 20], for kb1 = [0.3 0.45], for kb2 = [0.1 0.2 0.4], for wf = [0.5 1 2]
    cfgs = [cfgs; k g kb1 kb2 wf]; %#ok<AGROW>
end, end, end, end, end
nc = size(cfgs, 1); med = inf(nc,1); worst = inf(nc,1); nst = zeros(nc,1); nvx = zeros(nc,1); sat = zeros(nc,1);
parfor c = 1:nc
    o = evalset(@(Pr) simulate_blf(Pr, THb, X0b, cfgs(c,:)), P, Rset);
    med(c) = o.med; worst(c) = o.worst; nst(c) = o.nst; nvx(c) = o.nvx; sat(c) = o.sat;
end
okc = find(nvx == 0 & nst == 0);
fprintf('[%s] BLF: %d/%d designs admissible (zero state violations, all settled)\n', tag, numel(okc), nc);
[~, ifast] = min(med + 1e6*(nst > 0));
fprintf('[%s] fastest BLF regardless of violations: %s med %.2f worst %.2f viol runs %d sat %.1f%%\n', tag, mat2str(cfgs(ifast,:)), med(ifast), worst(ifast), nvx(ifast), sat(ifast));
if isempty(okc), bB = []; else, [~, i] = min(med(okc)); bB = cfgs(okc(i),:);
    fprintf('[%s] BLF admissible best: %s med %.2f worst %.2f sat %.2f%%\n', tag, mat2str(bB), med(okc(i)), worst(okc(i)), sat(okc(i))); end
% transfer
fprintf('[%s] transfer (state-violation %% | median ts):\n', tag);
for r = Rtr
    Pr = P; Pr.r = r; qL = simulate(Pr, THb, X0b, 'LBG', 20, Gt);
    s = sprintf('   r=%5.2f  Gam/F=%5.2f | LBG %5.1f%% %6.2f', r, interp1(Gt.v, Gt.G, r)/F, 100*mean(qL.maxviolx > 1e-6), median(qL.ts));
    if ~isempty(bB), qB = simulate_blf(Pr, THb, X0b, bB); s = [s sprintf(' | BLF %5.1f%% %6.2f', 100*mean(qB.maxviolx > 1e-6), median(qB.ts))]; end
    fprintf('%s\n', s);
end
save(['results/stress_' tag '.mat'], 'oL', 'cfgs', 'med', 'worst', 'nst', 'nvx', 'sat', 'bB', 'Rset', 'Rtr', 'X0set');
end

function o = evalset(fun, P, Rset)
ts = []; vx = []; v = []; sat = [];
for r = Rset
    Pr = P; Pr.r = r; q = fun(Pr);
    ts = [ts q.ts]; vx = [vx q.maxviolx]; sat = [sat q.satfrac]; %#ok<AGROW>
    if isfield(q, 'maxviol'), v = [v q.maxviol]; else, v = [v q.maxviolx]; end %#ok<AGROW>
end
o.med = median(ts); o.worst = max(ts); o.nst = sum(~isfinite(ts)); o.nvx = sum(vx > 1e-6); o.nv = sum(v > 1e-6); o.sat = 100*mean(sat);
end
