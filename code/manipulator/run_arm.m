% Physical example: single-link manipulator (params_arm.m). Certified table, envelope, LBG vs PF+, BLF, ERGn, CE.
P = params_arm(); F = P.R^2/(2*P.g);
start_pool();   % 4 workers were used for the reported timings
vn = round((-1.20:0.02:1.20)*100)/100;
if exist('results/Gamma_cert_arm.mat', 'file'), L = load('results/Gamma_cert_arm.mat'); Gt = L.Gt;
else, tic; Gt = build_cert_generic(P, vn, 3e5, 14, 'results/Gamma_cert_arm.mat'); fprintf('table %.0fs\n', toc); end
Pn = P; Pn.g = 0; Pn.R = 0;                                   % nominal (non-adaptive) threshold for ERGn
if exist('results/Gamma_cert_arm_nom.mat', 'file'), L = load('results/Gamma_cert_arm_nom.mat'); Gn = L.Gt;
else, Gn = build_cert_generic(Pn, vn, 3e5, 14, 'results/Gamma_cert_arm_nom.mat'); end
vq = [-1.2 -1.0 -0.6 -0.3 0 0.3 0.6 1.0 1.2];
fprintf('arm: F=%.5f, undecided %d, Gamma_c/F at %s = %s, L_Gamma %.3f\n', F, Gt.und, mat2str(vq), mat2str(interp1(Gt.v, Gt.G, vq)/F, 3), Gt.LG);
% envelope
[g1, g2] = ndgrid(-0.5:0.25:0.5, -0.6:0.3:0.6); X0c = [g1(:) g2(:)]'; adm = false(1, 25);
for j = 1:25
    [~,~,z1,z2] = prestab(X0c(1,j), X0c(2,j), P.th0(1), P.th0(2), X0c(1,j), P);
    adm(j) = 0.5*(z1^2 + z2^2) + F <= interp1(Gt.v, Gt.G, X0c(1,j));
end
X0set = X0c(:, adm);
Rall = [-1.0 -0.6 -0.3 0.3 0.6 1.0]; Rset = Rall(interp1(Gt.v, Gt.G, Rall) >= 1.5*F);
Rtr = [setdiff(Rall, Rset), -1.2 -0.9 -0.45 0.45 0.8 1.1 1.2];
a = linspace(0, 2*pi, 9); a(end) = []; THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
fprintf('admissible x0 %d/25, envelope %s, transfer %s, %d runs per setpoint\n', sum(adm), mat2str(Rset), mat2str(Rtr), size(THb,2));
% ---- LBG
oL = evalset(@(Pr) simulate(Pr, THb, X0b, 'LBG', 20, Gt), P, Rset);
% ---- CE ablation and ERGn
oC = evalset(@(Pr) simulate(Pr, THb, X0b, 'RERG', 20, Gt), P, Rset);
bestN = []; oN = [];
for k = [1 2 3], for kap = [20 100 500]
    Q = Pn; Q.k1 = k; Q.k2 = k;
    if k > 1, continue; end                                    % nominal table certified for k = 1 only
    o = evalset(@(Pr) simulate(Pr, THb, X0b, 'RERG', kap, Gn), Q, Rset);
    if isempty(oN) || (o.nvx == 0 && o.med < oN.med), oN = o; bestN = [k kap]; end
end, end
% ---- PF+ grid
pc = []; for k = [1 2 3], for g = [0 5 20], for lam = [0.25 0.5 1 2 4], for rho = [0.1 0.3 1 inf], pc = [pc; k g lam rho]; end, end, end, end %#ok<AGROW>
[pf, ipf] = tune(pc, @(c, Pr) simulate(setk(Pr, c), THb, X0b, 'PF', c(3:4), Gt), P, Rset);
% ---- BLF grid
bc = []; for k = [1 2 3], for g = [0 5 20], for kb1 = [0.2 0.3], for kb2 = [0.2 0.4 0.6], for wf = [0.5 1 2], bc = [bc; k g kb1 kb2 wf]; end, end, end, end, end %#ok<AGROW>
[bl, ibl] = tune(bc, @(c, Pr) simulate_blf(Pr, THb, X0b, c), P, Rset);
% ---- table
fprintf('\n=== arm, certified envelope ===\n%-5s %-20s %6s %6s %6s %5s %6s\n', 'meth', 'design', 'med', 'p90', 'worst', 'unst', 'Xviol');
pr('LBG', [1 20 20], oL); pr('CE', 20, oC); pr('ERGn', bestN, oN);
if ~isempty(ipf), pr('PF+', pc(ipf,:), pf); else, fprintf('PF+  no admissible design\n'); end
if ~isempty(ibl), pr('BLF', bc(ibl,:), bl); else, fprintf('BLF  no admissible design\n'); end
fprintf('\n=== transfer (state-violation %% | median ts) ===\n');
for r = Rtr
    Pr = P; Pr.r = r; q = simulate(Pr, THb, X0b, 'LBG', 20, Gt);
    s = sprintf('r=%5.2f Gam/F=%6.2f | LBG %5.1f%% %6.2f', r, interp1(Gt.v, Gt.G, r)/F, 100*mean(q.maxviolx > 1e-6), median(q.ts));
    if ~isempty(ipf), q = simulate(setk(Pr, pc(ipf,:)), THb, X0b, 'PF', pc(ipf,3:4), Gt); s = [s sprintf(' | PF+ %5.1f%% %6.2f', 100*mean(q.maxviolx > 1e-6), median(q.ts))]; end
    if ~isempty(ibl), q = simulate_blf(Pr, THb, X0b, bc(ibl,:)); s = [s sprintf(' | BLF %5.1f%% %6.2f', 100*mean(q.maxviolx > 1e-6), median(q.ts))]; end
    q = simulate(Pr, THb, X0b, 'RERG', 20, Gt); s = [s sprintf(' | CE %5.1f%% %6.2f', 100*mean(q.maxviolx > 1e-6), median(q.ts))];
    fprintf('%s\n', s);
end
save results/arm.mat oL oC oN bestN pf ipf pc bl ibl bc Rset Rtr X0set

function Pr = setk(Pr, c), Pr.k1 = c(1); Pr.k2 = c(1); Pr.g = c(2); end
function pr(n, d, o), fprintf('%-5s %-20s %6.2f %6.2f %6.2f %5d %6d\n', n, mat2str(d, 3), o.med, o.p90, o.worst, o.nst, o.nvx); end
function [best, ib] = tune(cfgs, fun, P, Rset)
nc = size(cfgs, 1); M = cell(nc, 1);
parfor c = 1:nc, M{c} = evalset(@(Pr) fun(cfgs(c,:), Pr), P, Rset); end
ib = []; best = [];
for c = 1:nc
    o = M{c};
    if o.nvx == 0 && o.nst == 0 && (isempty(best) || o.med < best.med), best = o; ib = c; end
end
end
function o = evalset(fun, P, Rset)
ts = []; vx = [];
for r = Rset, Pr = P; Pr.r = r; q = fun(Pr); ts = [ts q.ts]; vx = [vx q.maxviolx]; end %#ok<AGROW>
o.med = median(ts); o.p90 = prctile(ts, 90); o.worst = max(ts); o.nst = sum(~isfinite(ts)); o.nvx = sum(vx > 1e-6);
end
