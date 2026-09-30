% Extra baselines on the CERTIFIED envelope (Gamma_cert.mat, 12 x0 x 9 theta x 5 setpoints = 540 runs):
%   BLF  : BLF adaptive backstepping + 2nd-order prefilter + actuator saturation (simulate_blf.m)
%   ERGn : standard ERG (Garone & Nicotra 2016): non-adaptive pre-stabilizer, thetahat = th0, threshold for the
%          nominal model (R = 0), margin Gamma_nom(v) - V_z
%   CE   : certainty-equivalence adaptive ERG (ablation): adaptive pre-stabilizer, LBG's certified table, margin
%          Gamma_c(v) - V_z, i.e. the parameter-error energy is ignored
%   LBG  : k=1, g=20, kappa=20, certified table;  PF+ : prefilter design tuned earlier (envelope_cert_settle.mat)
% Selection rule for every tuned method: zero STATE violations and all runs settled on the envelope, then the
% smallest median settling time. Then all deployed designs are run on 8 transfer setpoints.
P = params_lbg(); P.c = [1.2; 0.8]; P.umax = 2.5;
C = load('results/Gamma_cert.mat'); Gc = C.Gt;
E = load('results/envelope_cert_settle.mat'); bp = E.bp; Rset = E.Rset; X0set = E.X0set;
a = linspace(0, 2*pi, 9); a(end) = []; THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
Rtr = [-0.9 -0.8 -0.45 0.4 0.6 0.75 0.8 0.85];
S = struct();
% ---------------- BLF ----------------
best = inf; tic;
for k = [1 2 3 5], for g = [0 5 20], for kb1 = [0.3 0.45], for kb2 = [0.2 0.4 0.6], for wf = [0.5 1 2 4]
    o = evalm(@(Pr) simulate_blf(Pr, THb, X0b, [k g kb1 kb2 wf]), P, Rset);
    if o.vx <= 1e-6 && o.nst == 0 && o.med < best, best = o.med; S.BLF = o; S.BLF.opt = [k g kb1 kb2 wf]; end
end, end, end, end, end
fprintf('BLF tuning done (%.0fs)\n', toc);
% ---------------- ERGn ----------------
best = inf;
for k = [1 2 3]
    Pn = P; Pn.k1 = k; Pn.k2 = k; Pn.g = 0; Pq = Pn; Pq.R = 0;
    Gn = make_table(Pq, -1.0:0.02:1.2, 3000);
    for kap = [5 20 100 500]
        o = evalm(@(Pr) simulate(Pr, THb, X0b, 'RERG', kap, Gn), Pn, Rset);
        fprintf('  ERGn k=%g kap=%g: med %.2f worst %.2f unsettled %d state-viol runs %d (max %.3f)\n', k, kap, o.med, o.worst, o.nst, o.nvx, o.vx);
        if o.vx <= 1e-6 && o.nst == 0 && o.med < best, best = o.med; S.ERGn = o; S.ERGn.opt = [k kap]; S.ERGn.Gt = Gn; end
        if ~isfield(S, 'ERGn_any') || o.med < S.ERGn_any.med, S.ERGn_any = o; S.ERGn_any.opt = [k kap]; S.ERGn_any.Gt = Gn; end
    end
end
% ---------------- CE ablation ----------------
Pc = P; Pc.k1 = 1; Pc.k2 = 1; Pc.g = 20;
for kap = [20 100]
    o = evalm(@(Pr) simulate(Pr, THb, X0b, 'RERG', kap, Gc), Pc, Rset);
    fprintf('  CE kap=%g: med %.2f worst %.2f unsettled %d state-viol runs %d (max %.3f) cmd>umax runs %d\n', kap, o.med, o.worst, o.nst, o.nvx, o.vx, o.nvu);
    if kap == 20, S.CE = o; S.CE.opt = kap; end
end
% ---------------- LBG and PF+ ----------------
Pl = P; Pl.k1 = 1; Pl.k2 = 1; Pl.g = 20;
S.LBG = evalm(@(Pr) simulate(Pr, THb, X0b, 'LBG', 20, Gc), Pl, Rset); S.LBG.opt = [1 20 20];
Pp = P; Pp.k1 = bp(1); Pp.k2 = bp(1); Pp.g = bp(2);
S.PF = evalm(@(Pr) simulate(Pr, THb, X0b, 'PF', bp(3:4), Gc), Pp, Rset); S.PF.opt = bp(1:4);
% ---------------- table ----------------
fprintf('\n=== certified envelope, %d runs per method ===\n', numel(Rset)*size(THb,2));
fprintf('%-6s %-22s %6s %6s %6s %5s %8s %8s %8s %7s\n', 'method', 'design', 'med', 'p90', 'worst', 'unst', 'Xviol%', 'maxXv', 'sat%', 'IAEmed');
names = {'LBG', 'PF', 'BLF', 'ERGn', 'CE'};
for m = names
    if ~isfield(S, m{1}), fprintf('%-6s no admissible design (best unconstrained shown below)\n', m{1}); continue; end
    o = S.(m{1}); prow(m{1}, o);
end
if ~isfield(S, 'ERGn'), prow('ERGn*', S.ERGn_any); end
% ---------------- transfer ----------------
fprintf('\n=== transfer setpoints (deployed designs, state-violation %% of 108 runs | median ts) ===\n');
fprintf('%6s', 'r'); for m = names, fprintf(' | %-14s', m{1}); end; fprintf('\n');
dep = S; if ~isfield(dep, 'ERGn'), dep.ERGn = S.ERGn_any; end
for r = Rtr
    fprintf('%6.2f', r);
    for m = names
        o = dep.(m{1}); Pr = P;
        switch m{1}
            case 'LBG', Pr = Pl; Pr.r = r; q = simulate(Pr, THb, X0b, 'LBG', 20, Gc);
            case 'PF',  Pr = Pp; Pr.r = r; q = simulate(Pr, THb, X0b, 'PF', bp(3:4), Gc);
            case 'BLF', Pr.r = r; q = simulate_blf(Pr, THb, X0b, o.opt);
            case 'ERGn', Pr.k1 = o.opt(1); Pr.k2 = o.opt(1); Pr.g = 0; Pr.r = r; q = simulate(Pr, THb, X0b, 'RERG', o.opt(2), o.Gt);
            case 'CE',  Pr = Pc; Pr.r = r; q = simulate(Pr, THb, X0b, 'RERG', 20, Gc);
        end
        fprintf(' | %5.1f%% %7.2f', 100*mean(q.maxviolx > 1e-6), median(q.ts));
    end
    fprintf('\n');
end
S = rmfield(S, intersect(fieldnames(S), {'ERGn_any'}));
save results/baselines2.mat S

function o = evalm(fun, Pp, Rset)
ts = []; vx = []; sat = []; iae = []; vu = [];
for r = Rset
    Pr = Pp; Pr.r = r; q = fun(Pr);
    ts = [ts q.ts]; vx = [vx q.maxviolx]; sat = [sat q.satfrac]; iae = [iae q.iae]; %#ok<AGROW>
    if isfield(q, 'maxviol'), vu = [vu q.maxviol > 1e-6 & q.maxviolx <= 1e-6]; else, vu = [vu false(size(q.ts))]; end %#ok<AGROW>
end
f = isfinite(ts);
o.med = median(ts); o.p90 = prctile(ts, 90); o.worst = max(ts); o.nst = sum(~f);
o.vx = max(vx); o.nvx = sum(vx > 1e-6); o.pvx = 100*mean(vx > 1e-6); o.sat = 100*mean(sat); o.iae = median(iae);
o.nvu = sum(vu);
end

function prow(name, o)
fprintf('%-6s %-22s %6.2f %6.2f %6.2f %5d %7.1f%% %8.3f %7.2f%% %7.2f\n', name, mat2str(o.opt, 3), o.med, o.p90, o.worst, o.nst, o.pvx, o.vx, o.sat, o.iae);
end
