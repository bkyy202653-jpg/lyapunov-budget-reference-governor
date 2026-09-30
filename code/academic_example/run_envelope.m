function run_envelope(Rset, tag, needSettle)
if nargin < 3, needSettle = false; end
if nargin < 1, Rset = [-0.9 -0.6 -0.3 0.3 0.5 0.7 0.85]; end
if nargin < 2, tag = ''; end
% Decisive test: PF+ must be (empirically) safe on the SAME envelope that LBG* certifies.
P = params_lbg(); P.c = [1.2; 0.8]; P.umax = 2.5;
T = load('results/tight_lbgtune.mat'); cfg = T.cfg; Gt = T.bestGt;
Pl = P; Pl.k1 = cfg(1); Pl.k2 = cfg(1); Pl.g = cfg(2);
[g1, g2] = ndgrid(-0.5:0.25:0.5, -0.6:0.3:0.6);
X0c = [g1(:) g2(:)]'; adm = false(1, size(X0c,2));
for j = 1:size(X0c,2)
    [~,~,z1,z2] = prestab(X0c(1,j), X0c(2,j), P.th0(1), P.th0(2), X0c(1,j), Pl);
    adm(j) = 0.5*(z1^2+z2^2) + 0.5*P.R^2/Pl.g <= interp1(Gt.v, Gt.G, X0c(1,j));
end
X0set = X0c(:, adm); fprintf('admissible x0 grid points: %d of %d\n', sum(adm), numel(adm));
a = linspace(0, 2*pi, 9); a(end) = [];
THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2));
THb = THt(:, I(:)); X0b = X0set(:, J(:));
    function [worst, safe, md] = evalset(Pp, mode, opt)
        worst = 0; safe = true; md = [];
        for r = Rset
            Pr = Pp; Pr.r = r; o = simulate(Pr, THb, X0b, mode, opt, Gt);
            safe = safe && all(o.maxviol <= 1e-6); worst = max(worst, max(o.ts)); md = [md o.ts]; %#ok<AGROW>
            if ~safe, return; end
        end
        md = median(md);
    end
[wL, sL, mL] = evalset(Pl, 'LBG', cfg(3));
fprintf('LBG*: safe=%d worst ts=%.2f median ts=%.2f\n', sL, wL, mL);
best = inf; bp = [];
for k = [1 2 3]
  for g = [0 5 20]
    Pp = P; Pp.k1 = k; Pp.k2 = k; Pp.g = g;
    for lam = [0.1 0.25 0.5 1 2 4]
      for rho = [0.02 0.05 0.1 0.2 0.3 0.5 1 inf]
        [w, sf, md] = evalset(Pp, 'PF', [lam rho]);
        if sf && (~needSettle || isfinite(w)) && (md < best || isempty(bp)), best = md; bp = [k g lam rho w md]; end
      end
    end
  end
end
if isempty(bp), fprintf('PF+: NO design is safe on the full envelope\n');
else, fprintf('PF+ envelope-safe: k=%g g=%g lam=%g rho=%g worst ts=%.2f median ts=%.2f\n', bp); end
save(['results/envelope' tag '.mat'], 'bp', 'wL', 'mL', 'sL', 'X0set', 'Rset')
end
