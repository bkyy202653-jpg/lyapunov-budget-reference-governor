function Gt = make_table(P, vg, Nq)
% Gamma(v) table with conservative 3-point lower envelope and 5% shrink.
if nargin < 3, Nq = 3000; end
G = compute_Gamma(P, vg, Nq, 1);
Glo = min([G; [G(2:end) G(end)]; [G(1) G(1:end-1)]], [], 1);
Gt.v = vg; Gt.G = 0.95*Glo; Gt.Dmax = 1;
end
