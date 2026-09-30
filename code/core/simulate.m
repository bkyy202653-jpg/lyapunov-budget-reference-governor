function out = simulate(P, TH, x0, mode, opt, Gt)
% Vectorized RK4 simulation (fixed step P.dt) over a batch of true parameters TH (2 x N) and initial states x0.
% Plant: academic example (default) or single-link manipulator (P.plant = 'arm'); the actuator always saturates.
% mode: 'LBG'    Lyapunov-budget governor, budget s' = -k1 z1^2 - k2 z2^2 + (dV/dv) v'  (paper, eq. of Lemma 2)
%       'NODISS' budget without the dissipation term, s' = (dV/dv) v' (not used in the paper)
%       'RERG'   margin Gamma(v) - V_z: ERG-N baseline (with P.g = 0 and a nominal table) or CE ablation
%                (adaptive loop, certified table)
%       'PF'     PF+ baseline, rate-limited first-order prefilter v' = sat(lam (r-v), rho)
% opt: kappa (LBG/NODISS/RERG) or [lam, rho] (PF).  Gt: threshold table struct(v, G, Dmax).
N = size(TH, 2);
if size(x0, 2) == 1, x0 = repmat(x0, 1, N); end
h0 = repmat(P.th0, 1, N);
v0 = x0(1, :);
[~, ~, z1, z2] = prestab(x0(1,:), x0(2,:), h0(1,:), h0(2,:), v0, P);
s0 = 0.5*(z1.^2 + z2.^2) + 0.5*P.R^2/max(P.g, eps);
X = [x0; h0; v0; s0];                 % state: x1 x2 h1 h2 v s
nt = round(P.T/P.dt);
dec = 10;
L = floor(nt/dec) + 1;
rec.x1 = zeros(L, N); rec.x2 = rec.x1; rec.u = rec.x1; rec.v = rec.x1;
rec.s = rec.x1; rec.V = rec.x1; rec.Gv = rec.x1; rec.t = (0:L-1)'*P.dt*dec;
maxviol = zeros(1, N); maxviolx = zeros(1, N); satfrac = zeros(1, N); minDelta = inf(1, N);
k = 1;
for i = 0:nt
    [dX, u, V, Gv] = rhs(X);
    if mod(i, dec) == 0
        rec.x1(k,:) = X(1,:); rec.x2(k,:) = X(2,:); rec.u(k,:) = u;
        rec.v(k,:) = X(5,:); rec.s(k,:) = X(6,:); rec.V(k,:) = V; rec.Gv(k,:) = Gv;
        k = k + 1;
    end
    viol = max([abs(X(1,:)) - P.c(1); abs(X(2,:)) - P.c(2); abs(u) - P.umax], [], 1);
    maxviol = max(maxviol, viol);
    maxviolx = max(maxviolx, max([abs(X(1,:)) - P.c(1); abs(X(2,:)) - P.c(2)], [], 1));
    satfrac = satfrac + (abs(u) > P.umax)/nt;
    if any(strcmp(mode, {'LBG', 'NODISS'})), minDelta = min(minDelta, Gv - X(6,:)); end
    if i == nt, break; end
    k1_ = dX;
    k2_ = rhs(X + 0.5*P.dt*k1_);
    k3_ = rhs(X + 0.5*P.dt*k2_);
    k4_ = rhs(X + P.dt*k3_);
    X = X + P.dt/6*(k1_ + 2*k2_ + 2*k3_ + k4_);
end
% settling time of x1 to the band
band = abs(rec.x1 - P.r) > P.tol;          % absolute settling band
ts = inf(1, N);
for j = 1:N
    last = find(band(:, j), 1, 'last');
    if isempty(last), ts(j) = 0; elseif last < L, ts(j) = rec.t(last + 1); end
end
out.ts = ts; out.maxviol = maxviol; out.maxviolx = maxviolx; out.satfrac = satfrac; out.minDelta = minDelta;
out.rec = rec; out.budgetOK = all(rec.s >= rec.V - 1e-9, 'all');
out.iae = sum(abs(rec.x1 - P.r), 1)*P.dt*dec;

    function [dX, u, V, Gv] = rhs(X)
        x1 = X(1,:); x2 = X(2,:); h1 = X(3,:); h2 = X(4,:); v = X(5,:); s = X(6,:);
        [u, hd, z1, z2, dVdv] = prestab(x1, x2, h1, h2, v, P);
        us = min(max(u, -P.umax), P.umax);           % plant always saturates
        switch mode
            case {'LBG', 'NODISS'}
                Gv = interp1(Gt.v, Gt.G, v, 'linear', 0);
                D = min(max(Gv - s, 0), Gt.Dmax);
                e = P.r - v;
                vd = opt*D.*e./max(abs(e), P.eta);
                if strcmp(mode, 'LBG')
                    sd = -P.k1*z1.^2 - P.k2*z2.^2 + dVdv.*vd;   % exact: s - V is constant
                else
                    sd = dVdv.*vd;                              % ablation: no dissipation credit
                end
            case 'RERG'   % robust ERG on non-adaptive pre-stabilizer (P.g = 0): DSM = Gamma_rob(v) - Vz
                Gv = interp1(Gt.v, Gt.G, v, 'linear', 0);
                Vz = 0.5*(z1.^2 + z2.^2);
                D = min(max(Gv - Vz, 0), Gt.Dmax);
                e = P.r - v;
                vd = opt*D.*e./max(abs(e), P.eta);
                sd = zeros(size(v));
            case 'PF'
                Gv = zeros(size(v));
                vd = min(max(opt(1)*(P.r - v), -opt(2)), opt(2));
                sd = zeros(size(v));
        end
        tt = TH - [h1; h2];
        if P.g > 0, V = 0.5*(z1.^2 + z2.^2) + 0.5*sum(tt.^2, 1)/P.g; else, V = 0.5*(z1.^2 + z2.^2); end
        if isfield(P, 'plant') && strcmp(P.plant, 'arm')
            dX = [x2; P.b*us - TH(1,:).*x2 - 10*TH(2,:).*sin(x1); hd; vd; sd];
        else
            dX = [x2 + TH(1,:).*x1.^2; us + TH(2,:).*sin(x1); hd; vd; sd];
        end
    end
end
