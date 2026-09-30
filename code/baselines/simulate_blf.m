function out = simulate_blf(P, TH, x0, opt)
% Baseline: BLF adaptive backstepping (Tee-Ge type, log-BLF on z1 and z2) tracking a 2nd-order prefiltered
% reference rf (rf(0) = x1(0), so z1(0) = 0). Plant x1' = x2 + th1 x1^2, x2' = sat(u) + th2 sin x1.
% opt = [k, g, kb1, kb2, wf]: k1 = k2 = k, adaptation gain g, BLF bounds kb1, kb2, prefilter bandwidth wf.
% |x1| <= |rf| + kb1 <= c1 needs kb1 <= c1 - |r|; kb2 is enlarged per run if |z2(0)| >= kb2 (BLF must start inside).
% V = log-BLFs + |thtilde|^2/(2g) gives V' = -k mu1 z1^2 - k mu2 z2^2 without saturation; input limits are NOT
% part of the BLF design, the actuator saturates.
k = opt(1); g = opt(2); kb1 = opt(3); kb2 = opt(4); wf = opt(5);
N = size(TH, 2);
if size(x0, 2) == 1, x0 = repmat(x0, 1, N); end
h0 = repmat(P.th0, 1, N);
isarm = isfield(P, 'plant') && strcmp(P.plant, 'arm');
if isarm, z20 = x0(2,:); else, z20 = x0(2,:) + h0(1,:).*x0(1,:).^2; end
KB2 = max(kb2, 1.1*abs(z20) + 0.01);
X = [x0; h0; x0(1,:); zeros(1, N)];            % x1 x2 h1 h2 rf rfd
nt = round(P.T/P.dt); dec = 10; L = floor(nt/dec) + 1;
rec.x1 = zeros(L, N); rec.u = rec.x1; rec.t = (0:L-1)'*P.dt*dec;
maxviolx = zeros(1, N); satfrac = zeros(1, N); broke = false(1, N); kk = 1;
for i = 0:nt
    [dX, u] = rhs(X);
    if mod(i, dec) == 0, rec.x1(kk,:) = X(1,:); rec.u(kk,:) = u; kk = kk + 1; end
    maxviolx = max(maxviolx, max([abs(X(1,:)) - P.c(1); abs(X(2,:)) - P.c(2)], [], 1));
    satfrac = satfrac + (abs(u) > P.umax)/nt;
    if i == nt, break; end
    k1_ = dX; k2_ = rhs(X + 0.5*P.dt*k1_); k3_ = rhs(X + 0.5*P.dt*k2_); k4_ = rhs(X + P.dt*k3_);
    X = X + P.dt/6*(k1_ + 2*k2_ + 2*k3_ + k4_);
    X(3:4,:) = P.th0 + (X(3:4,:) - P.th0).*min(1, P.R./max(vecnorm(X(3:4,:) - P.th0), eps));  % projection on Theta
    bad = ~all(isfinite(X), 1);
    broke = broke | bad; X(:, bad) = repmat([P.c(1)*2; 0; P.th0; 0; 0], 1, sum(bad));
end
band = abs(rec.x1 - P.r) > P.tol; ts = inf(1, N);
for j = 1:N
    last = find(band(:, j), 1, 'last');
    if isempty(last), ts(j) = 0; elseif last < L, ts(j) = rec.t(last + 1); end
end
ts(broke) = inf; maxviolx(broke) = inf;
out.ts = ts; out.maxviolx = maxviolx; out.satfrac = satfrac; out.rec = rec; out.KB2 = KB2;
out.iae = sum(abs(rec.x1 - P.r), 1)*P.dt*dec;

    function [dX, u] = rhs(X)
        x1 = X(1,:); x2 = X(2,:); h1 = X(3,:); h2 = X(4,:); rf = X(5,:); rfd = X(6,:);
        rfdd = wf^2*(P.r - rf) - 2*wf*rfd;
        z1 = x1 - rf;
        if isarm
            % x1' = x2, x2' = b u - th1 x2 - 10 th2 sin x1; alpha1 = -k z1 + rf'
            z2 = x2 + k*z1 - rfd;
            m1 = 1./max(kb1^2 - z1.^2, 1e-9); m2 = 1./max(KB2.^2 - z2.^2, 1e-9);
            a1d = -k*(x2 - rfd) + rfdd;
            u = (-k*z2 - (m1./m2).*z1 + h1.*x2 + 10*h2.*sin(x1) + a1d)/P.b;
            hd1 = -g*m2.*z2.*x2; hd2 = -10*g*m2.*z2.*sin(x1);
            if isfield(P, 'fixed1') && P.fixed1, hd1 = 0*hd1; end
            us = min(max(u, -P.umax), P.umax);
            dX = [x2; P.b*us - TH(1,:).*x2 - 10*TH(2,:).*sin(x1); hd1; hd2; rfd; rfdd];
            return
        end
        a1 = -k*z1 - h1.*x1.^2 + rfd;
        z2 = x2 - a1;
        m1 = 1./max(kb1^2 - z1.^2, 1e-9); m2 = 1./max(KB2.^2 - z2.^2, 1e-9);
        A = -k - 2*h1.*x1;
        hd1 = g*x1.^2.*(m1.*z1 - m2.*A.*z2);
        hd2 = g*m2.*z2.*sin(x1);
        u = -k*z2 - (m1./m2).*z1 - h2.*sin(x1) + A.*(x2 + h1.*x1.^2) - x1.^2.*hd1 + k*rfd + rfdd;
        us = min(max(u, -P.umax), P.umax);
        dX = [x2 + TH(1,:).*x1.^2; us + TH(2,:).*sin(x1); hd1; hd2; rfd; rfdd];
    end
end
