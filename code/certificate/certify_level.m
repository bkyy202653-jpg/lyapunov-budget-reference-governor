function [ok, info] = certify_level(vlo, vhi, gam, P, cap)
% Branch-and-bound interval certificate for one governor cell (n = 2 example, lbg/ plant).
% Proves: for all v in [vlo, vhi], theta in Theta, (z, thtilde) with |z|^2/2 + |thtilde|^2/(2g) <= gam,
%   |x1| <= (1-m)c1, |x2| <= (1-m)c2, |u| <= (1-m)umax.
% Exact reduction: the constraints depend on (z, h, v) only, h = thetahat = theta - thtilde, and the set of
% reachable h for given z is the ball |h - th0| <= R + sqrt(g(2 gam - |z|^2)).
% ok = 1 certified, 0 refuted (a feasible point violates), -1 undecided (box cap reached).
% Rounding: interval bounds are inflated by a relative 1e-12 (no directed rounding; the outward-rounded version is certify_level_rig.m).
if nargin < 5, cap = 2e6; end
lim = (1 - P.margin)*[P.c(1); P.c(2); P.umax];
rz = sqrt(2*gam); rh = P.R + sqrt(2*P.g*gam);
% box rows: [z1 z2 h1 h2 v], columns: boxes; L = lower, U = upper
L = [-rz; -rz; P.th0 - rh; vlo]; U = [rz; rz; P.th0 + rh; vhi];
fx1 = isfield(P, 'fixed1') && P.fixed1;          % theta1 known and not adapted: h1 = th0(1)
if fx1, L(3) = P.th0(1); U(3) = P.th0(1); end
W0 = max(U - L, 1e-12); nproc = 0; info.maxdepth = 0;
while ~isempty(L)
    nproc = nproc + size(L, 2);
    if nproc > cap, ok = -1; info.nproc = nproc; return; end
    % --- discard boxes that do not intersect the feasible set
    z2lo = sqlo(L(1,:), U(1,:)) + sqlo(L(2,:), U(2,:));
    hlo = sqrt((~fx1)*sqlo(L(3,:) - P.th0(1), U(3,:) - P.th0(1)) + sqlo(L(4,:) - P.th0(2), U(4,:) - P.th0(2)));
    keep = z2lo <= 2*gam & hlo <= P.R + sqrt(P.g*max(2*gam - z2lo, 0));
    L = L(:, keep); U = U(:, keep);
    if isempty(L), break; end
    % --- refutation at box centres that are feasible points
    C = (L + U)/2;
    inS = C(1,:).^2 + C(2,:).^2 <= 2*gam & ...
          vecnorm(C(3:4,:) - P.th0) <= P.R + sqrt(P.g*max(2*gam - C(1,:).^2 - C(2,:).^2, 0));
    vc = pointvals(C, P);
    if any(inS & any(abs(vc) > lim, 1)), ok = 0; info.nproc = nproc; return; end
    % --- interval enclosure
    [xl, xh] = boxvals(L, U, P);
    safe = all(max(abs(xl), abs(xh)) <= lim, 1);
    L = L(:, ~safe); U = U(:, ~safe);
    if isempty(L), break; end
    % --- split undecided boxes along the widest relative dimension
    [~, d] = max((U - L)./W0, [], 1);
    idx = sub2ind(size(L), d, 1:size(L, 2));
    mid = (L(idx) + U(idx))/2;
    L2 = L; U1 = U; U1(idx) = mid; L2(idx) = mid;
    L = [L, L2]; U = [U1, U];
end
ok = 1; info.nproc = nproc;
end

function y = sqlo(lo, hi)
y = min(lo.^2, hi.^2); y(lo <= 0 & hi >= 0) = 0;
end

function [lo, hi] = isq(a, b)
lo = min(a.^2, b.^2); hi = max(a.^2, b.^2); lo(a <= 0 & b >= 0) = 0;
end

function [lo, hi] = imul(a, b, c, d)
p = [a.*c; a.*d; b.*c; b.*d]; lo = min(p, [], 1); hi = max(p, [], 1);
end

function [lo, hi] = infl(lo, hi)
t = 1e-12*max(abs(lo), abs(hi)) + 1e-300; lo = lo - t; hi = hi + t;
end

function V = pointvals(C, P)
[xl, ~] = boxvals(C, C, P); V = xl;
end

function [lo, hi] = boxvals(L, U, P)
% interval enclosure of [x1; x2; u] (rows) over boxes; controller of lbg/prestab.m (phi1 = [x1^2;0], phi2 = [0;sin x1])
k1 = P.k1; k2 = P.k2; g = P.g;
if isfield(P, 'plant') && strcmp(P.plant, 'arm')
    % x1 = z1 + v, x2 = z2 - k1 z1, u = (-z1 - k2 z2 + (h1 - k1) x2 + 10 h2 sin x1)/b
    x1l = L(1,:) + L(5,:); x1h = U(1,:) + U(5,:);
    x2l = L(2,:) - k1*U(1,:); x2h = U(2,:) - k1*L(1,:);
    [p1l, p1h] = imul(L(3,:) - k1, U(3,:) - k1, x2l, x2h);
    snl = sin(max(x1l, -pi/2)); snh = sin(min(x1h, pi/2));
    bad = x1l < -pi/2 | x1h > pi/2; snl(bad) = -1; snh(bad) = 1;
    [p2l, p2h] = imul(L(4,:), U(4,:), snl, snh);
    ul = (-U(1,:) - k2*U(2,:) + p1l + 10*p2l)/P.b; uh = (-L(1,:) - k2*L(2,:) + p1h + 10*p2h)/P.b;
    [lo, hi] = infl([x1l; x2l; ul], [x1h; x2h; uh]);
    return
end
z1l = L(1,:); z1h = U(1,:); z2l = L(2,:); z2h = U(2,:);
h1l = L(3,:); h1h = U(3,:); h2l = L(4,:); h2h = U(4,:); vl = L(5,:); vh = U(5,:);
x1l = z1l + vl; x1h = z1h + vh;
[sl, sh] = isq(x1l, x1h);                          % x1^2
[q4l, q4h] = isq(sl, sh);                          % x1^4
% x2 = z2 - k1 z1 - h1 x1^2
[hsl, hsh] = imul(h1l, h1h, sl, sh);
x2l = z2l - k1*z1h - hsh; x2h = z2h - k1*z1l - hsl;
% A = -k1 - 2 h1 x1
[hxl, hxh] = imul(h1l, h1h, x1l, x1h);
Al = -k1 - 2*hxh; Ah = -k1 - 2*hxl;
% u = -z1 - k2 z2 + A (z2 - k1 z1) - h2 sin(x1) - g x1^4 (z1 - A z2)
wl = z2l - k1*z1h; wh = z2h - k1*z1l;              % z2 - k1 z1
[t1l, t1h] = imul(Al, Ah, wl, wh);
if all(x1l >= -pi/2 & x1h <= pi/2), snl = sin(x1l); snh = sin(x1h);
else, snl = sin(x1l); snh = sin(x1h); bad = ~(x1l >= -pi/2 & x1h <= pi/2); snl(bad) = -1; snh(bad) = 1; end
[t2l, t2h] = imul(h2l, h2h, snl, snh);
[azl, azh] = imul(Al, Ah, z2l, z2h);
el = z1l - azh; eh = z1h - azl;                    % z1 - A z2
[t3l, t3h] = imul(q4l, q4h, el, eh);
ul = -z1h - k2*z2h + t1l - t2h - g*t3h;
uh = -z1l - k2*z2l + t1h - t2l - g*t3l;
[lo, hi] = infl([x1l; x2l; ul], [x1h; x2h; uh]);
end
