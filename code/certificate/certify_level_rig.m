function [ok, info] = certify_level_rig(vlo, vhi, gam, P, cap)
% Rigorous version of certify_level.m: every floating-point operation in the interval enclosures and in the
% box-discarding test is rounded outward by one ulp (dn/up below). IEEE-754 round-to-nearest has error
% <= 0.5 ulp for + - * /, so each enclosure contains the exact value. sin is widened by 2 ulp, assuming the
% library sine is accurate to 1 ulp (stated as an assumption in the paper). Refutation at box centres is
% unchanged: it can only make a test fail, never make it pass.
% ok = 1 certified, 0 refuted, -1 undecided (box cap reached).
if nargin < 5, cap = 3e6; end
% decimal data (R, limits, margin) are not exactly representable: round them conservatively
Rr = up(up(P.R));                                   % >= the exact decimal radius
lim = dn((1 - P.margin)*[P.c(1); P.c(2); P.umax]*(1 - 1e-14));
rz = up(sqrt(up(2*gam))); rh = up(Rr + up(sqrt(up(up(2*P.g)*gam))));
L = [-rz; -rz; dn(P.th0 - rh); vlo]; U = [rz; rz; up(P.th0 + rh); vhi];
fx1 = isfield(P, 'fixed1') && P.fixed1;
if fx1, L(3) = P.th0(1); U(3) = P.th0(1); end
W0 = max(U - L, 1e-12); nproc = 0;
g2 = up(2*gam);
while ~isempty(L)
    nproc = nproc + size(L, 2);
    if nproc > cap, ok = -1; info.nproc = nproc; return; end
    % --- discard only boxes that provably miss the feasible set (lower bounds rounded down, radius rounded up)
    z2lo = dn(sqlo(L(1,:), U(1,:)) + sqlo(L(2,:), U(2,:)));
    h2lo = dn((~fx1)*sqlo(dn(L(3,:) - P.th0(1)), up(U(3,:) - P.th0(1))) + sqlo(dn(L(4,:) - P.th0(2)), up(U(4,:) - P.th0(2))));
    rad = up(Rr + up(sqrt(up(P.g*max(up(g2 - z2lo), 0)))));
    keep = z2lo <= g2 & dn(sqrt(max(h2lo, 0))) <= rad;
    L = L(:, keep); U = U(:, keep);
    if isempty(L), break; end
    % --- refutation at box centres (plain arithmetic is sufficient to reject)
    C = (L + U)/2;
    inS = C(1,:).^2 + C(2,:).^2 <= 2*gam & ...
          vecnorm(C(3:4,:) - P.th0) <= P.R + sqrt(P.g*max(2*gam - C(1,:).^2 - C(2,:).^2, 0));
    [vc, ~] = boxvals(C, C, P);
    if any(inS & any(abs(vc) > (1 - P.margin)*[P.c(1); P.c(2); P.umax], 1)), ok = 0; info.nproc = nproc; return; end
    % --- outward-rounded interval enclosure
    [xl, xh] = boxvals(L, U, P);
    safe = all(max(abs(xl), abs(xh)) <= lim, 1);
    L = L(:, ~safe); U = U(:, ~safe);
    if isempty(L), break; end
    [~, d] = max((U - L)./W0, [], 1);
    idx = sub2ind(size(L), d, 1:size(L, 2));
    mid = (L(idx) + U(idx))/2;
    L2 = L; U1 = U; U1(idx) = mid; L2(idx) = mid;
    L = [L, L2]; U = [U1, U];
end
ok = 1; info.nproc = nproc;
end

function y = dn(x), y = x - eps(x); end
function y = up(x), y = x + eps(x); end

function y = sqlo(lo, hi)                     % rounded-down lower bound of x^2 on [lo, hi]
y = dn(min(lo.^2, hi.^2)); y(lo <= 0 & hi >= 0) = 0;
end
function [l, h] = isq(a, b)
l = dn(min(a.^2, b.^2)); h = up(max(a.^2, b.^2)); l(a <= 0 & b >= 0) = 0;
end
function [l, h] = imul(a, b, c, d)
p = [a.*c; a.*d; b.*c; b.*d]; l = dn(min(p, [], 1)); h = up(max(p, [], 1));
end
function [l, h] = iadd(a, b, c, d), l = dn(a + c); h = up(b + d); end
function [l, h] = isub(a, b, c, d), l = dn(a - d); h = up(b - c); end
function [l, h] = iscal(k, a, b)              % k * [a, b], k a double constant
if k >= 0, l = dn(k*a); h = up(k*b); else, l = dn(k*b); h = up(k*a); end
end
function [l, h] = isin(a, b)                  % sin on [a, b] within [-pi/2, pi/2], else [-1, 1]
l = max(dn(dn(sin(a))), -1); h = min(up(up(sin(b))), 1);
bad = ~(a >= -1.5707963267 & b <= 1.5707963267); l(bad) = -1; h(bad) = 1;
end

function [lo, hi] = boxvals(L, U, P)
k1 = P.k1; k2 = P.k2; g = P.g;
z1l = L(1,:); z1h = U(1,:); z2l = L(2,:); z2h = U(2,:);
h1l = L(3,:); h1h = U(3,:); h2l = L(4,:); h2h = U(4,:); vl = L(5,:); vh = U(5,:);
[x1l, x1h] = iadd(z1l, z1h, vl, vh);
[k1zl, k1zh] = iscal(k1, z1l, z1h);
[k2zl, k2zh] = iscal(k2, z2l, z2h);
if isfield(P, 'plant') && strcmp(P.plant, 'arm')
    % x2 = z2 - k1 z1, u = (-z1 - k2 z2 + (h1 - k1) x2 + 10 h2 sin x1)/b
    [x2l, x2h] = isub(z2l, z2h, k1zl, k1zh);
    [a1l, a1h] = isub(h1l, h1h, k1, k1);
    [p1l, p1h] = imul(a1l, a1h, x2l, x2h);
    [snl, snh] = isin(x1l, x1h);
    [p2l, p2h] = imul(h2l, h2h, snl, snh); [p2l, p2h] = iscal(10, p2l, p2h);
    [sl, sh] = isub(-z1h, -z1l, k2zl, k2zh);
    [sl, sh] = iadd(sl, sh, p1l, p1h); [sl, sh] = iadd(sl, sh, p2l, p2h);
    ul = dn(sl/P.b); uh = up(sh/P.b);              % P.b > 0
    lo = [x1l; x2l; ul]; hi = [x1h; x2h; uh];
    return
end
[sl, sh] = isq(x1l, x1h);                          % x1^2
[q4l, q4h] = isq(sl, sh);                          % x1^4
[hsl, hsh] = imul(h1l, h1h, sl, sh);               % h1 x1^2
[wl, wh] = isub(z2l, z2h, k1zl, k1zh);             % z2 - k1 z1
[x2l, x2h] = isub(wl, wh, hsl, hsh);               % x2 = z2 - k1 z1 - h1 x1^2
[hxl, hxh] = imul(h1l, h1h, x1l, x1h);
[Al, Ah] = iscal(-2, hxl, hxh); Al = dn(Al - k1); Ah = up(Ah - k1);   % A = -k1 - 2 h1 x1
[t1l, t1h] = imul(Al, Ah, wl, wh);                 % A (z2 - k1 z1)
[snl, snh] = isin(x1l, x1h);
[t2l, t2h] = imul(h2l, h2h, snl, snh);             % h2 sin x1
[azl, azh] = imul(Al, Ah, z2l, z2h);
[el, eh] = isub(z1l, z1h, azl, azh);               % z1 - A z2
[t3l, t3h] = imul(q4l, q4h, el, eh); [t3l, t3h] = iscal(g, t3l, t3h);
[ul, uh] = isub(-z1h, -z1l, k2zl, k2zh);           % -z1 - k2 z2
[ul, uh] = iadd(ul, uh, t1l, t1h);
[ul, uh] = isub(ul, uh, t2l, t2h);
[ul, uh] = isub(ul, uh, t3l, t3h);
lo = [x1l; x2l; ul]; hi = [x1h; x2h; uh];
end
