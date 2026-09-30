function G = compute_Gamma(P, vgrid, Nq, seed)
% Offline admissible level Gamma(v):
%   largest gam such that for all theta in Theta and all (z, thtilde) with
%   1/2|z|^2 + 1/2|thtilde|^2/g <= gam, the pre-stabilizer satisfies
%   |x_i| <= (1-margin) c_i and |u| <= (1-margin) umax.
% Sampled evaluation (surface + interior of the 4-D ball, Theta boundary +
% interior); returns the table and the sample set used.
if nargin < 3, Nq = 6000; end
if nargin < 4, seed = 1; end
rng(seed);
q = randn(4, Nq);  q = q ./ vecnorm(q);            % unit sphere
ns = round(0.4*Nq);                                  % interior part
q(:, 1:ns) = q(:, 1:ns) .* rand(1, ns).^(1/4);
q = [q, eye(4), -eye(4)];
% Theta samples: boundary circle + two interior rings + centre
a = linspace(0, 2*pi, 25); a(end) = [];
TH = [P.th0 + P.R*[cos(a); sin(a)], P.th0 + 0.5*P.R*[cos(a(1:2:end)); sin(a(1:2:end))], P.th0];
nq = size(q, 2); nt = size(TH, 2);
Q  = repmat(q, 1, nt);
TT = kron(TH, ones(1, nq));
lim = 1 - P.margin;
G = zeros(size(vgrid));
for j = 1:numel(vgrid)
    v = vgrid(j);
    lo = 0; hi = 5;
    for it = 1:30
        gam = 0.5*(lo + hi);
        sc = sqrt(2*gam);
        z1 = sc*Q(1,:); z2 = sc*Q(2,:);
        tt = sc*sqrt(P.g)*Q(3:4,:);
        h = TT - tt;
        x1 = z1 + v;
        a1 = -P.k1*z1 - h(1,:).*x1.^2;
        x2 = z2 + a1;
        u = prestab(x1, x2, h(1,:), h(2,:), v, P);
        worst = max([abs(x1)/P.c(1); abs(x2)/P.c(2); abs(u)/P.umax], [], 'all');
        if worst <= lim, lo = gam; else, hi = gam; end
    end
    G(j) = lo;
end
end
