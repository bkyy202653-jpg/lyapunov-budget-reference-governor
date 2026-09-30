function [u, thhdot, z1, z2, dVdv] = prestab(x1, x2, h1, h2, v, P)
% Tuning-function adaptive backstepping with frozen setpoint v (vectorized).
% Returns u, parameter update, errors and dVz/dv = -(z1 + k1 z2).
k1 = P.k1; k2 = P.k2; g = P.g;
if isfield(P, 'plant') && strcmp(P.plant, 'arm')
    % x1' = x2, x2' = b u + phi2' th, phi2 = (-x2, -10 sin x1); phi1 = 0, so alpha1 = -k1 z1 has no thetahat
    z1 = x1 - v;
    z2 = x2 + k1*z1;
    u = (-z1 - k2*z2 + h1.*x2 + 10*h2.*sin(x1) - k1*x2)/P.b;
    thhdot = [-g*x2.*z2; -10*g*sin(x1).*z2];
    if isfield(P, 'fixed1') && P.fixed1, thhdot(1,:) = 0; end    % friction known: adapt payload only
    dVdv = -(z1 + k1*z2);
    return
end
z1 = x1 - v;
a1 = -k1.*z1 - h1.*x1.^2;
da1dx1 = -k1 - 2*h1.*x1;
z2 = x2 - a1;
% tau1 = phi1*z1, w2 = phi2 - da1dx1*phi1
tau2_1 = x1.^2.*z1 - da1dx1.*x1.^2.*z2;
tau2_2 = sin(x1).*z2;
w21 = -da1dx1.*x1.^2;  w22 = sin(x1);
% da1/dthhat = [-x1^2, 0]
u = -z1 - k2.*z2 - (w21.*h1 + w22.*h2) + da1dx1.*x2 - x1.^2.*(g.*tau2_1);
thhdot = [g.*tau2_1; g.*tau2_2];
dVdv = -(z1 + k1.*z2);
end
