function P = params_arm()
% Single-link manipulator with unknown viscous friction and unknown payload (gravity torque), normalised:
%   J q'' + B q' + m g l sin q = tau,  x1 = q, x2 = q', b = 1/J known,
%   x1' = x2,  x2' = b sat(u) - th1 x2 - 10 th2 sin x1,  th1 = B/J, th2 = m g l/(10 J).
% J = 0.1 kg m^2 (b = 10); th0 = (0.5, 1.0) <-> B = 0.05 N m s, m g l = 1.0 N m; ball radius R = 0.3
% (friction 0.02..0.08 N m s, gravity torque 0.7..1.3 N m, i.e. +-30 % payload uncertainty).
P.plant = 'arm';
P.b    = 10;
P.th0  = [0.5; 1.0];
P.R    = 0.3;
P.c    = [1.4; 1.2];      % joint limit |q| <= 1.4 rad (80 deg), speed limit |q'| <= 1.2 rad/s
P.umax = 1.5;             % torque limit 1.5 N m
P.r    = 0.9;
P.k1 = 1; P.k2 = 1; P.g = 20;
P.margin = 0.02;
P.eta = 0.01;
P.dt = 2e-3; P.T = 20;
P.tol = 0.02;
end
