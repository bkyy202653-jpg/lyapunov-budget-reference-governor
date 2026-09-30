% Symbolic check of the key identity for the n=2 LBG pre-stabilizer:
%   Vdot = -k1 z1^2 - k2 z2^2 + dVz/dv * vdot,  dVz/dv = -(z1 + k1 z2)
% System: x1' = x2 + th1*x1^2,  x2' = u + th2*sin(x1)
syms x1 x2 th1 th2 h1 h2 v vd k1 k2 g real
th  = [th1; th2];  thh = [h1; h2];  G = g*eye(2);
phi1 = [x1^2; 0];  phi2 = [0; sin(x1)];

z1  = x1 - v;
a1  = -k1*z1 - phi1.'*thh;
tau1 = phi1*z1;
da1dx1 = diff(a1, x1);  da1dth = [diff(a1,h1), diff(a1,h2)];  da1dv = diff(a1, v);
z2  = x2 - a1;
w2  = phi2 - da1dx1*phi1;
tau2 = tau1 + w2*z2;
u   = -z1 - k2*z2 - w2.'*thh + da1dx1*x2 + da1dth*(G*tau2);
thhdot = G*tau2;

x1dot = x2 + phi1.'*th;  x2dot = u + phi2.'*th;
% total derivatives along the closed loop, v time-varying with rate vd
z1dot = x1dot - vd;
z2dot = x2dot - (da1dx1*x1dot + da1dth*thhdot + da1dv*vd);
tt = th - thh;
Vdot = z1*z1dot + z2*z2dot - tt.'*(G\thhdot);
claim = -k1*z1^2 - k2*z2^2 - (z1 + k1*z2)*vd;
res = simplify(expand(Vdot - claim));
fprintf('Vdot - claim = %s\n', char(res));
fprintf('da1/dv = %s\n', char(da1dv));
