% Symbolic check: Vdot = -k1 z1^2 - k2 z2^2 - (z1 + k1 z2) vdot for the arm pre-stabilizer in prestab.m
syms x1 x2 th1 th2 h1 h2 v vd k1 k2 g b real
P.plant = 'arm'; P.k1 = k1; P.k2 = k2; P.g = g; P.b = b;
[u, hd, z1, z2] = prestab(x1, x2, h1, h2, v, P);
x1d = x2; x2d = b*u - th1*x2 - 10*th2*sin(x1);
z1d = x1d - vd; z2d = x2d + k1*z1d;
Vd = z1*z1d + z2*z2d - ([th1; th2] - [h1; h2]).'*hd/g;
fprintf('residual = %s\n', char(simplify(expand(Vd - (-k1*z1^2 - k2*z2^2 - (z1 + k1*z2)*vd)))));
