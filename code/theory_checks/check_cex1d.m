% Reproduce the 1-D counterexample (Proposition 2 and Appendix F of the paper): theta = theta0 = 0, target r = -3.
F = 0.005; epsS = 1e-12; r = -3; kap = 1; Dm = 1; eta = 0.1;
S = @(t) (t<=0).*0 + (t>=1).*1 + (t>0 & t<1).*(exp(-1./max(t,1e-300))./(exp(-1./max(t,1e-300)) + exp(-1./max(1-t,1e-300))));
psi = @(x) S(4*x - 1);
Gc = @(v) (v <= -2).*(F + epsS + (0.605 - F - epsS)*(v + 3)) + (v > -2)*0.605;
% budget from the exact identity s = F + z^2/2 + h^2/2 (theta = theta0, h(0) = 0), RK4, dt = 1e-3
f = @(y) rhs(y, psi, Gc, r, kap, Dm, eta);
dt = 1e-3; N = 100000; y = [1; 0; 1]; Y = zeros(N+1, 4); t = (0:N)'*dt;
for i = 1:N+1
    Y(i,:) = [y' F + 0.5*(y(1)-y(3))^2 + 0.5*y(2)^2];
    if i > N, break; end
    a = f(y); b = f(y + dt/2*a); c = f(y + dt/2*b); d = f(y + dt*c); y = y + dt/6*(a + 2*b + 2*c + d);
end
iT = find(Y(:,3) <= -2, 1);
fprintf('v reaches -2 at T = %.3f (Appendix F: T < 6.067); h(T) = %.3e (Appendix F: h(T) >= 5e-6)\n', t(iT), Y(iT,2));
fprintf('final: v = %.6f, h = %.3e, s - F = %.3e, 0.5 h^2 = %.3e, Gc(-3) - F = %.1e\n', Y(end,3), Y(end,2), Y(end,4) - F, 0.5*Y(end,2)^2, epsS);
fprintf('min over t>=T of s - Gc(v): %.3e (must stay <= 0), v stalls above -3: %d\n', max(Y(iT:end,4) - Gc(Y(iT:end,3))), Y(end,3) > -3 + 1e-6);
function dy = rhs(y, psi, Gc, r, kap, Dm, eta)
x = y(1); h = y(2); v = y(3); z = x - v; s = 0.005 + 0.5*z^2 + 0.5*h^2;
D = Gc(v) - s; vd = kap*min(Dm, max(0, D))*(r - v)/max(abs(r - v), eta);
u = -2*z - h*psi(x);
dy = [u; psi(x)*z; vd];
end
