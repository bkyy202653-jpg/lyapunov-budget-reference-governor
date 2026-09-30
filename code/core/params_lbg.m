function P = params_lbg()
% Parameters of the academic example (n = 2); the scripts override c, umax, k1, k2 and g (see docs/REPRODUCIBILITY.md).
% Plant: x1' = x2 + th1*x1^2,  x2' = sat(u) + th2*sin(x1)
P.th0 = [0.5; 1.0];     % centre of known parameter ball Theta
P.R   = 0.25;           % radius of Theta
P.c   = [1.2; 1.5];     % state bounds |x1|<=c1, |x2|<=c2
P.umax = 5;             % input bound
P.r   = 1.0;            % setpoint
P.k1  = 1;  P.k2 = 1;   % backstepping gains
P.g   = 5;              % adaptation gain Gamma = g*I
P.margin = 0.02;        % relative constraint margin used in Gamma(v)
P.eta = 0.01;           % governor direction smoothing
P.dt  = 2e-3;  P.T = 30;
P.tol = 0.02;           % settling band |x1-r| <= tol*|r|
end
