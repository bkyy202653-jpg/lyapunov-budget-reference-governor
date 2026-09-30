% Does the adaptation gain change the stall rate on the arm? Recertify per (k, g), run LBG on 6 setpoints, T = 40.
start_pool();   % 4 workers were used for the reported timings
P0 = params_arm(); vn = round((-1.20:0.04:1.20)*100)/100;
[g1, g2] = ndgrid(-0.5:0.25:0.5, -0.6:0.3:0.6); X0c = [g1(:) g2(:)]';
a = linspace(0, 2*pi, 9); a(end) = []; Rset = [-1 -0.6 -0.3 0.3 0.6 1];
for kg = [1 2; 1 5; 1 80; 2 20; 3 20]'
    P = P0; P.k1 = kg(1); P.k2 = kg(1); P.g = kg(2); F = P.R^2/(2*P.g);
    Gt = build_cert_generic(P, vn, 2e5, 12, sprintf('results/Gamma_cert_arm_k%d_g%d.mat', kg));
    adm = false(1, 25);
    for j = 1:25
        [~,~,z1,z2] = prestab(X0c(1,j), X0c(2,j), P.th0(1), P.th0(2), X0c(1,j), P);
        adm(j) = 0.5*(z1^2 + z2^2) + F <= interp1(Gt.v, Gt.G, X0c(1,j));
    end
    X0set = X0c(:, adm); if isempty(X0set), fprintf('k=%d g=%d: no admissible x0\n', kg); continue; end
    THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
    [I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
    ts = []; vx = []; vth = [];
    for r = Rset
        Pr = P; Pr.r = r; Pr.T = 40; o = simulate(Pr, THb, X0b, 'LBG', 20, Gt);
        ts = [ts o.ts]; vx = [vx o.maxviolx]; vth = [vth o.rec.V(end,:)/F]; %#ok<AGROW>
    end
    fprintf('k=%d g=%3d: adm x0 %2d, Gam/F at +-1: %.2f, stalls %d/%d, median ts(settled) %.2f, viol %d, median V_end/F %.2f\n', ...
        kg, sum(adm), interp1(Gt.v, Gt.G, 1)/F, sum(~isfinite(ts)), numel(ts), median(ts(isfinite(ts))), sum(vx > 1e-6), median(vth));
end
