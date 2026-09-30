function build_all_tables()
% Rebuild the certified tables of the recertification scenarios and of the manipulator from scratch, with the
% same calls as in run_stress_scen.m, run_arm_p1.m and run_arm.m (those scripts build a table only if the file
% is missing). The main academic table results/Gamma_cert.mat is built separately by build_cert_table.m (~8 h).
start_pool();
Pl = params_lbg(); Pl.c = [1.2; 0.8]; Pl.k1 = 1; Pl.k2 = 1; Pl.g = 20; Pl.margin = 0.02;
vs = round((-0.90:0.02:0.86)*100)/100;
for sc = {'S1', 1.8, 0.25; 'S2', 2.5, 0.40}'
    P = Pl; P.umax = sc{2}; P.R = sc{3};
    t = tic; build_cert_generic(P, vs, 3e5, 14, ['results/Gamma_cert_' sc{1} '.mat']);
    fprintf('Gamma_cert_%s: %.0f s\n', sc{1}, toc(t));
end
va = round((-1.20:0.02:1.20)*100)/100;
for fx = [true false]
    P = params_arm(); nm = 'arm';
    if fx, P.fixed1 = true; nm = 'arm_p1'; end
    t = tic; build_cert_generic(P, va, 3e5, 14, ['results/Gamma_cert_' nm '.mat']);
    fprintf('Gamma_cert_%s: %.0f s\n', nm, toc(t));
    Pn = P; Pn.g = 0; Pn.R = 0;                               % nominal threshold of the ERG-N baseline
    build_cert_generic(Pn, va, 3e5, 14, ['results/Gamma_cert_' nm '_nom.mat']);
end
end
