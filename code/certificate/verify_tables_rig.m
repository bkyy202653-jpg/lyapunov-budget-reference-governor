% Re-verify every certified table used in the paper with the outward-rounded certifier (certify_level_rig.m).
% Each cell level b_j is tested once; a cell that is not certified is shrunk by 0.1 % and retested (max 20 times).
% Rigorous tables are saved as results/<name>_rig.mat; the printed summary lists how many cells had to be shrunk.
start_pool();   % 4 workers were used for the reported timings
Pl = params_lbg(); Pl.c = [1.2; 0.8]; Pl.umax = 2.5; Pl.k1 = 1; Pl.k2 = 1; Pl.g = 20; Pl.margin = 0.02;
Pa = params_arm();
tabs = {'Gamma_cert',        Pl;
        'Gamma_cert_S1',     setf(Pl, 'umax', 1.8);
        'Gamma_cert_S2',     setf(Pl, 'R', 0.40);
        'Gamma_cert_arm_p1', setf(Pa, 'fixed1', true);
        'Gamma_cert_arm',    Pa};
for t = 1:size(tabs, 1)
    name = tabs{t,1}; P = tabs{t,2};
    S = load(['results/' name '.mat']); Gt = S.Gt; vn = Gt.v; b = Gt.b; N = numel(b);
    br = b; nshr = zeros(1, N); boxes = zeros(1, N); tic;
    parfor j = 1:N
        lev = b(j); k = 0; ok = 0; nb = 0;
        while k <= 20
            [ok, info] = certify_level_rig(vn(j), vn(j+1), lev, P, 3e6); nb = nb + info.nproc;
            if ok == 1, break; end
            lev = lev*0.999; k = k + 1;
        end
        if ok ~= 1, lev = 0; end
        br(j) = lev; nshr(j) = k; boxes(j) = nb;
    end
    g = [br(1), min(br(1:end-1), br(2:end)), br(end)];
    Gt.G = g; Gt.b = br; Gt.LG = max(abs(diff(g))./diff(vn));
    Gt.note = 'verified by certify_level_rig.m (outward rounding, sin within 1 ulp assumed)';
    save(['results/' name '_rig.mat'], 'Gt');
    fprintf('%-18s cells %3d | certified unchanged %3d | shrunk %2d (max %d steps) | failed %d | max rel. change %.2e | %.0fs\n', ...
        name, N, sum(nshr == 0), sum(nshr > 0 & br > 0), max(nshr), sum(br == 0), max(abs(br - b)./b), toc);
end

function P = setf(P, f, v), P.(f) = v; end
