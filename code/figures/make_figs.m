% Figures of the manuscript: Fig. 2 fig_gamma, Fig. 3 fig_response, Fig. 4 fig_ablation, Fig. 5 fig_transfer,
% Fig. 6 fig_arm (Fig. 1 is the TikZ block diagram in the LaTeX source). Uses the certified tables and the deployed
% designs. Output folder: variable figdir if it exists in the workspace, otherwise figures/ in the current folder.
if ~exist('figdir', 'var'), figdir = 'figures'; end
out = figdir; if ~exist(out, 'dir'), mkdir(out); end
set(groot, 'defaultAxesFontSize', 9, 'defaultLineLineWidth', 1.1, 'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', 'defaultLegendInterpreter', 'latex');
P = params_lbg(); P.c = [1.2; 0.8]; P.umax = 2.5; P.k1 = 1; P.k2 = 1; P.g = 20; P.margin = 0.02;
C = load('results/Gamma_cert.mat'); Gc = C.Gt; F = P.R^2/(2*P.g); lim = 1 - P.margin;
col = lines(4); gr = 0.45*[1 1 1]; fs = 8;
pl = @(s) title(s, 'FontWeight', 'normal', 'FontSize', 9);          % panel label

%% fig_response (Fig. 3): LBG responses, academic example, r = 0.5 from rest, three true parameters
TH = [P.th0, P.th0 + P.R*[1; 0], P.th0 + P.R*[-sqrt(0.5); sqrt(0.5)]];
Pr = P; Pr.r = 0.5; Pr.T = 15; o = simulate(Pr, TH, [0; 0], 'LBG', 20, Gc); t = o.rec.t;
thl = arrayfun(@(j) sprintf('$\\theta=(%.2f,%.2f)$', TH(1,j), TH(2,j)), 1:3, 'UniformOutput', false);
f = figure('Units', 'centimeters', 'Position', [2 2 17 10.5]);
subplot(2,2,1); hold on; h = gobjects(1, 3);
for j = 1:3, h(j) = plot(t, o.rec.x1(:,j), 'Color', col(j,:)); plot(t, o.rec.v(:,j), '--', 'Color', col(j,:)); end
yline(P.c(1), 'k:'); yline(Pr.r, 'k-.');
text(14.6, P.c(1) - 0.07, '$x_1=c_1$', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(14.6, Pr.r + 0.09, 'setpoint $r$', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(0.25, 0.78, '$v$ (dashed)', 'FontSize', fs, 'Color', gr); text(2.2, 0.12, '$x_1$ (solid)', 'FontSize', fs, 'Color', gr);
legend(h, thl, 'Location', 'southeast', 'FontSize', 7); ylabel('$x_1$, $v$'); xlabel('$t$ (s)'); box on; ylim([-0.1 1.3]); pl('(a) output and applied reference');
subplot(2,2,2); hold on; for j = 1:3, plot(t, o.rec.x2(:,j), 'Color', col(j,:)); end
yline(P.c(2), 'k:'); yline(-P.c(2), 'k:');
text(14.6, P.c(2) - 0.09, '$x_2=c_2$', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(14.6, -P.c(2) + 0.09, '$x_2=-c_2$', 'HorizontalAlignment', 'right', 'FontSize', fs);
ylabel('$x_2$'); xlabel('$t$ (s)'); box on; pl('(b) second state');
subplot(2,2,3); hold on; for j = 1:3, plot(t, o.rec.u(:,j), 'Color', col(j,:)); end
yline(P.umax, 'k:'); yline(-P.umax, 'k:');
text(14.6, P.umax - 0.25, '$u=u_{\max}$', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(14.6, -P.umax + 0.25, '$u=-u_{\max}$', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(7.5, 1.0, sprintf('$\\max|u|=%.2f$, no saturation', max(abs(o.rec.u(:)))), 'HorizontalAlignment', 'center', 'FontSize', fs);
ylabel('$u$'); xlabel('$t$ (s)'); box on; pl('(c) control input');
subplot(2,2,4); hold on; j = 2;
plot(t, o.rec.Gv(:,j), 'k'); plot(t, o.rec.s(:,j), 'Color', col(1,:)); plot(t, o.rec.V(:,j), '--', 'Color', col(2,:));
text(6, 0.03, '$\Delta=\Gamma_c(v)-s\ge0$', 'FontSize', fs);
text(6, 0.009, '$s-V=\mathrm{const}\in[0,F]$', 'FontSize', fs);
legend({'$\Gamma_c(v)$', '$s$ (budget)', '$V$ (unmeasurable)'}, 'Location', 'northeast'); ylabel('level'); xlabel('$t$ (s)'); box on;
pl(['(d) budget, ' thl{j}]);
exportgraphics(f, fullfile(out, 'fig_response.pdf'), 'ContentType', 'vector'); close(f);
fprintf('fig_response (Fig. 3): max|x2| %.3f max|u| %.3f min(Gv-s) %.2e min(s-V) %.2e\n', max(abs(o.rec.x2(:))), max(abs(o.rec.u(:))), min(o.rec.Gv(:)-o.rec.s(:)), min(o.rec.s(:)-o.rec.V(:)));

%% fig_gamma (Fig. 2): certified thresholds normalised by F (academic example S0/S1 and manipulator)
S1 = load('results/Gamma_cert_S1.mat'); A1 = load('results/Gamma_cert_arm_p1.mat'); Pa = params_arm(); Fa = Pa.R^2/(2*Pa.g);
Renv = [-0.6 -0.3 0.3 0.5 0.7]; Rarm = [-1 -0.6 -0.3 0.3 0.6 1];
f = figure('Units', 'centimeters', 'Position', [2 2 17 6]);
subplot(1,2,1); lo = [-0.9 0.85];
h1 = semilogy(Gc.v, Gc.G/F, 'Color', col(1,:)); hold on; h2 = semilogy(S1.Gt.v, S1.Gt.G/F, '--', 'Color', col(2,:));
h3 = semilogy(Renv, interp1(Gc.v, Gc.G, Renv)/F, 'o', 'Color', col(1,:), 'MarkerSize', 4, 'MarkerFaceColor', 'w');
hp = patch([lo(1) lo(2) lo(2) lo(1)], [0.1 0.1 1 1], 0.92*[1 1 1], 'EdgeColor', 'none', 'HandleVisibility', 'off'); uistack(hp, 'bottom'); yline(1, 'k:');
text(0, 0.35, '$\Gamma_c<F$: convergence not certified', 'HorizontalAlignment', 'center', 'FontSize', fs);
xlabel('$v$'); ylabel('$\Gamma_c(v)/F$'); box on; xlim(lo); ylim([0.1 1000]);
legend([h1 h2 h3], {'$u_{\max}=2.5$', '$u_{\max}=1.8$', 'tuning setpoints'}, 'Location', 'northeast');   % (FontSize 7 here breaks the first legend icon in vector export)
pl('(a) academic example');
subplot(1,2,2); lo = [-1.2 1.2];
semilogy(A1.Gt.v, A1.Gt.G/Fa, 'Color', col(3,:)); hold on;
semilogy(Rarm, interp1(A1.Gt.v, A1.Gt.G, Rarm)/Fa, 'o', 'Color', col(3,:), 'MarkerSize', 4, 'MarkerFaceColor', 'w');
hp = patch([lo(1) lo(2) lo(2) lo(1)], [0.5 0.5 1 1], 0.92*[1 1 1], 'EdgeColor', 'none', 'HandleVisibility', 'off'); uistack(hp, 'bottom'); yline(1, 'k:');
text(0, 0.75, '$\Gamma_c<F$', 'HorizontalAlignment', 'center', 'FontSize', fs);
xlabel('$v$ (rad)'); ylabel('$\Gamma_c(v)/F$'); box on; xlim(lo); ylim([0.5 1000]); pl('(b) manipulator');
exportgraphics(f, fullfile(out, 'fig_gamma.pdf'), 'ContentType', 'vector'); close(f);

%% fig_transfer (Fig. 5): transfer to untuned setpoints (S1: umax = 1.8): LBG vs envelope-tuned BLF; S0: LBG vs PF+
E = load('results/envelope_cert_settle.mat'); X0set = E.X0set;
a = linspace(0, 2*pi, 9); a(end) = []; THt = [P.th0 + P.R*[cos(a); sin(a)], P.th0];
[I, J] = ndgrid(1:size(THt,2), 1:size(X0set,2)); THb = THt(:, I(:)); X0b = X0set(:, J(:));
Rtr = [-0.9 -0.8 -0.45 0.4 0.6 0.75 0.8 0.85];
bp = E.bp; vL0 = zeros(size(Rtr)); vP0 = vL0;
for i = 1:numel(Rtr)
    Q = P; Q.r = Rtr(i); q = simulate(Q, THb, X0b, 'LBG', 20, Gc); vL0(i) = 100*mean(q.maxviolx > 1e-6);
    Q = P; Q.k1 = bp(1); Q.k2 = bp(1); Q.g = bp(2); Q.r = Rtr(i); q = simulate(Q, THb, X0b, 'PF', bp(3:4), Gc); vP0(i) = 100*mean(q.maxviolx > 1e-6);
end
T1 = load('results/stress_S1.mat'); P1 = P; P1.umax = 1.8; X1 = T1.X0set;
[I, J] = ndgrid(1:size(THt,2), 1:size(X1,2)); TH1 = THt(:, I(:)); X1b = X1(:, J(:));
R1 = T1.Rtr; vL1 = zeros(size(R1)); vB1 = vL1;
for i = 1:numel(R1)
    Q = P1; Q.r = R1(i); q = simulate(Q, TH1, X1b, 'LBG', 20, S1.Gt); vL1(i) = 100*mean(q.maxviolx > 1e-6);
    q = simulate_blf(Q, TH1, X1b, T1.bB); vB1(i) = 100*mean(q.maxviolx > 1e-6);
end
f = figure('Units', 'centimeters', 'Position', [2 2 17 6]);
subplot(1,2,1); barlab(Rtr, [vL0; vP0], {'LBG (0\% at every setpoint)', 'PF$^+$ (tuned for $u_{\max}=2.5$)'}, 30, fs); pl('(a) $u_{\max}=2.5$, 108 runs per setpoint');
subplot(1,2,2); barlab(R1, [vL1; vB1], {'LBG (0\% at every setpoint)', 'BLF (re-tuned for $u_{\max}=1.8$)'}, 55, fs); pl('(b) $u_{\max}=1.8$, 72 runs per setpoint');
exportgraphics(f, fullfile(out, 'fig_transfer.pdf'), 'ContentType', 'vector'); close(f);
fprintf('fig_transfer (Fig. 5) umax=2.5: LBG %s PF %s | umax=1.8 (r=%s): LBG %s BLF %s\n', mat2str(vL0,3), mat2str(vP0,3), mat2str(R1), mat2str(vL1,3), mat2str(vB1,3));

%% fig_arm (Fig. 6): manipulator (friction known, payload adapted), r = 1.0 rad from rest
Pa.fixed1 = true; THa = [Pa.th0(1)*[1 1 1]; Pa.th0(2) + Pa.R*[-1 0 1]];
Q = Pa; Q.r = 1.0; Q.T = 15; o = simulate(Q, THa, [0; 0], 'LBG', 20, A1.Gt); t = o.rec.t;
mgl = arrayfun(@(j) sprintf('$mgl=%.1f$ N\\,m', THa(2,j)*10*0.1), 1:3, 'UniformOutput', false);   % mgl = 10 J th2, J = 0.1
f = figure('Units', 'centimeters', 'Position', [2 2 17 10.5]);
subplot(2,2,1); hold on; h = gobjects(1, 3);
for j = 1:3, h(j) = plot(t, o.rec.x1(:,j), 'Color', col(j,:)); plot(t, o.rec.v(:,j), '--', 'Color', col(j,:)); end
yline(Pa.c(1), 'k:'); yline(Q.r, 'k-.');
text(14.6, Pa.c(1) - 0.08, 'joint limit', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(0.4, 0.78, '$v$ (dashed)', 'FontSize', fs, 'Color', gr); text(2.6, 0.12, '$q$ (solid)', 'FontSize', fs, 'Color', gr);
legend(h, mgl, 'Location', 'southeast', 'FontSize', 7); ylabel('$q$, $v$ (rad)'); xlabel('$t$ (s)'); box on; pl('(a) joint angle and applied reference');
subplot(2,2,2); hold on; for j = 1:3, plot(t, o.rec.x2(:,j), 'Color', col(j,:)); end
yline(Pa.c(2), 'k:'); yline(-Pa.c(2), 'k:');
text(14.6, Pa.c(2) - 0.12, 'speed limit', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(14.6, -Pa.c(2) + 0.12, 'speed limit', 'HorizontalAlignment', 'right', 'FontSize', fs);
ylabel('$\dot q$ (rad/s)'); xlabel('$t$ (s)'); box on; pl('(b) joint speed');
subplot(2,2,3); hold on; for j = 1:3, plot(t, o.rec.u(:,j), 'Color', col(j,:)); end
yline(Pa.umax, 'k:'); yline(-Pa.umax, 'k:');
text(14.6, Pa.umax - 0.17, 'torque limit', 'HorizontalAlignment', 'right', 'FontSize', fs);
text(14.6, -Pa.umax + 0.17, 'torque limit', 'HorizontalAlignment', 'right', 'FontSize', fs);
ylabel('$\tau$ (N\,m)'); xlabel('$t$ (s)'); box on; pl('(c) motor torque');
subplot(2,2,4); for j = 1:3, semilogy(t, o.rec.Gv(:,j) - o.rec.s(:,j), 'Color', col(j,:)); hold on; end
text(7.5, 0.05, sprintf('$\\min\\Delta=%.2f\\times10^{-3}>0$', 1e3*min(o.rec.Gv(:)-o.rec.s(:))), 'HorizontalAlignment', 'center', 'FontSize', fs);
ylabel('$\Delta=\Gamma_c(v)-s$'); xlabel('$t$ (s)'); box on; ylim([1e-3 1]); pl('(d) budget margin (log scale)');
exportgraphics(f, fullfile(out, 'fig_arm.pdf'), 'ContentType', 'vector'); close(f);
fprintf('fig_arm (Fig. 6): max|qdot| %.3f max|tau| %.3f min Delta %.2e ts %s\n', max(abs(o.rec.x2(:))), max(abs(o.rec.u(:))), min(o.rec.Gv(:)-o.rec.s(:)), mat2str(o.ts, 3));

%% fig_ablation (Fig. 4): certainty-equivalence margin Gamma_c - V_z vs budget margin, same run
bestj = []; bestr = [];
for r = [0.7 0.5 -0.6 0.3 -0.3]
    Q = P; Q.r = r; q = simulate(Q, THb, X0b, 'RERG', 20, Gc);
    [mv, j] = max(q.maxviolx); if mv > 1e-3 && isfinite(mv), bestj = j; bestr = r; break; end
end
if ~isempty(bestj)
    Q = P; Q.r = bestr; Q.T = 15;
    qc = simulate(Q, THb(:,bestj), X0b(:,bestj), 'RERG', 20, Gc); ql = simulate(Q, THb(:,bestj), X0b(:,bestj), 'LBG', 20, Gc); t = qc.rec.t;
    iv = find(abs(qc.rec.x2) > P.c(2), 1);                       % first recorded sample outside the bound
    f = figure('Units', 'centimeters', 'Position', [2 2 17 6]);
    subplot(1,2,1); plot(t, qc.rec.x2, 'Color', col(2,:)); hold on; plot(t, ql.rec.x2, 'Color', col(1,:));
    yline(P.c(2), 'k:'); yline(-P.c(2), 'k:');
    if ~isempty(iv)
        plot(t(iv), qc.rec.x2(iv), 'x', 'Color', col(2,:), 'MarkerSize', 7, 'LineWidth', 1.3, 'HandleVisibility', 'off');
        text(t(iv) + 0.4, qc.rec.x2(iv) - 0.25, sprintf('violation at $t=%.2f$ s', t(iv)), 'FontSize', fs);
    end
    text(14.6, P.c(2) + 0.15, '$x_2=\pm c_2$', 'HorizontalAlignment', 'right', 'FontSize', fs);
    ylim([-2 2]); ylabel('$x_2$'); xlabel('$t$ (s)'); legend({'CE ablation', 'LBG'}, 'Location', 'northeast', 'FontSize', 7); box on;
    pl('(a) second state');
    subplot(1,2,2); plot(t, qc.rec.v, 'Color', col(2,:)); hold on; plot(t, ql.rec.v, 'Color', col(1,:)); yline(bestr, 'k-.');
    text(14.6, bestr + 0.05, sprintf('setpoint $r=%.1f$', bestr), 'HorizontalAlignment', 'right', 'FontSize', fs);
    ylim([0 1]); ylabel('$v$'); xlabel('$t$ (s)'); box on; pl('(b) applied reference');
    exportgraphics(f, fullfile(out, 'fig_ablation.pdf'), 'ContentType', 'vector'); close(f);
    fprintf('fig_ablation (Fig. 4): r=%.2f theta=%s x0=%s CE maxviol %.3f LBG maxviol %.3f\n', bestr, mat2str(THb(:,bestj)',3), mat2str(X0b(:,bestj)',3), qc.maxviolx, ql.maxviolx);
end

function barlab(R, Y, names, ytop, fs)
% grouped bars (LBG first) with the value of every nonzero bar printed above it
b = bar(categorical(string(R), string(R)), Y'); hold on;
for k = 1:numel(b)
    for i = 1:numel(R)
        if Y(k,i) > 0, text(b(k).XEndPoints(i), Y(k,i) + 0.02*ytop, sprintf('%.1f', Y(k,i)), 'HorizontalAlignment', 'center', 'FontSize', 6.5); end
    end
end
ylim([0 ytop]); ylabel('violating runs (\%)'); xlabel('transfer setpoint $r$'); legend(names, 'Location', 'northwest', 'FontSize', 7); box on;
end
