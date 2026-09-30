function reproduce_figures()
% Figures 2-6 of the paper (make_figs.m, ~5 min). Fig. 1 is a block diagram drawn in the LaTeX source.
% Output: <root>/output/figures/figures/fig_*.pdf and the printed checks in <root>/output/figures/logs/make_figs.log.
here = fileparts(mfilename('fullpath')); addpath(genpath(fullfile(fileparts(here), 'code')));
old = pwd; back = onCleanup(@() cd(old));
wd = lbg_workdir('figures', {'certified_tables/Gamma_cert.mat', 'certified_tables/Gamma_cert_S1.mat', ...
    'certified_tables/Gamma_cert_arm_p1.mat', 'reference_results/envelope_cert_settle.mat', 'reference_results/stress_S1.mat'});
run_script('make_figs', 'make_figs', fullfile(wd, 'figures'));
fprintf('Figures: %s (compare the printed values with results/logs/make_figs.log)\n', fullfile(wd, 'figures'));
end
