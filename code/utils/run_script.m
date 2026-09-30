function run_script(lbg_logname__, lbg_cmd__, figdir) %#ok<INUSD>
% Run an experiment and write its console output to logs/<logname>.log in the current folder.
% cmd: name of a script (run in the workspace of this function) or a function handle, e.g.
%      @() run_stress_scen('S1', 1.8, 0.25). figdir (optional) is visible to scripts as a variable
%      (used by make_figs.m). Local names end in '__' so that scripts cannot overwrite them.
if ~exist('logs', 'dir'), mkdir('logs'); end
lbg_lf__ = fullfile('logs', [lbg_logname__ '.log']); if exist(lbg_lf__, 'file'), delete(lbg_lf__); end
diary(lbg_lf__); lbg_cleanup__ = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('=== %s (%s) ===\n', lbg_logname__, char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'))); lbg_t0__ = tic;
if ischar(lbg_cmd__) || isstring(lbg_cmd__), run(char(lbg_cmd__)); else, lbg_cmd__(); end
fprintf('=== %s done in %.0f s ===\n', lbg_logname__, toc(lbg_t0__));
end
