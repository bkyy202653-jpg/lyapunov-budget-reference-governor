function reproduce_manipulator(mode)
% Single-link manipulator with unknown payload.
%   reproduce_manipulator('quick')  (default, ~5 min): re-evaluates the deployed designs (LBG, CE, ERG-N, PF+,
%       BLF) for p = 1 and the LBG for p = 2 and compares them with the reference results (eval_arm_deployed.m).
%   reproduce_manipulator('full')   (~1-2 h with 4 workers): rebuilds the certified tables and repeats the tuning
%       run_arm_p1.m     friction known, payload adapted (paper results)
%       run_arm.m        friction also unknown (p = 2, illustration of Theorem 2(c))
%       diag_arm_stall.m stalled runs end with s = Gamma_c(v_inf) (Theorem 2(a))
%       arm_gain_scan.m  stall rate for other (k, g); each (k, g) is recertified
%   and then runs the quick evaluation on the new results.
% Output: <root>/output/manipulator_<mode>/.
if nargin < 1, mode = 'quick'; end
here = fileparts(mfilename('fullpath')); addpath(genpath(fullfile(fileparts(here), 'code')));
old = pwd; back = onCleanup(@() cd(old));
switch mode
    case 'quick'
        lbg_workdir('manipulator_quick', {'certified_tables/Gamma_cert_arm_p1.mat', 'certified_tables/Gamma_cert_arm_p1_nom.mat', ...
            'certified_tables/Gamma_cert_arm.mat', 'reference_results/arm_p1.mat', 'reference_results/arm.mat'});
    case 'full'
        lbg_workdir('manipulator_full', {});
        run_script('run_arm_p1', 'run_arm_p1');
        run_script('run_arm', 'run_arm');
        run_script('diag_arm_stall', 'diag_arm_stall');
        run_script('arm_gain_scan', 'arm_gain_scan');
    otherwise
        error('mode must be ''quick'' or ''full''');
end
run_script('eval_arm_deployed', 'eval_arm_deployed');
fprintf('Logs: %s\n', fullfile(pwd, 'logs'));
end
