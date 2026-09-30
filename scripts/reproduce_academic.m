function reproduce_academic(mode)
% Academic example: Table 1, transfer setpoints, CE ablation, tighter actuator (S1) and larger uncertainty (S2).
%   reproduce_academic('quick')  (default, ~10 min): re-evaluates the deployed designs of every method on the
%       envelope and on the transfer setpoints and compares them with the reference results
%       (eval_academic_deployed.m). The tuning searches are skipped.
%   reproduce_academic('full')   (~15 h; run_envelope.m alone took 11 h single thread): repeats the tuning searches
%       run_envelope.m   PF+ grid (432 designs)                     -> results/envelope_cert_settle.mat
%       run_baselines2.m BLF grid (288), ERG-N, CE, Table 1, transfer -> results/baselines2.mat
%       run_stress_scen  S1 (umax = 1.8) and S2 (R = 0.4): recertification, BLF grid (162), transfer
%       check_k500.m     integration-step check for a large governor gain (kappa = 500)
%   and then runs the quick evaluation on the new results. The main table Gamma_cert.mat is taken from
%   results/certified_tables (rebuild it with reproduce_certificate('full')).
% Output: <root>/output/academic_<mode>/.
if nargin < 1, mode = 'quick'; end
here = fileparts(mfilename('fullpath')); addpath(genpath(fullfile(fileparts(here), 'code')));
old = pwd; back = onCleanup(@() cd(old));
switch mode
    case 'quick'
        lbg_workdir('academic_quick', {'certified_tables/Gamma_cert.mat', 'certified_tables/Gamma_cert_S1.mat', ...
            'certified_tables/Gamma_cert_S2.mat', 'reference_results/envelope_cert_settle.mat', ...
            'reference_results/baselines2.mat', 'reference_results/stress_S1.mat', 'reference_results/stress_S2.mat'});
    case 'full'
        lbg_workdir('academic_full', {'certified_tables/Gamma_cert.mat', 'reference_results/tight_lbgtune.mat'});
        run_script('run_envelope', @() run_envelope([-0.6 -0.3 0.3 0.5 0.7], '_cert_settle', true));
        run_script('run_baselines2', 'run_baselines2');
        run_script('stress_S1', @() run_stress_scen('S1', 1.8, 0.25));
        run_script('stress_S2', @() run_stress_scen('S2', 2.5, 0.40));
        run_script('check_k500', 'check_k500');
    otherwise
        error('mode must be ''quick'' or ''full''');
end
run_script('eval_academic_deployed', 'eval_academic_deployed');
fprintf('Logs: %s\n', fullfile(pwd, 'logs'));
end
