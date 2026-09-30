function reproduce_certificate(mode)
% Interval certificate and theory checks.
%   reproduce_certificate('quick')  (default, ~9 min with 4 workers)
%       - symbolic check of the exact identity (Lemma 1) for both examples      sym_check.m, sym_check_arm.m
%       - counterexample of Proposition 2                                         check_cex1d.m
%       - outward-rounded re-verification of all 591 cells of the five tables     verify_tables_rig.m
%         used in the paper, and comparison with the shipped *_rig.mat tables
%       - refutation of two known infeasible levels and tightness test            test_certify4.m
%         (multistart fmincon finds violations at 1.02 x the certified level; needs Optimization Toolbox)
%   reproduce_certificate('full')   (~9 h)
%       - rebuilds every certified table from scratch (build_cert_table.m, build_all_tables.m), compares them
%         with the shipped tables, then runs the quick checks on the rebuilt tables
% Output: <root>/output/certificate_<mode>/ (results/, logs/). The shipped files in <root>/results are not changed.
if nargin < 1, mode = 'quick'; end
here = fileparts(mfilename('fullpath')); addpath(genpath(fullfile(fileparts(here), 'code')));
old = pwd; back = onCleanup(@() cd(old));
tabs = {'Gamma_cert', 'Gamma_cert_S1', 'Gamma_cert_S2', 'Gamma_cert_arm_p1', 'Gamma_cert_arm'};
switch mode
    case 'quick'
        lbg_workdir('certificate_quick', strcat('certified_tables/', tabs, '.mat'));
    case 'full'
        lbg_workdir('certificate_full', {'reference_results/tight_lbgtune.mat'});
        run_script('build_cert_table', 'build_cert_table');                 % Gamma_cert.mat, 175 cells, ~8 h
        run_script('build_all_tables', @build_all_tables);                  % S1, S2, arm_p1, arm (+ nominal)
        okb = compare_tables([tabs, {'Gamma_cert_arm_p1_nom', 'Gamma_cert_arm_nom'}]);
        fprintf('rebuilt tables identical to the shipped ones: %d\n', okb);
    otherwise
        error('mode must be ''quick'' or ''full''');
end
run_script('sym_check', 'sym_check');
run_script('sym_check_arm', 'sym_check_arm');
run_script('check_cex1d', 'check_cex1d');
run_script('test_certify4', 'test_certify4');
run_script('verify_tables_rig', 'verify_tables_rig');
ok = compare_tables(strcat(tabs, '_rig'));
fprintf('\nCertificate check (%s): outward-rounded tables identical to the shipped *_rig.mat: %d\n', mode, ok);
fprintf('Logs: %s\n', fullfile(pwd, 'logs'));
end
