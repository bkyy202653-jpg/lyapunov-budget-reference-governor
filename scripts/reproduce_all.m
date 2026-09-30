function reproduce_all(mode)
% Run every reproduction step. reproduce_all('quick') (default, ~25 min with 4 workers) or reproduce_all('full')
% (~25-30 h; the single-thread PF+ grid and certificate of Gamma_cert.mat take ~11 h and ~8 h).
if nargin < 1, mode = 'quick'; end
here = fileparts(mfilename('fullpath')); addpath(here);
reproduce_certificate(mode);
reproduce_academic(mode);
reproduce_manipulator(mode);
reproduce_figures();
end
