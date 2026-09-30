function start_pool()
% Start a parallel pool with up to 4 workers if the Parallel Computing Toolbox is installed and licensed.
% Without it, parfor loops run serially and give identical results (every parfor iteration is independent).
if isempty(ver('parallel')) || ~license('test', 'Distrib_Computing_Toolbox'), return; end
if isempty(gcp('nocreate')), parpool(min(4, feature('numcores'))); end
end
