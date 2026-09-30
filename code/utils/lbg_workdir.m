function wd = lbg_workdir(name, inputs)
% Create the working folder <root>/output/<name>/ with a results/ subfolder, copy the listed reference files
% (paths relative to <root>/results) into it, and make it the current folder.
% All experiment scripts read and write 'results/<file>.mat' relative to the current folder, so the shipped
% reference results in <root>/results are never overwritten.
root = lbg_root();
wd = fullfile(root, 'output', name); rs = fullfile(wd, 'results');
if ~exist(rs, 'dir'), mkdir(rs); end
for i = 1:numel(inputs)
    src = fullfile(root, 'results', inputs{i});
    if ~exist(src, 'file'), error('lbg_workdir:missing', 'Reference file not found: %s', src); end
    copyfile(src, rs);
end
cd(wd);
end
