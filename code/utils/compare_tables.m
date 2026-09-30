function ok = compare_tables(names, refsub)
% Compare the threshold tables results/<name>.mat in the current folder with the shipped reference tables in
% <root>/results/<refsub>/<name>.mat (fields v, G, b of Gt). Prints one line per table.
if nargin < 2, refsub = 'certified_tables'; end
ok = true;
for i = 1:numel(names)
    A = load(fullfile('results', [names{i} '.mat'])); B = load(fullfile(lbg_root(), 'results', refsub, [names{i} '.mat']));
    d = max(abs(A.Gt.G - B.Gt.G)./max(B.Gt.G, realmin));
    same = isequal(A.Gt.v, B.Gt.v) && isequal(A.Gt.G, B.Gt.G) && isequal(A.Gt.b, B.Gt.b);
    ok = ok && same;
    fprintf('%-22s cells %3d | identical to reference: %d | max rel. difference of G %.1e\n', names{i}, numel(A.Gt.b), same, d);
end
end
