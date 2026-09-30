function root = lbg_root()
% Absolute path of the repository root (this file is in <root>/code/utils).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
