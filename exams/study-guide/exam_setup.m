% EXAM_SETUP  Nap bo ham toolkit/ va toan bo chuong trinh cua 18 cau bai tap lon
% (solutions/qNN-*) vao duong dan. Chay mot lan dau gio thi:
%   >> cd exams/study-guide
%   >> exam_setup
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'toolkit'));
S = fullfile(here, '..', '..', 'solutions');
d = dir(fullfile(S, 'q*'));
for k = 1:numel(d)
  if d(k).isdir, addpath(fullfile(S, d(k).name)); end
end
fprintf('Da nap toolkit/ va %d thu muc solutions/q*\n', numel(d));
