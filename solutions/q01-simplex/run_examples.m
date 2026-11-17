% RUN_EXAMPLES  Cac vi du minh hoa thuat toan don hinh (cau 1).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

%% Vi du 1: bai toan co phuong an toi uu duy nhat
%   min f = -3x1 - 2x2
%   x1 +  x2 + x3           = 4
%   x1 + 3x2      + x4      = 6
%   x1                 + x5 = 3
%   x >= 0
fprintf('==================== VI DU 1 ====================\n');
c1 = [-3; -2; 0; 0; 0];
A1 = [1 1 1 0 0;
      1 3 0 1 0;
      1 0 0 0 1];
b1 = [4; 6; 3];
[x1, f1, flag1, out1] = simplex(c1, A1, b1);
fprintf('x* = %s,  f* = %g,  exitflag = %d,  so buoc lap = %d\n\n', ...
  mat2str(x1', 4), f1, flag1, out1.iterations);
check_solution(c1, A1, b1, f1, isOctave);

%% Vi du 2: ham muc tieu khong bi chan duoi
%   min f = -x1 - x2
%   -x1 +  x2 + x3      = 1
%    x1 - 2x2      + x4 = 2
%   x >= 0
fprintf('==================== VI DU 2 ====================\n');
c2 = [-1; -1; 0; 0];
A2 = [-1  1 1 0;
       1 -2 0 1];
b2 = [1; 2];
[x2, f2, flag2, out2] = simplex(c2, A2, b2);
fprintf('exitflag = %d (khong bi chan). Huong giam vo han d = %s\n', ...
  flag2, mat2str(out2.direction'));
fprintf('Kiem tra: A*d = %s, c''*d = %g < 0\n\n', mat2str((A2 * out2.direction)'), c2' * out2.direction);

%% Vi du 3: bai toan suy bien gay xoay vong (vi du cua Beale, 1955)
%   min f = -3/4 x4 + 20 x5 - 1/2 x6 + 6 x7
%   x1 + 1/4 x4 -  8 x5 -     x6 + 9 x7 = 0
%   x2 + 1/2 x4 - 12 x5 - 1/2 x6 + 3 x7 = 0
%   x3                  +     x6        = 1
%   x >= 0
fprintf('==================== VI DU 3 ====================\n');
c3 = [0; 0; 0; -3/4; 20; -1/2; 6];
A3 = [1 0 0 1/4  -8   -1  9;
      0 1 0 1/2 -12 -1/2  3;
      0 0 1 0     0    1  0];
b3 = [0; 0; 1];

fprintf('--- (a) Quy tac Dantzig ---\n');
[~, ~, flag3a, out3a] = simplex(c3, A3, b3, [], struct('display', false));
fprintf('exitflag = %d sau %d buoc lap (co so cuoi: %s)\n\n', ...
  flag3a, out3a.iterations, mat2str(out3a.basis));

fprintf('--- (b) Quy tac Bland ---\n');
[x3, f3, flag3b, out3b] = simplex(c3, A3, b3, [], struct('rule', 'bland'));
fprintf('x* = %s,  f* = %g,  exitflag = %d,  so buoc lap = %d\n\n', ...
  mat2str(x3', 4), f3, flag3b, out3b.iterations);
check_solution(c3, A3, b3, f3, isOctave);

%% Hinh minh hoa vi du 1: mien chap nhan va duong di cua thuat toan
fig = figure('visible', 'off');
hold on; grid on; box on;
xs = linspace(0, 4.5, 200);
fill([0 3 3 0], [0 0 1 2], [0.85 0.92 1], 'EdgeColor', 'none');
plot(xs, 4 - xs, 'b-', xs, (6 - xs) / 3, 'r-', [3 3], [0 4.5], 'm-');
p = out1.xHist(:, 1:2);
plot(p(:, 1), p(:, 2), 'k-o', 'LineWidth', 1.8, 'MarkerFaceColor', 'k');
for k = 1:size(p, 1)
  text(p(k, 1) + 0.08, p(k, 2) + 0.15, sprintf('B%d', k - 1));
end
axis([0 4.5 0 3]);
xlabel('x_1'); ylabel('x_2');
legend('Mien chap nhan', 'x_1 + x_2 = 4', 'x_1 + 3x_2 = 6', 'x_1 = 3', ...
  'Duong di don hinh', 'Location', 'northeast');
title('Vi du 1: thuat toan don hinh di qua cac dinh');
if ~exist('figures', 'dir'), mkdir('figures'); end
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_simplex_path.png'));
close(fig);
fprintf('Da luu hinh figures/example1_simplex_path.png\n');
