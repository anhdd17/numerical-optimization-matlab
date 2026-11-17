% RUN_EXAMPLES  Cac vi du minh hoa thuat toan don hinh hai pha (cau 2).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q01-simplex'));       % simplex.m, check_solution.m
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

%% Vi du 1: can bien gia, co phuong an toi uu
%   min f = 4x1 + x2
%   3x1 +  x2 - x3           = 3
%   4x1 + 3x2      - x4      = 6
%    x1 + 2x2           + x5 = 4
%   x >= 0
fprintf('==================== VI DU 1 ====================\n');
c1 = [4; 1; 0; 0; 0];
A1 = [3 1 -1  0 0;
      4 3  0 -1 0;
      1 2  0  0 1];
b1 = [3; 6; 4];
[x1, f1, flag1, out1] = two_phase_simplex(c1, A1, b1);
fprintf('x* = %s,  f* = %g,  exitflag = %d,  so buoc lap = %d + %d\n\n', ...
  mat2str(x1', 4), f1, flag1, out1.phase1.iterations, out1.phase2.iterations);
check_solution(c1, A1, b1, f1, isOctave);

%% Vi du 2: bai toan khong co phuong an chap nhan duoc
%   min f = x1 + x2
%   x1 + 2x2 - x3      = 6
%   x1 +  x2      + x4 = 2
%   x >= 0
fprintf('==================== VI DU 2 ====================\n');
c2 = [1; 1; 0; 0];
A2 = [1 2 -1 0;
      1 1  0 1];
b2 = [6; 2];
[~, ~, flag2, out2] = two_phase_simplex(c2, A2, b2);
fprintf('exitflag = %d,  w* = %g > 0\n\n', flag2, out2.phase1.w);

%% Vi du 3: rang buoc thua (rank A = 2 < m = 3)
%   min f = x1 + 2x2 + 3x3
%    x1 + x2 +  x3 = 3
%    x1 - x2 + 2x3 = 1
%   2x1      + 3x3 = 4     (= dong 1 + dong 2)
%   x >= 0
fprintf('==================== VI DU 3 ====================\n');
c3 = [1; 2; 3];
A3 = [1  1 1;
      1 -1 2;
      2  0 3];
b3 = [3; 1; 4];
[x3, f3, flag3, out3] = two_phase_simplex(c3, A3, b3);
fprintf('x* = %s,  f* = %g,  exitflag = %d,  rang buoc bi loai: %s\n\n', ...
  mat2str(x3', 4), f3, flag3, mat2str(out3.removedRows));
check_solution(c3, A3, b3, f3, isOctave);

%% Hinh minh hoa vi du 1: duong di qua hai pha trong mat phang (x1, x2)
fig = figure('visible', 'off');
hold on; grid on; box on;
V = [0.4 1.8; 0.6 1.2; 1.5 0; 4 0];          % cac dinh cua mien chap nhan
fill(V(:, 1), V(:, 2), [0.85 0.92 1], 'EdgeColor', 'none');
xs = linspace(0, 4.5, 200);
plot(xs, 3 - 3 * xs, 'b-', xs, (6 - 4 * xs) / 3, 'r-', xs, (4 - xs) / 2, 'm-');
p = out1.xHist(:, 1:2);
n1 = out1.phase1.iterations + 1;              % so diem cua pha 1
plot(p(1:n1, 1), p(1:n1, 2), 'k--o', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
plot(p(n1:end, 1), p(n1:end, 2), 'k-o', 'LineWidth', 1.8, 'MarkerFaceColor', 'k');
plot(x1(1), x1(2), 'p', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0], ...
  'MarkerEdgeColor', 'k');
axis([0 4.5 0 3.2]);
xlabel('x_1'); ylabel('x_2');
legend('Mien chap nhan', '3x_1 + x_2 = 3', '4x_1 + 3x_2 = 6', 'x_1 + 2x_2 = 4', ...
  'Pha 1', 'Pha 2', 'Toi uu', 'Location', 'northeast');
title('Vi du 1: duong di cua thuat toan don hinh hai pha');
if ~exist('figures', 'dir'), mkdir('figures'); end
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_two_phase_path.png'));
close(fig);
fprintf('Da luu hinh figures/example1_two_phase_path.png\n');
