% RUN_EXAMPLES  Cac vi du minh hoa phuong phap ham phat - bai toan M (cau 3).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q01-simplex'));       % simplex.m, check_solution.m
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

%% Vi du 1: can bien gia, co phuong an toi uu (cung bai toan voi vi du 1 cau 2)
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
[x1, f1, flag1, out1] = big_m_simplex(c1, A1, b1);
fprintf('x* = %s,  f* = %g,  exitflag = %d,  so buoc lap = %d\n\n', ...
  mat2str(x1', 4), f1, flag1, out1.iterations);
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
[~, ~, flag2, out2] = big_m_simplex(c2, A2, b2);
fprintf('exitflag = %d,  bien gia x%d = %g > 0\n\n', flag2, out2.artificial, out2.xa);

%% Vi du 3: ham muc tieu khong bi chan duoi
%   min f = -x1 - 2x2
%    x1 - x2 - x3      = 1
%   -x1 + x2      + x4 = 2
%   x >= 0
fprintf('==================== VI DU 3 ====================\n');
c3 = [-1; -2; 0; 0];
A3 = [ 1 -1 -1 0;
      -1  1  0 1];
b3 = [1; 2];
[x3, f3, flag3, out3] = big_m_simplex(c3, A3, b3);
fprintf('exitflag = %d,  f = %g,  bien gia x%d = %g\n\n', ...
  flag3, f3, out3.artificial, out3.xa);

%% Vi du 4: chon gia tri M cu the cho vi du 1
fprintf('==================== VI DU 4 ====================\n');
fprintf('%10s | %-24s | %8s | %8s | %s\n', 'M', 'x*', 'f*', 'exitflag', 'bien gia');
for M = [0.5 1 1.5 2 10 1e3 1e8]
  [x, f, flag, out] = big_m_simplex(c1, A1, b1, struct('M', M, 'display', false));
  fprintf('%10g | %-24s | %8g | %8d | %s\n', M, mat2str(x', 4), f, flag, ...
    mat2str(out.xa', 4));
end
fprintf('\n');

%% Hinh minh hoa vi du 1: duong di cua phuong phap M trong mat phang (x1, x2)
fig = figure('visible', 'off');
hold on; grid on; box on;
V = [0.4 1.8; 0.6 1.2; 1.5 0; 4 0];          % cac dinh cua mien chap nhan
fill(V(:, 1), V(:, 2), [0.85 0.92 1], 'EdgeColor', 'none');
xs = linspace(0, 4.5, 200);
plot(xs, 3 - 3 * xs, 'b-', xs, (6 - 4 * xs) / 3, 'r-', xs, (4 - xs) / 2, 'm-');
p = out1.xHist(:, 1:2);
plot(p(:, 1), p(:, 2), 'k-o', 'LineWidth', 1.8, 'MarkerFaceColor', 'k');
for k = 1:size(p, 1)
  text(p(k, 1) + 0.08, p(k, 2) + 0.12, sprintf('B%d', k - 1));
end
plot(x1(1), x1(2), 'p', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0], ...
  'MarkerEdgeColor', 'k');
axis([0 4.5 0 3.2]);
xlabel('x_1'); ylabel('x_2');
legend('Mien chap nhan', '3x_1 + x_2 = 3', '4x_1 + 3x_2 = 6', 'x_1 + 2x_2 = 4', ...
  'Duong di (bai toan M)', 'Toi uu', 'Location', 'northeast');
title('Vi du 1: duong di cua phuong phap ham phat');
if ~exist('figures', 'dir'), mkdir('figures'); end
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_big_m_path.png'));
close(fig);
fprintf('Da luu hinh figures/example1_big_m_path.png\n');
