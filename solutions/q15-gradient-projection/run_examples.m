% RUN_EXAMPLES  Cac vi du minh hoa phuong phap chieu gradient Rosen va bien the
%               metric bien thien Murtagh-Sargent (cau 15).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q04-kkt-conditions'));    % solve_nlp.m, kkt_check.m
addpath(fullfile('..', 'q05-lagrange-duality'));  % solve_qp.m
addpath(fullfile('..', 'q06-bisection-search'));  % bisection_min.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: quy hoach toan phuong, rang buoc tuyen tinh (giai tay duoc)
%   min 2x1^2 + 2x2^2 - 2x1x2 - 4x1 - 6x2
%   s.t. x1 + x2 <= 2,  x1 + 5x2 <= 5,  -x1 <= 0,  -x2 <= 0,   x0 = (0, 0)
fprintf('==================== VI DU 1 ====================\n');
P1.f    = @(x) 2*x(1)^2 + 2*x(2)^2 - 2*x(1)*x(2) - 4*x(1) - 6*x(2);
P1.grad = @(x) [4*x(1) - 2*x(2) - 4; 4*x(2) - 2*x(1) - 6];
P1.g    = @(x) [x(1) + x(2) - 2; x(1) + 5*x(2) - 5; -x(1); -x(2)];
P1.jacg = @(x) [1 1; 1 5; -1 0; 0 -1];
x01 = [0; 0];
fprintf('-- Rosen (H = I), tim kiem chinh xac --\n');
[x1r, f1r, o1r] = gradient_projection(P1, x01, struct('lineSearch', 'exact', 'display', true));
fprintf('x* = %s = (35/31, 24/31)? sai so %.1e;  f* = %.8g;  lambda = %s\n', mat2str(x1r', 8), ...
  norm(x1r - [35; 24] / 31), f1r, mat2str(o1r.lambda', 5));
fprintf('-- Murtagh-Sargent (H cap nhat BFGS), tim kiem chinh xac --\n');
[x1m, f1m, o1m] = gradient_projection(P1, x01, struct('lineSearch', 'exact', 'metric', 'bfgs', 'display', true));
fprintf('-- Rosen, Armijo (buoc thu 1, mo rong/quay lui) --\n');
[x1a, f1a, o1a] = gradient_projection(P1, x01, struct('display', true));
[xq1, fq1] = solve_qp([4 -2; -2 4], [-4; -6], P1.jacg(0), [2; 5; 0; 0], [], []);
fprintf('Solver qp: x* = %s, f* = %.8g\n', mat2str(xq1', 8), fq1);
kkt_check(struct('gradf', P1.grad, 'g', P1.g, 'jacg', P1.jacg), [35; 24] / 31);

%% Vi du 2: rang buoc phi tuyen (vi du 1 cua cau 13)
%   min (x1-3)^2 + (x2-2)^2   s.t.  x1^2 + x2^2 <= 5,  x1 + 2x2 <= 4,   x0 = (0, 0)
fprintf('==================== VI DU 2 ====================\n');
P2.f    = @(x) (x(1) - 3)^2 + (x(2) - 2)^2;
P2.grad = @(x) [2 * (x(1) - 3); 2 * (x(2) - 2)];
P2.g    = @(x) [x(1)^2 + x(2)^2 - 5; x(1) + 2 * x(2) - 4];
P2.jacg = @(x) [2 * x(1), 2 * x(2); 1, 2];
[x2, f2, o2] = gradient_projection(P2, [0; 0], struct('display', true));
fprintf('x* = %s, f* = %.8g, lambda = %s (cau 13: (2, 1), 2, (1/3, 2/3))\n\n', mat2str(x2', 8), f2, mat2str(o2.lambda', 5));

%% Vi du 3: ham Rosenbrock tren hinh tron don vi - bien cong, can phuc hoi
%   min 100 (x2 - x1^2)^2 + (1 - x1)^2   s.t.  x1^2 + x2^2 <= 1,   x0 = (0, 0)
fprintf('==================== VI DU 3 ====================\n');
P3.f    = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
P3.grad = @(x) [-400 * x(1) * (x(2) - x(1)^2) - 2 * (1 - x(1)); 200 * (x(2) - x(1)^2)];
P3.g    = @(x) x(1)^2 + x(2)^2 - 1;
P3.jacg = @(x) [2 * x(1), 2 * x(2)];
[xs3, fs3] = solve_nlp(P3, [0; 0]);
fprintf('Solver sqp: x* = %s, f* = %.10g\n', mat2str(xs3', 10), fs3);
[x3r, f3r, o3r] = gradient_projection(P3, [0; 0], struct('maxIter', 5000));
[x3m, f3m, o3m] = gradient_projection(P3, [0; 0], struct('metric', 'bfgs', 'display', true));
fprintf('%-28s %6s %8s %16s %12s %10s\n', 'Phuong phap', 'buoc', 'exitflag', 'f', '||x - x_sqp||', 'lambda');
fprintf('%-28s %6d %8d %16.10g %12.1e %10.6g\n', 'Rosen (H = I)', o3r.iterations, o3r.exitflag, f3r, norm(x3r - xs3), o3r.lambda);
fprintf('%-28s %6d %8d %16.10g %12.1e %10.6g\n', 'Murtagh-Sargent (BFGS)', o3m.iterations, o3m.exitflag, f3m, norm(x3m - xs3), o3m.lambda);
k3 = find(o3r.X(1, :).^2 + o3r.X(2, :).^2 >= 1 - 1e-8, 1);
fprintf('Rosen: %d buoc dau trong mien (W rong, tuc la gradient voi Armijo), %d buoc tren bien\n', ...
  k3 - 1, o3r.iterations - k3 + 1);
fprintf('Do lech khoi bien |x1^2 + x2^2 - 1| lon nhat sau phuc hoi (Murtagh-Sargent): %.1e\n', ...
  max(abs(sum(o3m.X(:, 11:end).^2, 1) - 1)));
r3 = kkt_check(struct('gradf', P3.grad, 'g', P3.g, 'jacg', P3.jacg), x3m);

%% Vi du 4: danh muc Markowitz n = 50 (vi du 3 cua cau 14)
fprintf('==================== VI DU 4 ====================\n');
n = 50;  i = (1:n)';
mu  = 0.05 + 0.1 * sin(i).^2;
Fm  = cos(i * (1:3)) / sqrt(n);
Sig = Fm * Fm' + diag(0.01 + 0.01 * (1 + cos(3 * i)));
P4.f = @(x) x' * Sig * x - mu' * x;   P4.grad = @(x) 2 * Sig * x - mu;
P4.g = @(x) -x;                       P4.jacg = @(x) -eye(n);
P4.h = @(x) sum(x) - 1;               P4.jach = @(x) ones(1, n);
[xq4, fq4] = solve_qp(2 * Sig, -mu, [], [], ones(1, n), 1, zeros(n, 1));
fprintf('Solver qp: f* = %.12g\n', fq4);
fprintf('%-34s %6s %8s %12s %14s %10s\n', 'Phuong phap (x0 = (1/50, ...))', 'buoc', 'exitflag', 'f - f_qp', '||x - x_qp||', 'tai san');
lab = {'Rosen, chinh xac', 'Rosen, Armijo', 'Murtagh-Sargent, chinh xac', 'Murtagh-Sargent, Armijo'};
cfg = {struct('lineSearch', 'exact'), struct(), struct('lineSearch', 'exact', 'metric', 'bfgs'), struct('metric', 'bfgs')};
for j = 1:4
  op = cfg{j};  op.maxIter = 5000;
  [x, fx, o] = gradient_projection(P4, ones(n, 1) / n, op);
  fprintf('%-34s %6d %8d %12.1e %14.1e %10d\n', lab{j}, o.iterations, o.exitflag, fx - fq4, norm(x - xq4), nnz(x > 1e-10));
  if j == 2, o4 = o; end
end
fprintf('Rosen Armijo: |W| tai k = 0, 1, 2, 5, cuoi: %s (36 rang buoc x_i >= 0 tich cuc tai nghiem)\n', ...
  mat2str([nnz(o4.X(:, 1) <= 1e-8), nnz(o4.X(:, 2) <= 1e-8), nnz(o4.X(:, 3) <= 1e-8), nnz(o4.X(:, 6) <= 1e-8), nnz(o4.X(:, end) <= 1e-8)]));
fprintf('(Frank-Wolfe, cau 14: 6183 buoc; co buoc lui: 90 buoc)\n\n');

%% Hinh 1: vi du 1
fig = figure('visible', 'off', 'Position', [100 100 720 560]);
hold on; grid on; box on;
fill([0 2 1.25 0], [0 0 0.75 1], [0.85 0.92 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 1.2);
[X1, X2] = meshgrid(linspace(-0.3, 2.3, 300), linspace(-0.3, 1.6, 300));
Z = 2 * X1.^2 + 2 * X2.^2 - 2 * X1 .* X2 - 4 * X1 - 6 * X2;
contour(X1, X2, Z, [-7.1 -6.8 -6 -5 -4 -3 -2 -1], 'LineColor', [0.6 0.6 0.6]);
ha = plot(o1a.X(1, :), o1a.X(2, :), '.--', 'Color', [0.9 0.5 0], 'LineWidth', 1.1, 'MarkerSize', 10);
hr = plot(o1r.X(1, :), o1r.X(2, :), 'r.-', 'LineWidth', 1.6, 'MarkerSize', 16);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
plot(35/31, 24/31, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0]);
text(0.05, -0.12, 'x_0 = (0, 0)', 'FontSize', 10);
text(-0.28, 1.08, 'x_1 = (0, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
text(1.15, 0.88, 'x^* = (35/31, 24/31)', 'FontSize', 10, 'BackgroundColor', 'w');
text(1.6, 0.6, 'x_1 + x_2 \leq 2', 'FontSize', 9, 'Rotation', -45);
text(0.35, 1.04, 'x_1 + 5x_2 \leq 5', 'FontSize', 9, 'Rotation', -11);
axis equal; axis([-0.3 2.3 -0.3 1.5]);
xlabel('x_1'); ylabel('x_2');
legend([hr, ha], {sprintf('Rosen, chinh xac (%d buoc)', o1r.iterations), ...
  sprintf('Rosen, Armijo (%d buoc)', o1a.iterations)}, 'Location', 'northeast');
title('Vi du 1: chieu gradient tren da dien');
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_linear.png'));
close(fig);

%% Hinh 2: vi du 3
fig = figure('visible', 'off', 'Position', [100 100 1150 460]);
subplot(1, 2, 1); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-0.2, 1.1, 300), linspace(-0.2, 1.1, 300));
Z = 100 * (X2 - X1.^2).^2 + (1 - X1).^2;
contour(X1, X2, Z, [0.05 0.2 0.5 1 2 4 8 16 30], 'LineColor', [0.65 0.65 0.65]);
t = linspace(0, pi / 2, 200);
plot(cos(t), sin(t), 'b-', 'LineWidth', 1.4);
hr = plot(o3r.X(1, :), o3r.X(2, :), 'r.-', 'LineWidth', 0.8, 'MarkerSize', 6);
hm = plot(o3m.X(1, :), o3m.X(2, :), '.-', 'Color', [0 0.5 0], 'LineWidth', 1.3, 'MarkerSize', 12);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
plot(xs3(1), xs3(2), 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(0.82, 0.72, 'x^* = (0.7864, 0.6177)', 'FontSize', 10, 'BackgroundColor', 'w');
text(0.05, 0.98, 'x_1^2 + x_2^2 = 1', 'FontSize', 10, 'Color', 'b');
axis equal; axis([-0.1 1.1 -0.1 1.1]);
xlabel('x_1'); ylabel('x_2');
legend([hr, hm], {sprintf('Rosen (%d buoc)', o3r.iterations), sprintf('Murtagh-Sargent (%d buoc)', o3m.iterations)}, ...
  'Location', 'southoutside');
title('Vi du 3: Rosenbrock tren hinh tron');
subplot(1, 2, 2);
semilogy(0:o3r.iterations, max(o3r.F - fs3, 1e-16), 'r-', 'LineWidth', 1.2); hold on; grid on;
semilogy(0:o3m.iterations, max(o3m.F - fs3, 1e-16), '.-', 'Color', [0 0.5 0], 'LineWidth', 1.3, 'MarkerSize', 12);
set(gca, 'YTick', 10.^(-14:2:0));
xlabel('buoc lap k'); ylabel('f(x_k) - f^*');
legend('Rosen (H = I)', 'Murtagh-Sargent (BFGS)', 'Location', 'south');
title('Sai so gia tri ham');
print(fig, '-dpng', '-r150', fullfile('figures', 'example3_rosenbrock_disc.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
