% RUN_EXAMPLES  Cac vi du minh hoa quy hoach loi: dieu kien chinh quy (Slater)
%               va dieu kien KKT (cau 13).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q04-kkt-conditions'));    % kkt_check.m, solve_nlp.m, draw_arrow.m
addpath(fullfile('..', 'q05-lagrange-duality'));  % solve_lp.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: chieu mot diem len tap loi - KKT la dieu kien can va du
%   min (x1-3)^2 + (x2-2)^2   s.t.  g1 = x1^2 + x2^2 - 5 <= 0,  g2 = x1 + 2 x2 - 4 <= 0
fprintf('==================== VI DU 1 ====================\n');
P1.f     = @(x) (x(1) - 3)^2 + (x(2) - 2)^2;
P1.gradf = @(x) [2 * (x(1) - 3); 2 * (x(2) - 2)];
P1.g     = @(x) [x(1)^2 + x(2)^2 - 5; x(1) + 2 * x(2) - 4];
P1.jacg  = @(x) [2 * x(1), 2 * x(2); 1, 2];
[ok1, xh1, s1] = slater_check(P1, [0; 0]);
fprintf('Slater: %d, s* = %.4g, diem Slater %s, g(xhat) = %s\n', ok1, s1, mat2str(xh1', 4), mat2str(P1.g(xh1)', 4));
fprintf('--- Diem x = (2, 1) (giao cua duong tron va duong thang) ---\n');
r1 = kkt_check(P1, [2; 1]);
fprintf('--- Diem x = (-0.4, 2.2) (giao diem thu hai) ---\n');
r1b = kkt_check(P1, [-0.4; 2.2]);
[xs1, fs1] = solve_nlp(P1, [0; 0]);
fprintf('Solver: x* = %s, f* = %.6g\n', mat2str(xs1', 6), fs1);
% Dieu kien bien phan: grad f(x*)'(x - x*) >= 0 voi moi x thuoc D (luoi diem)
[G1, G2] = meshgrid(linspace(-2.5, 2.5, 501));
in = (G1.^2 + G2.^2 <= 5) & (G1 + 2 * G2 <= 4);
vi = 2 * (2 - 3) * (G1(in) - 2) + 2 * (1 - 2) * (G2(in) - 1);
fprintf('min grad f(x*)''(x - x*) tren %d diem luoi cua D = %.4g (>= 0)\n', nnz(in), min(vi) + 0);
% Chan tren cua nhan tu nho diem Slater xhat = (0, 0)
xh = [0; 0];
bnd = (P1.f(xh) - fs1) / min(-P1.g(xh));
fprintf('Chan nhan tu: sum lambda <= (f(xhat) - f*)/min(-g(xhat)) = (%g - %g)/%g = %.4g; thuc te sum lambda = %.4g\n\n', ...
  P1.f(xh), fs1, min(-P1.g(xh)), bnd, sum(r1.lambda));

%% Vi du 2: Slater khong thoa - tap chap nhan duoc chi gom mot diem
%   g1 = (x1-1)^2 + x2^2 - 1 <= 0,  g2 = (x1+1)^2 + x2^2 - 1 <= 0  =>  D = {(0, 0)}
fprintf('==================== VI DU 2 ====================\n');
P2.g    = @(x) [(x(1) - 1)^2 + x(2)^2 - 1; (x(1) + 1)^2 + x(2)^2 - 1];
P2.jacg = @(x) [2 * (x(1) - 1), 2 * x(2); 2 * (x(1) + 1), 2 * x(2)];
[ok2, xh2, s2] = slater_check(P2, [0.3; 0.4]);
fprintf('Slater: %d, s* = %.2e (tai %s)\n', ok2, s2, mat2str(xh2' + 0, 3));
fprintf('--- (a) f = x2: x* = (0, 0) la nghiem duy nhat nhung KHONG la diem KKT ---\n');
P2a = P2;  P2a.f = @(x) x(2);  P2a.gradf = @(x) [0; 1];
r2a = kkt_check(P2a, [0; 0]);
fprintf('--- (b) f = x1: x* = (0, 0) la diem KKT, nhan tu khong duy nhat, khong bi chan ---\n');
P2b = P2;  P2b.f = @(x) x(1);  P2b.gradf = @(x) [1; 0];
r2b = kkt_check(P2b, [0; 0]);
for l2 = [0 1 10]
  fprintf('   lambda = (%g, %g): ||grad_x L|| = %.1e\n', l2 + 0.5, l2, ...
    norm([1; 0] + (l2 + 0.5) * [-2; 0] + l2 * [2; 0]));
end
fprintf('\n');

%% Vi du 3: Slater dang tinh - rang buoc affine khong can thoa chat
%   min (x1-2)^2 + (x2-2)^2  s.t.  g1 = x1 + x2 - 1 <= 0,  g2 = 1 - x1 - x2 <= 0,
%                                 g3 = x1^2 + x2^2 - 1 <= 0
fprintf('==================== VI DU 3 ====================\n');
P3.f     = @(x) (x(1) - 2)^2 + (x(2) - 2)^2;
P3.gradf = @(x) [2 * (x(1) - 2); 2 * (x(2) - 2)];
P3.g     = @(x) [x(1) + x(2) - 1; 1 - x(1) - x(2); x(1)^2 + x(2)^2 - 1];
P3.jacg  = @(x) [1, 1; -1, -1; 2 * x(1), 2 * x(2)];
[ok3a, ~, s3a] = slater_check(P3, [0; 0]);
[ok3b, xh3, s3b] = slater_check(P3, [0; 0], [true; true; false]);
fprintf('Slater dang chat (moi g_i < 0): %d, s* = %.3g\n', ok3a, s3a + 0);
fprintf('Slater dang tinh (g1, g2 affine): %d, s* = %.3g, diem %s, g(xhat) = %s\n', ok3b, s3b, ...
  mat2str(xh3', 4), mat2str(P3.g(xh3)' + 0, 4));
r3 = kkt_check(P3, [0.5; 0.5]);
[xs3, fs3] = solve_nlp(P3, [0; 0]);
fprintf('Solver: x* = %s, f* = %.6g\n', mat2str(xs3', 6), fs3);
fprintf('Moi lambda = (t + 3, t, 0), t >= 0, deu la nhan tu: t = 5 cho ||grad_x L|| = %.1e\n\n', ...
  norm(P3.gradf([0.5; 0.5]) + P3.jacg([0.5; 0.5])' * [8; 5; 0]));

%% Vi du 4: phan bo cong suat (water-filling), n = 8 bien
%   min -sum log(alpha_i + x_i)  s.t.  x >= 0, sum x = P
fprintf('==================== VI DU 4 ====================\n');
alpha = [0.3; 0.8; 1.2; 0.5; 2.0; 0.9; 0.4; 1.5];
Ptot = 2;
[x4, nu4, lam4, L4] = water_filling(alpha, Ptot);
fprintf('alpha  = %s, P = %g\n', mat2str(alpha'), Ptot);
fprintf('Muc nuoc L = 1/nu = %.6g, nu = %.6g\n', L4, nu4);
fprintf('x*     = %s, sum = %.6g\n', mat2str(x4', 6), sum(x4));
fprintf('lambda = %s\n', mat2str(lam4', 4));
fprintf('Kiem tra KKT: min lambda = %.3g, max |lambda_i x_i| = %.1e, ||grad L|| = %.1e\n', ...
  min(lam4), max(abs(lam4 .* x4)), norm(-1 ./ (alpha + x4) - lam4 + nu4));
P4.f = @(x) -sum(log(alpha + x));
P4.g = @(x) -x;
P4.h = @(x) sum(x) - Ptot;
[xs4, fs4] = solve_nlp(P4, Ptot / 8 * ones(8, 1));
fprintf('Solver: x* = %s, f* = %.8g (KKT: %.8g), ||x_solver - x_KKT|| = %.1e\n\n', ...
  mat2str(xs4', 5), fs4, P4.f(x4), norm(xs4 - x4));

%% Vi du 5: bai toan loi khong tron - KKT voi duoi vi phan
%   min ||x||_1  s.t.  a'x = 1,  a = (1, -3, 2)
fprintf('==================== VI DU 5 ====================\n');
a5 = [1; -3; 2];
[~, k5] = max(abs(a5));
x5 = zeros(3, 1);  x5(k5) = 1 / a5(k5);
mu5 = -sign(x5(k5)) / a5(k5);
fprintf('Giai tay: x* = %s, f* = %.6g, mu = %.6g\n', mat2str(x5', 4), norm(x5, 1), mu5);
fprintf('Can co -mu*a_i thuoc duoi vi phan cua |x_i|: -mu*a = %s (|.| <= 1 khi x_i = 0, = sign(x_i) khi x_i ~= 0)\n', ...
  mat2str(-mu5 * a5', 4));
% Doi chieu bang quy hoach tuyen tinh: bien (x, t), min sum t, -t <= x <= t, a'x = 1
c5 = [0; 0; 0; 1; 1; 1];
A5 = [eye(3), -eye(3); -eye(3), -eye(3)];
[z5, f5, st5] = solve_lp(c5, A5, zeros(6, 1), [a5', 0, 0, 0], 1);
fprintf('Quy hoach tuyen tinh: x = %s, f* = %.6g (status %d)\n\n', mat2str(z5(1:3)' + 0, 4), f5, st5);

%% Hinh 1: vi du 1 - hinh hoc cua KKT
fig = figure('visible', 'off', 'Position', [100 100 700 620]);
hold on; grid on; box on; axis equal;
t = linspace(0, 2 * pi, 400);
% mien chap nhan duoc: hinh tron giao nua mat phang
[Xg, Yg] = meshgrid(linspace(-2.5, 3.5, 600));
a1 = atan2(1, 2);  a2 = atan2(2.2, -0.4);      % goc cua (2, 1) va (-0.4, 2.2)
arc = linspace(a2, a1 + 2 * pi, 400);             % cung tron nam trong nua mat phang x1 + 2x2 <= 4
px = [sqrt(5) * cos(arc), -0.4];  py = [sqrt(5) * sin(arc), 2.2];
hD = fill(px, py, [0.85 0.92 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 1.3);
contour(Xg, Yg, (Xg - 3).^2 + (Yg - 2).^2, [0.5 2 4.5 8 12.5 18], 'LineColor', [0.6 0.6 0.6]);
plot(sqrt(5) * cos(t), sqrt(5) * sin(t), ':', 'Color', [0.2 0.4 0.8]);
plot([-2.5 3.5], [(4 + 2.5) / 2, (4 - 3.5) / 2], ':', 'Color', [0.2 0.4 0.8]);
xs = [2; 1];
hf = draw_arrow(xs, -0.5 * P1.gradf(xs), [0.85 0 0], 0.15);
hg1 = draw_arrow(xs, 0.3 * [4; 2], [0 0.5 0], 0.15);
draw_arrow(xs, 0.6 * [1; 2], [0 0.5 0], 0.15);
plot(3, 2, 'ks', 'MarkerFaceColor', 'k');
plot(2, 1, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0]);
plot(-0.4, 2.2, 'ko', 'MarkerFaceColor', 'w');
plot(0, 0, 'k^', 'MarkerFaceColor', [0.5 0.5 0.5]);
text(3.08, 2.15, '(3, 2)', 'FontSize', 10);
text(1.75, 0.68, 'x^* = (2, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
text(-1.25, 2.55, '(-0.4, 2.2): \lambda_1 < 0', 'FontSize', 10, 'BackgroundColor', 'w');
text(0.08, -0.2, 'diem Slater (0, 0)', 'FontSize', 10);
text(3.2, 1.65, '\nabla g_1', 'FontSize', 10, 'Color', [0 0.5 0]);
text(2.3, 2.3, '\nabla g_2', 'FontSize', 10, 'Color', [0 0.5 0]);
text(2.75, 1.05, '-\nabla f', 'FontSize', 10, 'Color', [0.85 0 0]);
axis([-2.5 3.5 -2.5 3]);
xlabel('x_1'); ylabel('x_2');
title('Vi du 1: -\nabla f(x^*) = \lambda_1 \nabla g_1 + \lambda_2 \nabla g_2, \lambda = (1/3, 2/3)');
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_projection.png'));
close(fig);

%% Hinh 2: vi du 2 - hai hinh tron tiep xuc
fig = figure('visible', 'off', 'Position', [100 100 700 460]);
hold on; grid on; box on; axis equal;
h1 = plot(1 + cos(t), sin(t), 'b-', 'LineWidth', 1.3);
plot(-1 + cos(t), sin(t), 'b-', 'LineWidth', 1.3);
plot(0, 0, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0]);
draw_arrow([0; 0], [-0.8; 0], [0 0.5 0], 0.12);
draw_arrow([0; 0], [0.8; 0], [0 0.5 0], 0.12);
draw_arrow([0; 0], [0; -0.8], [0.85 0 0], 0.12);
text(-0.1, 0.17, '\nabla g_1 = (-2, 0)', 'FontSize', 10, 'Color', [0 0.5 0], 'HorizontalAlignment', 'right');
text(0.1, 0.17, '\nabla g_2 = (2, 0)', 'FontSize', 10, 'Color', [0 0.5 0]);
text(0.05, -1.1, '-\nabla f = (0, -1)  (f = x_2)', 'FontSize', 10, 'Color', [0.85 0 0]);
text(0.6, 0.6, 'g_1 \leq 0', 'FontSize', 11);
text(-1.35, 0.6, 'g_2 \leq 0', 'FontSize', 11);
axis([-2.2 2.2 -1.25 1.2]);
xlabel('x_1'); ylabel('x_2');
title('Vi du 2: D = \{(0, 0)\}, Slater khong thoa');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_no_slater.png'));
close(fig);

%% Hinh 3: vi du 4 - water-filling
fig = figure('visible', 'off', 'Position', [100 100 800 420]);
hold on; box on;
hb = bar(1:8, [alpha, x4], 'stacked');
set(hb(1), 'FaceColor', [0.6 0.6 0.6]);
set(hb(2), 'FaceColor', [0.3 0.6 1]);
hl = plot([0.4 8.6], [L4 L4], 'r--', 'LineWidth', 1.5);
text(8.7, L4, sprintf('L = %.2f', L4), 'Color', 'r', 'FontSize', 10);
for i = 1:8
  if x4(i) > 0.15
    text(i, alpha(i) + x4(i) / 2, sprintf('%.3f', x4(i)), 'HorizontalAlignment', 'center', 'FontSize', 9);
  elseif x4(i) > 0
    text(i, L4 + 0.08, sprintf('%.3f', x4(i)), 'HorizontalAlignment', 'center', 'FontSize', 9);
  else
    text(i, alpha(i) + 0.08, sprintf('\\lambda = %.3f', lam4(i)), 'HorizontalAlignment', 'center', 'FontSize', 9);
  end
end
xlim([0.4 10]); ylim([0 2.3]);
set(gca, 'XTick', 1:8);
xlabel('kenh i'); ylabel('\alpha_i + x_i');
legend([hb(1), hb(2), hl], {'\alpha_i', 'x_i^*', 'muc nuoc'}, 'Location', 'northwest');
title(sprintf('Vi du 4: water-filling, P = %g', Ptot));
print(fig, '-dpng', '-r150', fullfile('figures', 'example4_water_filling.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
