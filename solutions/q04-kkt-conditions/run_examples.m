% RUN_EXAMPLES  Cac vi du minh hoa dieu kien chinh quy va dieu kien KKT (cau 4).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;

%% Vi du 1: chi co rang buoc dang thuc (nhan tu Lagrange co dien)
%   min x1 + x2   s.t.  h(x) = x1^2 + x2^2 - 2 = 0
fprintf('==================== VI DU 1 ====================\n');
P1.f     = @(x) x(1) + x(2);
P1.gradf = @(x) [1; 1];
P1.h     = @(x) x(1)^2 + x(2)^2 - 2;
P1.jach  = @(x) [2*x(1), 2*x(2)];
for xc = {[1; 1], [-1; -1]}
  r1 = kkt_check(P1, xc{1});
  % Dieu kien cap 2: Hessian cua L = 2*mu*I tren khong gian tiep xuc
  if r1.mu > 0
    kind = 'xac dinh duong: cuc tieu dia phuong';
  else
    kind = 'xac dinh am: cuc dai dia phuong';
  end
  fprintf('  Hessian cua L = %g*I  =>  %s\n\n', 2 * r1.mu, kind);
end
[xs, fs] = solve_nlp(P1, [0.5; -0.2]);
fprintf('  Solver: x* = %s, f* = %g\n\n', mat2str(xs', 5), fs);

%% Vi du 2: chi co rang buoc bat dang thuc
%   min (x1-2)^2 + (x2-1)^2   s.t.  g1 = x1^2 - x2 <= 0,  g2 = x1 + x2 - 2 <= 0
fprintf('==================== VI DU 2 ====================\n');
P2.f     = @(x) (x(1) - 2)^2 + (x(2) - 1)^2;
P2.gradf = @(x) [2*(x(1) - 2); 2*(x(2) - 1)];
P2.g     = @(x) [x(1)^2 - x(2); x(1) + x(2) - 2];
P2.jacg  = @(x) [2*x(1), -1; 1, 1];
fprintf('--- Diem x = (1, 1) ---\n');
r2 = kkt_check(P2, [1; 1]);
fprintf('--- Diem x = (0, 0) (chap nhan duoc nhung khong phai diem KKT) ---\n');
kkt_check(P2, [0; 0]);
[xs, fs] = solve_nlp(P2, [0; 0]);
fprintf('  Solver: x* = %s, f* = %g\n\n', mat2str(xs', 5), fs);

%% Vi du 3: rang buoc hon hop dang thuc va bat dang thuc
%   min x1^2 + x2^2 + x3^2   s.t.  h = x1 + x2 + x3 - 3 = 0,  g = 2 - x1 <= 0
fprintf('==================== VI DU 3 ====================\n');
P3.f     = @(x) sum(x.^2);
P3.gradf = @(x) 2 * x(:);
P3.g     = @(x) 2 - x(1);
P3.jacg  = @(x) [-1, 0, 0];
P3.h     = @(x) sum(x) - 3;
P3.jach  = @(x) [1, 1, 1];
kkt_check(P3, [2; 0.5; 0.5]);
[xs, fs] = solve_nlp(P3, [0; 0; 0]);
fprintf('  Solver: x* = %s, f* = %g\n\n', mat2str(xs', 5), fs);

%% Vi du 4: MFCQ dung, LICQ sai - nhan tu khong duy nhat
%   min x1 + 2x2   s.t.  -x1 <= 0,  -x2 <= 0,  -x1 - x2 <= 0
fprintf('==================== VI DU 4 ====================\n');
P4.f     = @(x) x(1) + 2*x(2);
P4.gradf = @(x) [1; 2];
P4.g     = @(x) [-x(1); -x(2); -x(1) - x(2)];
P4.jacg  = @(x) [-1 0; 0 -1; -1 -1];
kkt_check(P4, [0; 0]);
fprintf('  Ho nhan tu: lambda = (1-t, 2-t, t), t in [0, 1]. Kiem tra:\n');
for t = [0 0.5 1]
  lam = [1 - t; 2 - t; t];
  fprintf('    t = %.1f: lambda = %s, ||grad_x L|| = %g\n', t, mat2str(lam'), ...
    norm(P4.gradf([0; 0]) + P4.jacg([0; 0])' * lam));
end
fprintf('\n');

%% Vi du 5: khong thoa dieu kien chinh quy - cuc tieu nhung khong phai diem KKT
%   min x1   s.t.  g1 = x2 - x1^3 <= 0,  g2 = -x2 <= 0      (x* = (0, 0))
fprintf('==================== VI DU 5 ====================\n');
P5.f     = @(x) x(1);
P5.gradf = @(x) [1; 0];
P5.g     = @(x) [x(2) - x(1)^3; -x(2)];
P5.jacg  = @(x) [-3*x(1)^2, 1; 0, -1];
kkt_check(P5, [0; 0]);

%% Vi du 6: Slater khong thoa (bai toan loi)
%   min x   s.t.  g = x^2 <= 0      (mien chap nhan chi gom x = 0)
fprintf('==================== VI DU 6 ====================\n');
P6.f     = @(x) x;
P6.gradf = @(x) 1;
P6.g     = @(x) x^2;
P6.jacg  = @(x) 2*x;
kkt_check(P6, 0);

%% Hinh 1: vi du 2 - duong muc, mien chap nhan, gradient tai x* = (1, 1)
if ~exist('figures', 'dir'), mkdir('figures'); end
fig = figure('visible', 'off');
hold on; grid on; box on;
t = linspace(-2, 1, 200);                     % mien: x1^2 <= x2 <= 2 - x1
hF = fill([t, fliplr(t)], [t.^2, 2 - fliplr(t)], [0.85 0.92 1], 'EdgeColor', 'none');
[X1, X2] = meshgrid(linspace(-1.5, 2.5, 300), linspace(-0.5, 3, 300));
[~, hC] = contour(X1, X2, (X1 - 2).^2 + (X2 - 1).^2, [0.25 1 2 4 6.5], ...
  'LineColor', [0.65 0.65 0.65]);
t = linspace(-1.5, 2.5, 200);
hG1 = plot(t, t.^2, 'b-', 'LineWidth', 1.2);
hG2 = plot(t, 2 - t, 'r-', 'LineWidth', 1.2);
xs = [1; 1];
s = 0.35;                                     % ty le mui ten
gf = P2.gradf(xs); J = P2.jacg(xs);
draw_arrow(xs, -s * gf, [0 0 0]);
draw_arrow(xs, s * J(1, :), [0 0 1]);
draw_arrow(xs, s * J(2, :), [1 0 0]);
plot(2, 1, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(1, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(1.45, 1.13, '-\nabla f', 'FontSize', 12);
text(1.55, 0.5, '\nabla g_1', 'Color', 'b', 'FontSize', 12);
text(1.05, 1.5, '\nabla g_2', 'Color', 'r', 'FontSize', 12);
text(2.08, 0.85, '(2, 1)', 'FontSize', 10);
text(0.55, 1.15, 'x^*', 'FontSize', 12);
axis equal; axis([-1.5 2.5 -0.5 3]);
xlabel('x_1'); ylabel('x_2');
legend([hF, hC, hG1, hG2], {'Mien chap nhan', 'Duong muc f', ...
  'g_1 = 0: x_2 = x_1^2', 'g_2 = 0: x_1 + x_2 = 2'}, 'Location', 'northwest');
title('Vi du 2: -\nabla f(x^*) = \lambda_1\nabla g_1 + \lambda_2\nabla g_2,  \lambda = (2/3, 2/3)');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_kkt_geometry.png'));
close(fig);

%% Hinh 2: vi du 5 - mien chap nhan co diem nhon tai (0, 0)
fig = figure('visible', 'off');
hold on; grid on; box on;
t = linspace(0, 1.2, 200);
hF = fill([t, fliplr(t)], [t.^3, zeros(size(t))], [0.85 0.92 1], 'EdgeColor', 'none');
tt = linspace(-0.6, 1.2, 200);
hG1 = plot(tt, tt.^3, 'b-', 'LineWidth', 1.2);
hG2 = plot([-0.6 1.2], [0 0], 'r-', 'LineWidth', 1.2);
draw_arrow([0 0], [0 0.35], [0 0 1], 0.05);
draw_arrow([0 0], [0 -0.35], [1 0 0], 0.05);
draw_arrow([0 0], [-0.35 0], [0 0 0], 0.05);
plot(0, 0, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(0.04, 0.3, '\nabla g_1 = (0, 1)', 'Color', 'b', 'FontSize', 11);
text(0.04, -0.3, '\nabla g_2 = (0, -1)', 'Color', 'r', 'FontSize', 11);
text(-0.58, 0.08, '-\nabla f = (-1, 0)', 'FontSize', 11);
axis([-0.6 1.2 -0.5 0.8]);
xlabel('x_1'); ylabel('x_2');
legend([hF, hG1, hG2], {'Mien chap nhan', 'g_1 = 0: x_2 = x_1^3', 'g_2 = 0: x_2 = 0'}, ...
  'Location', 'northwest');
title('Vi du 5: -\nabla f(x^*) khong nam trong non sinh boi \nabla g_1, \nabla g_2');
print(fig, '-dpng', '-r150', fullfile('figures', 'example5_no_kkt.png'));
close(fig);
fprintf('Da luu hinh figures/example2_kkt_geometry.png, figures/example5_no_kkt.png\n');
