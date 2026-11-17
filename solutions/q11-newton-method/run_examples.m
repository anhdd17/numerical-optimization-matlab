% RUN_EXAMPLES  Cac vi du minh hoa thuat toan Newton (cau 11).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q07-golden-section-search'));   % golden_section.m
addpath(fullfile('..', 'q08-gauss-seidel'));            % bracket_minimum.m
addpath(fullfile('..', 'q10-gradient-momentum'));       % gradient_descent.m, momentum_gradient.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: ham bac hai n = 10 va ham khong bac hai n = 10
%   f(x) = 1/2 x'Ax - b'x,  A = tridiag(-1, 2, -1),  b = (1, ..., 1)
fprintf('==================== VI DU 1 ====================\n');
n = 10;
A = 2 * eye(n) - diag(ones(n - 1, 1), 1) - diag(ones(n - 1, 1), -1);
b = ones(n, 1);
fq = @(x) 0.5 * x' * A * x - b' * x;
gq = @(x) A * x - b;
hq = @(x) A;
xs = A \ b;
[x1, ~, o1] = newton_method(fq, gq, hq, zeros(n, 1), struct('tol', 1e-8));
fprintf('Bac hai n = 10 (kappa = %.2f): Newton %d buoc, ||x - x*|| = %.1e, ||grad|| sau buoc 1 = %.1e\n', ...
  cond(A), o1.iterations, norm(x1 - xs), o1.G(2));
for kap = [1e2 1e4 1e6]
  D = diag(logspace(0, log10(kap), n));
  [xk, ~, ok] = newton_method(@(x) 0.5 * x' * D * x - b' * x, @(x) D * x - b, @(x) D, 5 * ones(n, 1), struct('tol', 1e-8));
  fprintf('  bac hai duong cheo kappa = %g: Newton %d buoc, ||x - x*|| = %.1e\n', kap, ok.iterations, norm(xk - D \ b));
end
% Ham khong bac hai: them sum x_i^4 / 4 (loi ngat)
f1 = @(x) 0.5 * x' * A * x - b' * x + sum(x.^4) / 4;
g1 = @(x) A * x - b + x.^3;
h1 = @(x) A + diag(3 * x.^2);
x01 = 2 * ones(n, 1);
fprintf('\nKhong bac hai n = 10: f(x) = 1/2 x''Ax - b''x + sum x_i^4/4, x0 = (2, ..., 2)\n');
[x1n, f1n, o1n] = newton_method(f1, g1, h1, x01, struct('tol', 1e-12, 'display', true));
x1ref = x1n;
[~, ~, o1g] = gradient_descent(f1, g1, x01, struct('step', 'armijo', 'tol', 1e-8, 'maxIter', 1e5));
[~, ~, o1e] = gradient_descent(f1, g1, x01, struct('step', 'exact', 'tol', 1e-8, 'maxIter', 1e5));
fprintf('x* = %s\n', mat2str(x1n', 8));
fprintf('Sai so ||x_k - x*||: %s\n', mat2str(sqrt(sum((o1n.X - x1ref).^2, 1)), 3));
fprintf('Gradient (Armijo) can %d buoc, gradient (chinh xac) %d buoc de ||grad|| <= 1e-8\n\n', ...
  o1g.iterations, o1e.iterations);

%% Vi du 2: ham Rosenbrock - Newton thuan tuy va Newton giam buoc
%   R(x) = 100 (x2 - x1^2)^2 + (1 - x1)^2,  x0 = (-1.2, 1)
fprintf('==================== VI DU 2 ====================\n');
R  = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
RG = @(x) [-400 * x(1) * (x(2) - x(1)^2) - 2 * (1 - x(1)); 200 * (x(2) - x(1)^2)];
RH = @(x) [1200 * x(1)^2 - 400 * x(2) + 2, -400 * x(1); -400 * x(1), 200];
x02 = [-1.2; 1];
fprintf('-- Newton thuan tuy (alpha = 1) --\n');
[x2p, ~, o2p] = newton_method(R, RG, RH, x02, struct('variant', 'pure', 'display', true));
e2p = sqrt(sum((o2p.X - 1).^2, 1));
fprintf('e_k = ||x_k - x*||: %s\n', mat2str(e2p, 3));
fprintf('e_{k+1} / e_k^2:    %s\n', mat2str(e2p(2:end) ./ e2p(1:end-1).^2, 3));
fprintf('-- Newton giam buoc (Armijo) --\n');
[x2d, ~, o2d] = newton_method(R, RG, RH, x02, struct('variant', 'damped', 'display', true));
e2d = sqrt(sum((o2d.X - 1).^2, 1));
fprintf('Cac buoc alpha_k: %s\n', mat2str(o2d.alphas, 3));
fprintf('e_k cua 5 buoc cuoi: %s\n', mat2str(e2d(end-4:end), 3));
fprintf('R(x_k) co giam don dieu? thuan tuy: %d, giam buoc: %d\n', all(diff(o2p.F) < 0), all(diff(o2d.F) < 0));
[~, ~, o2g] = gradient_descent(R, RG, x02, struct('step', 'armijo', 'tol', 1e-10, 'maxIter', 1e6));
fprintf('Doi chieu: gradient Armijo can %d buoc de ||grad|| <= 1e-10\n\n', o2g.iterations);

%% Vi du 3: f(x) = sqrt(1 + x^2) - Newton thuan tuy phan ky
%   f'(x) = x / sqrt(1 + x^2), f''(x) = (1 + x^2)^(-3/2)  =>  x_{k+1} = -x_k^3
fprintf('==================== VI DU 3 ====================\n');
f3 = @(x) sqrt(1 + x^2);
g3 = @(x) x / sqrt(1 + x^2);
h3 = @(x) (1 + x^2)^(-1.5);
for x0 = [0.5 0.9 1 1.1]
  [x3, ~, o3] = newton_method(f3, g3, h3, x0, struct('maxIter', 12));
  fprintf('x0 = %-4g: %2d buoc, exitflag = %2d, x_1..x_4 = %s, x cuoi = %.4g\n', x0, o3.iterations, ...
    o3.exitflag, mat2str(o3.X(2:min(5, end)), 4), x3);
end
[x3d, ~, o3d] = newton_method(f3, g3, h3, 1.1, struct('variant', 'damped', 'display', true));
[x3e, ~, o3e] = newton_method(f3, g3, h3, 10, struct('variant', 'damped'));
fprintf('Giam buoc tu x0 = 10: %d buoc, alpha = %s, x = %.2e\n\n', o3e.iterations, mat2str(o3e.alphas, 3), x3e);

%% Vi du 4: ham khong loi - Hessian khong xac dinh duong
%   f(x) = x1^4/4 - x1^2/2 + x2^2/2: cuc tieu (+-1, 0), diem yen ngua (0, 0)
fprintf('==================== VI DU 4 ====================\n');
f4 = @(x) x(1)^4 / 4 - x(1)^2 / 2 + x(2)^2 / 2;
g4 = @(x) [x(1)^3 - x(1); x(2)];
h4 = @(x) [3 * x(1)^2 - 1, 0; 0, 1];
x04 = [0.3; 1];
fprintf('Hessian tai x0: %s\n', mat2str(h4(x04), 4));
fprintf('-- Newton thuan tuy --\n');
[x4p, f4p, o4p] = newton_method(f4, g4, h4, x04, struct('display', true));
fprintf('-- Newton giam buoc --\n');
[x4d, f4d, o4d] = newton_method(f4, g4, h4, x04, struct('variant', 'damped', 'display', true));
d4 = -h4(x4d) \ g4(x4d);
fprintf('Tai x_%d: grad = %s, d = %s, grad''d = %.3e > 0\n', o4d.iterations, mat2str(g4(x4d)', 4), ...
  mat2str(d4', 4), g4(x4d)' * d4);
fprintf('-- Newton cai bien (H + tau I) --\n');
[x4m, f4m, o4m] = newton_method(f4, g4, h4, x04, struct('variant', 'modified', 'display', true));
fprintf('Ket qua: thuan tuy -> %s (f = %.4g, tri rieng H = %s); giam buoc: exitflag = %d; cai bien -> %s (f = %.4g)\n\n', ...
  mat2str(x4p' + 0, 4), f4p, mat2str(eig(h4(x4p))', 3), o4d.exitflag, mat2str(x4m' + 0, 6), f4m);

%% Vi du 5: ham Bazaraa - Hessian suy bien tai x*, hoi tu tuyen tinh
%   F(x) = (x1 - 2)^4 + (x1 - 2x2)^2,  x0 = (0, 3),  x* = (2, 1)
fprintf('==================== VI DU 5 ====================\n');
F = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
G = @(x) [4*(x(1) - 2)^3 + 2*(x(1) - 2*x(2)); -4*(x(1) - 2*x(2))];
H = @(x) [12*(x(1) - 2)^2 + 2, -4; -4, 8];
fprintf('Hessian tai x* = (2, 1): %s, det = %g\n', mat2str(H([2; 1])), det(H([2; 1])));
[x5, F5, o5] = newton_method(F, G, H, [0; 3], struct('tol', 1e-8, 'display', true, 'maxIter', 6));
[x5, F5, o5] = newton_method(F, G, H, [0; 3], struct('tol', 1e-8));
e5 = sqrt(sum((o5.X - [2; 1]).^2, 1));
fprintf('tol = 1e-8: %d buoc, x = %s, F = %.2e, ||x - x*|| = %.2e\n', o5.iterations, mat2str(x5', 6), F5, e5(end));
fprintf('Ti so e_{k+1}/e_k: %s\n', mat2str(e5(2:end) ./ e5(1:end-1), 4));
fprintf('x1_k - 2: %s\n', mat2str(o5.X(1, 1:6) - 2, 5));
[x5g, ~, o5g] = gradient_descent(F, G, [0; 3], struct('step', 'exact', 'tol', 1e-3, 'maxIter', 1e5));
fprintf('Doi chieu (cau 10): gradient tim kiem chinh xac, tol 1e-3: %d buoc, ||x - x*|| = %.3f\n\n', ...
  o5g.iterations, norm(x5g - [2; 1]));

%% Hinh 1: vi du 2 - duong di tren duong muc Rosenbrock va toc do hoi tu
fig = figure('visible', 'off', 'Position', [100 100 1150 460]);
subplot(1, 2, 1); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-1.6, 1.6, 400), linspace(-0.6, 2.6, 400));
Z = 100 * (X2 - X1.^2).^2 + (1 - X1).^2;
contour(X1, X2, Z, [0.5 2 5 10 25 50 100 200 400], 'LineColor', [0.7 0.7 0.7]);
hp = plot(o2p.X(1, :), o2p.X(2, :), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12);
hd = plot(o2d.X(1, :), o2d.X(2, :), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 12);
plot(x02(1), x02(2), 'ko', 'MarkerFaceColor', 'k');
plot(1, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(-1.5, 1.2, 'x_0', 'FontSize', 10);
text(1.08, 0.9, 'x^*', 'FontSize', 10);
axis([-1.6 1.6 -0.6 2.6]);
xlabel('x_1'); ylabel('x_2');
legend([hp, hd], {sprintf('Newton thuan tuy (%d buoc)', o2p.iterations), ...
  sprintf('Newton giam buoc (%d buoc)', o2d.iterations)}, 'Location', 'northeast');
text(0.8, -0.45, sprintf('x_2 = (%.2f, %.2f) \\downarrow', o2p.X(1, 3), o2p.X(2, 3)), 'FontSize', 9, ...
  'Color', 'r', 'BackgroundColor', 'w', 'HorizontalAlignment', 'center');
title('Vi du 2: ham Rosenbrock');
subplot(1, 2, 2);
semilogy(0:o2p.iterations, max(e2p, 1e-17), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12); hold on; grid on;
semilogy(0:o2d.iterations, max(e2d, 1e-17), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 12);
semilogy(0:o2g.iterations, sqrt(sum(([o2g.X] - 1).^2, 1)), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
xlim([0 40]); ylim([1e-17 1e2]);
xlabel('buoc lap k'); ylabel('||x_k - x^*||');
legend('Newton thuan tuy', 'Newton giam buoc', sprintf('Gradient Armijo (%d buoc)', o2g.iterations), ...
  'Location', 'southwest');
title('Sai so (thang log)');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_rosenbrock.png'));
close(fig);

%% Hinh 2: vi du 3 - Newton cho f(x) = sqrt(1 + x^2) la Newton giai f'(x) = 0
fig = figure('visible', 'off', 'Position', [100 100 1150 440]);
for p = 1:2
  subplot(1, 2, p); hold on; grid on; box on;
  t = linspace(-4, 4, 400);
  plot(t, t ./ sqrt(1 + t.^2), 'k-', 'LineWidth', 1.4);
  plot([-4 4], [0 0], 'k-', 'LineWidth', 0.5);
  x0 = 0.9;  if p == 2, x0 = 1.1; end
  xk = x0;
  for k = 1:3
    yk = g3(xk);  xn = -xk^3;
    plot([xk xk], [0 yk], 'k:', 'LineWidth', 1);
    plot([xk xn], [yk 0], 'r-', 'LineWidth', 1.3);
    plot(xk, yk, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
    if k == 1, dx = 0.25; else, dx = 0.22 * sign(xk); end   % lech khoi tiep tuyen den tu x_{k-1}
    text(xk + dx, -0.12 * sign(yk), sprintf('x_%d', k - 1), 'HorizontalAlignment', 'center', 'FontSize', 10);
    xk = xn;
    if abs(xk) > 4, break; end
  end
  axis([-4 4 -1.1 1.1]);
  xlabel('x'); ylabel('f''(x) = x / (1 + x^2)^{1/2}');
  if p == 1
    title('x_0 = 0.9: x_{k+1} = -x_k^3 \rightarrow 0');
  else
    title('x_0 = 1.1: |x_k| tang, phan ky');
  end
end
print(fig, '-dpng', '-r150', fullfile('figures', 'example3_divergence.png'));
close(fig);

%% Hinh 3: vi du 4 - ham khong loi
fig = figure('visible', 'off', 'Position', [100 100 760 520]);
hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-1.6, 1.6, 300), linspace(-1.2, 1.4, 300));
Z = X1.^4 / 4 - X1.^2 / 2 + X2.^2 / 2;
contour(X1, X2, Z, [-0.24 -0.2 -0.12 -0.05 0 0.05 0.15 0.3 0.5 0.8], 'LineColor', [0.7 0.7 0.7]);
hp = plot(o4p.X(1, :), o4p.X(2, :), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 13);
hd = plot(o4d.X(1, :), o4d.X(2, :), 'b--', 'LineWidth', 1.2);
plot(o4d.X(1, end), o4d.X(2, end), 'bs', 'MarkerSize', 9, 'MarkerFaceColor', 'b');
hm = plot(o4m.X(1, :), o4m.X(2, :), '.-', 'Color', [0 0.5 0], 'LineWidth', 1.2, 'MarkerSize', 13);
plot(x04(1), x04(2), 'ko', 'MarkerFaceColor', 'k');
plot([-1 1], [0 0], 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
plot(0, 0, 'kx', 'MarkerSize', 12, 'LineWidth', 2);
text(0.35, 1.05, 'x_0', 'FontSize', 10);
text(0.05, -0.12, 'yen ngua (0, 0)', 'FontSize', 10, 'BackgroundColor', 'w');
text(0.82, -0.17, 'cuc tieu (1, 0)', 'FontSize', 10, 'BackgroundColor', 'w');
text(-1.25, -0.17, '(-1, 0)', 'FontSize', 10, 'BackgroundColor', 'w');
axis([-1.6 1.6 -1.2 1.4]);
xlabel('x_1'); ylabel('x_2');
legend([hp, hd, hm], {'Newton thuan tuy \rightarrow yen ngua', ...
  'Newton giam buoc: dung (d khong la huong giam)', 'Newton cai bien \rightarrow cuc tieu'}, ...
  'Location', 'south');
title('Vi du 4: f = x_1^4/4 - x_1^2/2 + x_2^2/2');
print(fig, '-dpng', '-r150', fullfile('figures', 'example4_nonconvex.png'));
close(fig);

%% Hinh 4: vi du 1 va 5 - so sanh toc do hoi tu
fig = figure('visible', 'off', 'Position', [100 100 1100 440]);
subplot(1, 2, 1);
e1n = sqrt(sum((o1n.X - x1ref).^2, 1));
semilogy(0:o1n.iterations, max(e1n, 1e-17), 'r.-', 'LineWidth', 1.3, 'MarkerSize', 12); hold on; grid on;
semilogy(0:o1g.iterations, sqrt(sum((o1g.X - x1ref).^2, 1)), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 8);
xlim([0 50]); ylim([1e-17 1e2]);
xlabel('buoc lap k'); ylabel('||x_k - x^*||');
legend(sprintf('Newton (%d buoc)', o1n.iterations), sprintf('Gradient Armijo (%d buoc)', o1g.iterations), ...
  'Location', 'northeast');
title('Vi du 1: n = 10, khong bac hai');
subplot(1, 2, 2);
semilogy(0:o5.iterations, e5, 'r.-', 'LineWidth', 1.3, 'MarkerSize', 12); hold on; grid on;
semilogy(0:o5.iterations, e5(2) * (2/3).^((0:o5.iterations) - 1), 'k--', 'LineWidth', 1);
semilogy(0:o5g.iterations, sqrt(sum((o5g.X - [2; 1]).^2, 1)), 'b-', 'LineWidth', 1.2);
xlim([0 60]); ylim([5e-4 50]);
xlabel('buoc lap k'); ylabel('||x_k - x^*||');
legend('Newton', 'C (2/3)^k', 'Gradient, tim kiem chinh xac', 'Location', 'northeast');
title('Vi du 5: ham Bazaraa (Hessian suy bien tai x^*)');
print(fig, '-dpng', '-r150', fullfile('figures', 'examples1_5_convergence.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
