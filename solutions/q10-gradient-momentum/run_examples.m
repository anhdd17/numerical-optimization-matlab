% RUN_EXAMPLES  Cac vi du minh hoa thuat toan gradient va gradient quan tinh (cau 10).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q07-golden-section-search'));   % golden_section.m
addpath(fullfile('..', 'q08-gauss-seidel'));            % bracket_minimum.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: ham bac hai hai bien, kappa = 10 - zic-zac va quan tinh
%   f(x) = 1/2 (x1^2 + 10 x2^2),  L = 10, mu = 1
fprintf('==================== VI DU 1 ====================\n');
A1 = diag([1 10]);  L1 = 10;  mu1 = 1;  kap1 = L1 / mu1;
f1 = @(x) 0.5 * x' * A1 * x;
g1 = @(x) A1 * x;
x01 = [10; 1];
aGD = 2 / (L1 + mu1);
aHB = 4 / (sqrt(L1) + sqrt(mu1))^2;
bHB = ((sqrt(kap1) - 1) / (sqrt(kap1) + 1))^2;
fprintf('Buoc toi uu cho gradient: alpha = 2/(L+mu) = %.4f\n', aGD);
fprintf('Tham so toi uu cho quan tinh: alpha = %.4f, beta = %.4f\n', aHB, bHB);
[~, ~, o1a] = gradient_descent(f1, g1, x01, struct('step', 'fixed', 'alpha', aGD, 'display', true, 'maxIter', 3));
[~, ~, o1a] = gradient_descent(f1, g1, x01, struct('step', 'fixed', 'alpha', aGD));
[~, ~, o1b] = gradient_descent(f1, g1, x01, struct('step', 'exact'));
[~, ~, o1c] = momentum_gradient(f1, g1, x01, struct('alpha', aHB, 'beta', bHB));
[~, ~, o1d] = gradient_descent(f1, g1, x01, struct('step', 'fixed', 'alpha', 0.21, 'maxIter', 30));
[~, ~, o1e] = gradient_descent(f1, g1, [1; 1], struct('step', 'exact'));
[~, ~, o1f] = gradient_descent(f1, g1, [1; 1], struct('step', 'fixed', 'alpha', aGD));
fprintf('\nSo buoc de ||grad f|| <= 1e-6:\n');
fprintf('  x0 = (10, 1): gradient buoc 2/(L+mu): %d; tim kiem chinh xac: %d; quan tinh: %d\n', ...
  o1a.iterations, o1b.iterations, o1c.iterations);
fprintf('  x0 = (1, 1):  gradient buoc 2/(L+mu): %d; tim kiem chinh xac: %d\n', ...
  o1f.iterations, o1e.iterations);
fprintf('  Cac buoc tim kiem chinh xac tu (10, 1): %s (= 2/(L+mu))\n', mat2str(o1b.alphas(1:4), 6));
fprintf('  Buoc qua lon alpha = 0.21 > 2/L = 0.2: ||x_30|| = %.3e (phan ky)\n\n', norm(o1d.X(:, end)));

%% Vi du 2: so buoc theo so dieu kien kappa
fprintf('==================== VI DU 2 ====================\n');
fprintf('%8s | %10s %10s | %10s %10s | %8s\n', 'kappa', 'gradient', 'ly thuyet', 'quan tinh', 'ly thuyet', 'ti so');
kaps = [10 100 1000];
itG = zeros(size(kaps));  itH = itG;
for j = 1:numel(kaps)
  kap = kaps(j);
  A = diag([1 kap]);  f = @(x) 0.5 * x' * A * x;  g = @(x) A * x;
  x0 = [kap; 1];  g0 = norm(g(x0));
  [~, ~, oG] = gradient_descent(f, g, x0, struct('alpha', 2 / (kap + 1), 'maxIter', 1e6));
  a = 4 / (sqrt(kap) + 1)^2;  b = ((sqrt(kap) - 1) / (sqrt(kap) + 1))^2;
  [~, ~, oH] = momentum_gradient(f, g, x0, struct('alpha', a, 'beta', b, 'maxIter', 1e6));
  qG = (kap - 1) / (kap + 1);  qH = (sqrt(kap) - 1) / (sqrt(kap) + 1);
  itG(j) = oG.iterations;  itH(j) = oH.iterations;
  fprintf('%8d | %10d %10.0f | %10d %10.0f | %8.1f\n', kap, oG.iterations, ...
    ceil(log(1e-6 / g0) / log(qG)), oH.iterations, ceil(log(1e-6 / g0) / log(qH)), ...
    oG.iterations / oH.iterations);
end
% n = 10: ham bac hai cua cau 8, vi du 3
n = 10;
A3 = 2 * eye(n) - diag(ones(n - 1, 1), 1) - diag(ones(n - 1, 1), -1);
b3 = ones(n, 1);
f3 = @(x) 0.5 * x' * A3 * x - b3' * x;  g3 = @(x) A3 * x - b3;
ev = eig(A3);  L3 = max(ev);  mu3 = min(ev);  kap3 = L3 / mu3;
[x3g, ~, o3g] = gradient_descent(f3, g3, zeros(n, 1), struct('alpha', 2 / (L3 + mu3), 'tol', 1e-8, 'maxIter', 1e6));
[x3h, ~, o3h] = momentum_gradient(f3, g3, zeros(n, 1), struct('alpha', 4 / (sqrt(L3) + sqrt(mu3))^2, ...
  'beta', ((sqrt(kap3) - 1) / (sqrt(kap3) + 1))^2, 'tol', 1e-8, 'maxIter', 1e6));
fprintf('n = 10, A = tridiag(-1, 2, -1): kappa = %.2f\n', kap3);
fprintf('  gradient: %d buoc, ||x - x*|| = %.1e;  quan tinh: %d buoc, ||x - x*|| = %.1e\n\n', ...
  o3g.iterations, norm(x3g - A3 \ b3), o3h.iterations, norm(x3h - A3 \ b3));

%% Vi du 3: ham Bazaraa - gradient voi tim kiem chinh xac theo tia
%   F(x) = (x1 - 2)^4 + (x1 - 2x2)^2,  x0 = (0, 3)
fprintf('==================== VI DU 3 ====================\n');
F = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
G = @(x) [4*(x(1) - 2)^3 + 2*(x(1) - 2*x(2)); -4*(x(1) - 2*x(2))];
[~, ~, o3s] = gradient_descent(F, G, [0; 3], struct('step', 'exact', 'display', true, 'maxIter', 6));
for t = [1e-2 1e-3]
  [x, fx, o] = gradient_descent(F, G, [0; 3], struct('step', 'exact', 'tol', t, 'maxIter', 1e5));
  fprintf('tol = %g: %d buoc, x = %s, F = %.2e, ||x - x*|| = %.2e, %d lan tinh F\n', ...
    t, o.iterations, mat2str(x', 6), fx, norm(x - [2; 1]), o.nfev);
end
fprintf('\n');

%% Vi du 4: ham Rosenbrock
%   R(x) = 100 (x2 - x1^2)^2 + (1 - x1)^2,  x0 = (-1.2, 1)
fprintf('==================== VI DU 4 ====================\n');
R  = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
RG = @(x) [-400 * x(1) * (x(2) - x(1)^2) - 2 * (1 - x(1)); 200 * (x(2) - x(1)^2)];
x04 = [-1.2; 1];
opt4 = struct('tol', 1e-4, 'maxIter', 1e5);
[x4a, ~, o4a] = gradient_descent(R, RG, x04, setfield(setfield(opt4, 'step', 'fixed'), 'alpha', 1e-3));
[x4b, ~, o4b] = gradient_descent(R, RG, x04, setfield(opt4, 'step', 'armijo'));
[x4c, ~, o4c] = momentum_gradient(R, RG, x04, setfield(setfield(opt4, 'alpha', 1e-3), 'beta', 0.9));
[x4d, ~, o4d] = gradient_descent(R, RG, x04, setfield(setfield(opt4, 'step', 'fixed'), 'alpha', 2.5e-3));
fprintf('%-34s %8s %12s %10s %10s\n', 'Phuong phap (tol ||grad|| = 1e-4)', 'buoc', '||x - x*||', 'tinh f', 'tinh grad');
fprintf('%-34s %8d %12.2e %10d %10d\n', 'Gradient, alpha = 1e-3', o4a.iterations, norm(x4a - 1), o4a.nfev, o4a.ngev);
fprintf('%-34s %8d %12.2e %10d %10d\n', 'Gradient, Armijo', o4b.iterations, norm(x4b - 1), o4b.nfev, o4b.ngev);
fprintf('%-34s %8d %12.2e %10s %10d\n', 'Quan tinh, alpha = 1e-3, beta = 0.9', o4c.iterations, norm(x4c - 1), '-', o4c.ngev);
H4 = [1200 * 1 - 400 * 1 + 2, -400; -400, 200];   % Hessian cua R tai x* = (1, 1)
fprintf('Hessian tai x*: tri rieng %s, 2/L = %.6f, kappa = %.0f\n', mat2str(eig(H4)', 6), ...
  2 / max(eig(H4)), max(eig(H4)) / min(eig(H4)));
fprintf('Gradient, alpha = 2.5e-3 > 2/L: exitflag = %d sau %d buoc; 2 diem cuoi %s va %s, R = %s\n\n', ...
  o4d.exitflag, o4d.iterations, mat2str(o4d.X(:, end-1)', 4), mat2str(o4d.X(:, end)', 4), ...
  mat2str(o4d.F(end-1:end), 4));

%% Hinh 1: vi du 1 - duong di tren duong muc
fig = figure('visible', 'off', 'Position', [100 100 1150 420]);
[X1, X2] = meshgrid(linspace(-11, 11, 300), linspace(-1.6, 1.6, 300));
Z = 0.5 * (X1.^2 + 10 * X2.^2);
paths = {o1a.X, o1c.X, o1d.X};
ttl = {sprintf('Gradient, \\alpha = 2/(L+\\mu): %d buoc', o1a.iterations), ...
       sprintf('Quan tinh (toi uu): %d buoc', o1c.iterations), ...
       'Gradient, \alpha = 0.21 > 2/L: phan ky'};
col = {'b', [0 0.5 0], 'r'};
for p = 1:3
  subplot(1, 3, p); hold on; grid on; box on;
  contour(X1, X2, Z, [0.5 2 5 10 20 35 55], 'LineColor', [0.6 0.6 0.6]);
  P = paths{p};
  if p == 3, P = P(:, 1:13); end
  plot(P(1, :), P(2, :), '.-', 'Color', col{p}, 'LineWidth', 1.2, 'MarkerSize', 10);
  plot(x01(1), x01(2), 'ko', 'MarkerFaceColor', 'k');
  plot(0, 0, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.8 0]);
  axis([-11 11 -1.6 1.6]);
  xlabel('x_1'); ylabel('x_2');
  title(ttl{p});
end
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_paths.png'));
close(fig);

%% Hinh 2: vi du 2 va 4 - toc do hoi tu
fig = figure('visible', 'off', 'Position', [100 100 1100 440]);
subplot(1, 2, 1);
kap = 100;  A = diag([1 kap]);  f = @(x) 0.5 * x' * A * x;  g = @(x) A * x;
[~, ~, oG] = gradient_descent(f, g, [kap; 1], struct('alpha', 2 / (kap + 1), 'maxIter', 1e6));
[~, ~, oH] = momentum_gradient(f, g, [kap; 1], struct('alpha', 4 / (sqrt(kap) + 1)^2, ...
  'beta', ((sqrt(kap) - 1) / (sqrt(kap) + 1))^2, 'maxIter', 1e6));
semilogy(0:oG.iterations, oG.G, 'b-', 'LineWidth', 1.4); hold on; grid on;
semilogy(0:oH.iterations, oH.G, '-', 'Color', [0 0.5 0], 'LineWidth', 1.4);
xlabel('buoc lap k'); ylabel('||\nabla f(x_k)||');
legend('Gradient, \alpha = 2/(L+\mu)', 'Quan tinh (toi uu)', 'Location', 'northeast');
title('Vi du 2: ham bac hai, \kappa = 100');
subplot(1, 2, 2);
semilogy(0:o4a.iterations, o4a.G, 'b-', 'LineWidth', 1.2); hold on; grid on;
semilogy(0:o4b.iterations, o4b.G, 'm-', 'LineWidth', 1.2);
semilogy(0:o4c.iterations, o4c.G, '-', 'Color', [0 0.5 0], 'LineWidth', 1.2);
xlabel('buoc lap k'); ylabel('||\nabla R(x_k)||');
legend('Gradient, \alpha = 10^{-3}', 'Gradient, Armijo', 'Quan tinh, \beta = 0.9', 'Location', 'northeast');
title('Vi du 4: ham Rosenbrock');
print(fig, '-dpng', '-r150', fullfile('figures', 'examples2_4_convergence.png'));
close(fig);

%% Hinh 3: vi du 3 - zic-zac cua gradient voi tim kiem chinh xac
fig = figure('visible', 'off');
hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-0.5, 3.5, 300), linspace(0, 3.5, 300));
Z = (X1 - 2).^4 + (X1 - 2 * X2).^2;
contour(X1, X2, Z, [0.01 0.05 0.1 0.37 1 3 8 16 30 52], 'LineColor', [0.6 0.6 0.6]);
[~, ~, o3p] = gradient_descent(F, G, [0; 3], struct('step', 'exact', 'maxIter', 25));
plot(o3p.X(1, :), o3p.X(2, :), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 10);
plot(0, 3, 'ko', 'MarkerFaceColor', 'k');
plot(2, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(0.1, 3.2, 'x_0 = (0, 3)', 'FontSize', 10);
text(1.2, 0.55, 'x^* = (2, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
axis equal; axis([-0.5 3.5 0 3.5]);
xlabel('x_1'); ylabel('x_2');
title('Vi du 3: gradient + tim kiem chinh xac, 25 buoc dau');
print(fig, '-dpng', '-r150', fullfile('figures', 'example3_zigzag.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
