% RUN_EXAMPLES  Cac vi du minh hoa giai thuat tim kiem truc tiep Hooke-Jeeves (cau 9).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q08-gauss-seidel'));            % gauss_seidel_min.m (so sanh)
addpath(fullfile('..', 'q07-golden-section-search'));   % golden_section.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: ham Bazaraa (vi du 2 cua cau 8) - buoc theo mau vuot thung lung
%   F(x) = (x1 - 2)^4 + (x1 - 2x2)^2,  x0 = (0, 3),  x* = (2, 1)
fprintf('==================== VI DU 1 ====================\n');
F1 = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
fprintf('--- Hooke-Jeeves, Delta0 = 0.5, beta = 0.5 ---\n');
[x1, f1, o1] = hooke_jeeves(F1, [0; 3], struct('delta', 0.5, 'display', true));
fprintf('Ket qua: x = %s, F = %.2e, %d lan tinh F\n\n', mat2str(x1', 6), f1, o1.nfev);
fprintf('%-28s %-22s %12s %12s %10s\n', 'Phuong phap', 'x', 'F(x)', '||x - x*||', 'tinh F');
for d = [0.5 0.3 0.7]
  [x, fx, o] = hooke_jeeves(F1, [0; 3], struct('delta', d));
  fprintf('%-28s %-22s %12.2e %12.2e %10d\n', sprintf('Hooke-Jeeves, Delta0 = %.1f', d), ...
    mat2str(x', 6), fx, norm(x - [2; 1]), o.nfev);
end
[xg, fg, og] = gauss_seidel_min(F1, [0; 3], struct('maxIter', 500));
fprintf('%-28s %-22s %12.2e %12.2e %10d\n\n', 'Gauss-Seidel (cau 8), 500 vong', ...
  mat2str(xg', 6), fg, norm(xg - [2; 1]), og.nfev);

%% Vi du 2: ham Rosenbrock - so sanh voi Nelder-Mead va Gauss-Seidel
%   R(x) = 100 (x2 - x1^2)^2 + (1 - x1)^2,  x0 = (-1.2, 1),  x* = (1, 1)
fprintf('==================== VI DU 2 ====================\n');
R = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
x0 = [-1.2; 1];
[x2, f2, o2] = hooke_jeeves(R, x0, struct('delta', 0.5, 'tol', 1e-6));
[xn, fn, ~, on] = fminsearch(R, x0, optimset('TolX', 1e-8, 'TolFun', 1e-12, ...
  'MaxFunEvals', 1e5, 'MaxIter', 1e5));
[xg2, fg2, og2] = gauss_seidel_min(R, x0, struct('maxIter', 2000, 'ftol', 0));
fprintf('%-30s %-26s %12s %12s %10s\n', 'Phuong phap', 'x', 'R(x)', '||x - x*||', 'tinh R');
fprintf('%-30s %-26s %12.2e %12.2e %10d\n', 'Hooke-Jeeves', mat2str(x2', 8), f2, norm(x2 - [1; 1]), o2.nfev);
fprintf('%-30s %-26s %12.2e %12.2e %10d\n', 'Nelder-Mead (fminsearch)', mat2str(xn', 8), fn, norm(xn - [1; 1]), on.funcCount);
fprintf('%-30s %-26s %12.2e %12.2e %10d\n\n', 'Gauss-Seidel, 2000 vong', mat2str(xg2', 8), fg2, norm(xg2 - [1; 1]), og2.nfev);

%% Vi du 3: n = 10 bien
fprintf('==================== VI DU 3 ====================\n');
n = 10;
% (a) ham bac hai cua cau 8, vi du 3: f = 1/2 x'Ax - b'x, A = tridiag(-1, 2, -1), b = 1
A3 = 2 * eye(n) - diag(ones(n - 1, 1), 1) - diag(ones(n - 1, 1), -1);
b3 = ones(n, 1);
f3 = @(x) 0.5 * x' * A3 * x - b3' * x;
xs3 = A3 \ b3;
[x3, ~, o3] = hooke_jeeves(f3, zeros(n, 1), struct('delta', 0.7, 'tol', 1e-8));   % Delta0 = 1 se roi dung nghiem nguyen tren luoi
fprintf('(a) bac hai: ||x - x*|| = %.2e, %d lan tinh f, %d diem co so, %d lan giam buoc\n', ...
  norm(x3 - xs3), o3.nfev, o3.iterations, o3.reductions);
% (b) ham khong tron, tach bien: f = sum |x_i - i| + sum (x_i - i)^2
c3 = (1:n)';
g3 = @(x) sum(abs(x - c3)) + sum((x - c3).^2);
[x3b, fx3b, o3b] = hooke_jeeves(g3, zeros(n, 1), struct('delta', 0.7, 'tol', 1e-8));
fprintf('(b) khong tron, tach bien: ||x - x*|| = %.2e, f = %.2e, %d lan tinh f\n\n', ...
  norm(x3b - c3), fx3b, o3b.nfev);

%% Vi du 4: ham khong tron, khong tach bien (vi du 4 cua cau 8)
%   f(x) = 2|x1 - x2| + (x1 + x2 - 3)^2,  x* = (1.5, 1.5), f* = 0
fprintf('==================== VI DU 4 ====================\n');
f4 = @(x) 2 * abs(x(1) - x(2)) + (x(1) + x(2) - 3)^2;
cases = {[0; 0], 0.5; [0; 0], 0.3; [0; 0], 0.7; [0; 0], 1.0; [0; 0.3], 0.3; [0.1; 0.7], 0.3};
stuck = zeros(2, size(cases, 1));
fprintf('%-12s %6s | %-22s %10s\n', 'x0', 'Delta0', 'Hooke-Jeeves dung tai', 'f');
for c = 1:size(cases, 1)
  [x, fx] = hooke_jeeves(f4, cases{c, 1}, struct('delta', cases{c, 2}));
  stuck(:, c) = x;
  fprintf('%-12s %6.1f | %-22s %10.4g\n', mat2str(cases{c, 1}'), cases{c, 2}, mat2str(x', 6), fx);
end
[xn4, fn4] = fminsearch(f4, [0; 0]);
fprintf('Nelder-Mead (fminsearch) tu (0, 0): x = %s, f = %.2e\n\n', mat2str(xn4', 6), fn4);

%% Hinh 1: vi du 1 - Hooke-Jeeves va Gauss-Seidel tren duong muc
fig = figure('visible', 'off', 'Position', [100 100 1100 470]);
[X1, X2] = meshgrid(linspace(-0.5, 3.5, 300), linspace(0, 3.5, 300));
Z = (X1 - 2).^4 + (X1 - 2 * X2).^2;
for p = 1:2
  subplot(1, 2, p); hold on; grid on; box on;
  contour(X1, X2, Z, [0.1 0.37 1 3 8 16 30 52], 'LineColor', [0.6 0.6 0.6]);
  if p == 1
    [~, ~, oa] = hooke_jeeves(F1, [0; 3], struct('delta', 0.3));
    B = oa.base;
    plot(B(1, :), B(2, :), 'b.-', 'LineWidth', 1.4, 'MarkerSize', 12);
    ttl = sprintf('Hooke-Jeeves (\\Delta_0 = 0.3): %d lan tinh F', oa.nfev);
  else
    P = og.path(:, 1:2 * 20 + 1);
    plot(P(1, :), P(2, :), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 10);
    ttl = 'Gauss-Seidel (cau 8): 20 vong dau';
  end
  plot(0, 3, 'ko', 'MarkerFaceColor', 'k');
  plot(2, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
  text(0.1, 3.2, 'x_0 = (0, 3)', 'FontSize', 10);
  text(1.25, 0.55, 'x^* = (2, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
  axis equal; axis([-0.5 3.5 0 3.5]);
  xlabel('x_1'); ylabel('x_2');
  title(ttl);
end
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_pattern_vs_coordinate.png'));
close(fig);

%% Hinh 2: vi du 2 - duong di cua Hooke-Jeeves tren ham Rosenbrock
fig = figure('visible', 'off', 'Position', [100 100 1100 470]);
subplot(1, 2, 1); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-1.6, 1.6, 300), linspace(-0.4, 1.6, 300));
Z = 100 * (X2 - X1.^2).^2 + (1 - X1).^2;
contour(X1, X2, Z, [0.5 2 5 15 40 100 250], 'LineColor', [0.6 0.6 0.6]);
B = o2.base;
plot(B(1, :), B(2, :), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 9);
plot(x0(1), x0(2), 'ko', 'MarkerFaceColor', 'k');
plot(1, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(-1.5, 1.25, 'x_0 = (-1.2, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
text(0.55, 1.35, 'x^* = (1, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
axis([-1.6 1.6 -0.4 1.6]);
xlabel('x_1'); ylabel('x_2');
title(sprintf('Hooke-Jeeves: %d diem co so', size(B, 2)));
subplot(1, 2, 2);
semilogy(o2.nfevHist, o2.fHist, 'b.-', 'LineWidth', 1.2); hold on; grid on;
semilogy(o2.nfevHist, o2.deltaHist, 'k--', 'LineWidth', 1.2);
xlabel('so lan tinh R'); ylabel('gia tri');
legend('R(x_b) tai diem co so', '\Delta', 'Location', 'southwest');
title('Gia tri ham va buoc \Delta');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_rosenbrock.png'));
close(fig);

%% Hinh 3: vi du 4 - cac diem dung tren duong gay x1 = x2
fig = figure('visible', 'off');
hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-0.3, 3, 300), linspace(-0.3, 3, 300));
Z = 2 * abs(X1 - X2) + (X1 + X2 - 3).^2;
contour(X1, X2, Z, [0.04 0.25 0.5 1 1.5 2.5 4 6 9], 'LineColor', [0.6 0.6 0.6]);
plot([-0.3 3], [-0.3 3], 'k:');
plot(stuck(1, :), stuck(2, :), 'rv', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
plot(1.5, 1.5, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(1.45, 0.6, 'cac diem Hooke-Jeeves dung lai', 'Color', 'r', 'FontSize', 10, 'BackgroundColor', 'w');
text(0.35, 2.3, 'x^* = (1.5, 1.5)', 'FontSize', 10, 'BackgroundColor', 'w');
plot([1.5 0.85], [1.5 2.22], 'k-');
axis equal; axis([-0.3 3 -0.3 3]);
xlabel('x_1'); ylabel('x_2');
title('Vi du 4: diem dung phu thuoc \Delta_0 va x_0');
print(fig, '-dpng', '-r150', fullfile('figures', 'example4_nonsmooth.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
