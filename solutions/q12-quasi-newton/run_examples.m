% RUN_EXAMPLES  Cac vi du minh hoa thuat toan tua Newton BFGS va DFP (cau 12).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q06-bisection-search'));        % bisection_min.m
addpath(fullfile('..', 'q07-golden-section-search'));   % golden_section.m (cho gradient_descent)
addpath(fullfile('..', 'q08-gauss-seidel'));            % bracket_minimum.m
addpath(fullfile('..', 'q10-gradient-momentum'));       % gradient_descent.m
addpath(fullfile('..', 'q11-newton-method'));           % newton_method.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: ham bac hai n = 10 - ket thuc huu han va H_n = A^{-1}
%   f(x) = 1/2 x'Ax - b'x,  A = tridiag(-1, 2, -1)
fprintf('==================== VI DU 1 ====================\n');
n = 10;
A = 2 * eye(n) - diag(ones(n - 1, 1), 1) - diag(ones(n - 1, 1), -1);
Ainv = inv(A);
methods = {'bfgs', 'dfp'};
bs = {ones(n, 1), (1:n)'};
bname = {'b = (1,...,1)', 'b = (1,2,...,10)'};
o1 = cell(2, 2);
for j = 1:2
  b = bs{j};
  fq = @(x) 0.5 * x' * A * x - b' * x;  gq = @(x) A * x - b;
  for m = 1:2
    [x, ~, o] = quasi_newton(fq, gq, zeros(n, 1), struct('method', methods{m}, 'lineSearch', 'exact', 'tol', 1e-10));
    S = o.S;  C = S' * A * S;
    fprintf('%s, %s, tim kiem chinh xac: %d buoc, ||x - x*|| = %.1e, ||H - A^-1|| = %.1e, max |s_i''A s_j| (i~=j) = %.1e\n', ...
      bname{j}, upper(methods{m}), o.iterations, norm(x - A \ b), norm(o.H - Ainv), max(max(abs(C - diag(diag(C))))));
    if j == 1, fprintf('   ||grad f(x_k)||: %s\n', mat2str(o.G, 3)); end
    o1{j, m} = o;
  end
end
% Theo doi ||H_k - A^{-1}|| qua tung buoc (b = 1..10): chay lai voi maxIter = k
b = bs{2};  fq = @(x) 0.5 * x' * A * x - b' * x;  gq = @(x) A * x - b;
Herr = zeros(2, n + 1);
for m = 1:2
  for k = 0:n
    [~, ~, o] = quasi_newton(fq, gq, zeros(n, 1), struct('method', methods{m}, 'lineSearch', 'exact', ...
      'tol', 0, 'maxIter', k, 'scaleH0', false));
    Herr(m, k + 1) = norm(o.H - Ainv);
  end
  fprintf('%s (H0 = I): ||H_k - A^-1||, k = 0..10: %s\n', upper(methods{m}), mat2str(Herr(m, :), 3));
end
% Tim kiem Wolfe va so sanh voi cac cau truoc (b = 1,...,1)
b = bs{1};  fq = @(x) 0.5 * x' * A * x - b' * x;  gq = @(x) A * x - b;
for m = 1:2
  [x, ~, o] = quasi_newton(fq, gq, zeros(n, 1), struct('method', methods{m}, 'tol', 1e-8));
  fprintf('%s, b = (1,...,1), tim kiem Wolfe, tol 1e-8: %d buoc, %d lan tinh f, %d lan tinh grad, ||x - x*|| = %.1e\n', ...
    upper(methods{m}), o.iterations, o.nfev, o.ngev, norm(x - A \ b));
end
fprintf('\n');

%% Vi du 2: ham Rosenbrock - BFGS va DFP voi tim kiem Wolfe
%   R(x) = 100 (x2 - x1^2)^2 + (1 - x1)^2,  x0 = (-1.2, 1)
fprintf('==================== VI DU 2 ====================\n');
R  = @(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
RG = @(x) [-400 * x(1) * (x(2) - x(1)^2) - 2 * (1 - x(1)); 200 * (x(2) - x(1)^2)];
RH = @(x) [1200 * x(1)^2 - 400 * x(2) + 2, -400 * x(1); -400 * x(1), 200];
x02 = [-1.2; 1];
opt2 = struct('tol', 1e-10, 'maxIter', 5000);
fprintf('-- BFGS, Wolfe (c1 = 1e-4, c2 = 0.9), 12 buoc dau --\n');
quasi_newton(R, RG, x02, setfield(setfield(opt2, 'display', true), 'maxIter', 12));
[x2b, ~, o2b] = quasi_newton(R, RG, x02, opt2);
[x2d, ~, o2d] = quasi_newton(R, RG, x02, setfield(opt2, 'method', 'dfp'));
[x2n, ~, o2n] = newton_method(R, RG, RH, x02, struct('variant', 'damped'));
[x2g, ~, o2g] = gradient_descent(R, RG, x02, struct('step', 'armijo', 'tol', 1e-10, 'maxIter', 1e6));
e2b = sqrt(sum((o2b.X - 1).^2, 1));  e2d = sqrt(sum((o2d.X - 1).^2, 1));
fprintf('%-28s %7s %9s %9s %12s\n', 'Phuong phap (tol 1e-10)', 'buoc', 'tinh f', 'tinh grad', '||x - x*||');
fprintf('%-28s %7d %9d %9d %12.1e\n', 'BFGS, Wolfe', o2b.iterations, o2b.nfev, o2b.ngev, norm(x2b - 1));
fprintf('%-28s %7d %9d %9d %12.1e\n', 'DFP, Wolfe', o2d.iterations, o2d.nfev, o2d.ngev, norm(x2d - 1));
fprintf('%-28s %7d %9s %9s %12.1e  (+ %d Hessian)\n', 'Newton giam buoc (cau 11)', o2n.iterations, '-', '-', ...
  norm(x2n - 1), o2n.iterations);
fprintf('%-28s %7d %9d %9d %12.1e\n', 'Gradient Armijo (cau 10)', o2g.iterations, o2g.nfev, o2g.ngev, norm(x2g - 1));
fprintf('BFGS: ti so e_{k+1}/e_k o 8 buoc cuoi: %s\n', mat2str(e2b(end-7:end) ./ e2b(end-8:end-1), 3));
fprintf('BFGS: ty le buoc alpha = 1: %.2f; H cuoi co tri rieng %s, nghich dao Hessian tai x*: %s\n', ...
  mean(o2b.alphas == 1), mat2str(eig(o2b.H)', 4), mat2str(eig(inv(RH([1; 1])))', 4));
fprintf('DFP: R(x_k) tai k = 50, 100, 200, 400, 800: %s\n', mat2str(o2d.F([51 101 201 401 801]), 3));
% Do nhay voi c2 va voi viec co dinh co H0
fprintf('So buoc theo c2 va cach chon H0:\n%8s %14s %14s %14s %14s\n', 'c2', 'BFGS, H0 co', 'BFGS, H0 = I', 'DFP, H0 co', 'DFP, H0 = I');
c2s = [0.9 0.5 0.1 0.01];
tab2 = zeros(numel(c2s), 4);
for i = 1:numel(c2s)
  c = 0;
  for m = 1:2
    for sc = [true false]
      c = c + 1;
      [~, ~, o] = quasi_newton(R, RG, x02, struct('method', methods{m}, 'c2', c2s(i), 'scaleH0', sc, ...
        'tol', 1e-10, 'maxIter', 5000));
      tab2(i, c) = o.iterations;
    end
  end
  fprintf('%8g %14d %14d %14d %14d\n', c2s(i), tab2(i, :));
end
[~, ~, o2be] = quasi_newton(R, RG, x02, setfield(opt2, 'lineSearch', 'exact'));
[~, ~, o2de] = quasi_newton(R, RG, x02, setfield(setfield(opt2, 'lineSearch', 'exact'), 'method', 'dfp'));
fprintf('Tim kiem chinh xac: BFGS %d buoc, DFP %d buoc, max |x_k(BFGS) - x_k(DFP)| = %.1e\n\n', ...
  o2be.iterations, o2de.iterations, max(max(abs(o2be.X - o2de.X(:, 1:size(o2be.X, 2))))));

%% Vi du 3: ham Bazaraa - Hessian suy bien tai x*
%   F(x) = (x1 - 2)^4 + (x1 - 2x2)^2,  x0 = (0, 3)
fprintf('==================== VI DU 3 ====================\n');
F = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
G = @(x) [4*(x(1) - 2)^3 + 2*(x(1) - 2*x(2)); -4*(x(1) - 2*x(2))];
[~, ~, ~] = quasi_newton(F, G, [0; 3], struct('display', true, 'maxIter', 6));
for m = 1:2
  [x, Fx, o] = quasi_newton(F, G, [0; 3], struct('method', methods{m}, 'tol', 1e-8));
  e = sqrt(sum((o.X - [2; 1]).^2, 1));
  fprintf('%s, tol 1e-8: %d buoc, %d lan tinh f, %d lan tinh grad, x = %s, F = %.2e, ||x - x*|| = %.2e\n', ...
    upper(methods{m}), o.iterations, o.nfev, o.ngev, mat2str(x', 6), Fx, e(end));
  if m == 1
    fprintf('   ti so e_{k+1}/e_k o 6 buoc cuoi: %s\n', mat2str(e(end-5:end) ./ e(end-6:end-1), 3));
    o3b = o;
  end
end
fprintf('\n');

%% Vi du 4: ham khong bac hai n = 10 (vi du 1 cua cau 11)
fprintf('==================== VI DU 4 ====================\n');
b = ones(n, 1);
f4 = @(x) 0.5 * x' * A * x - b' * x + sum(x.^4) / 4;
g4 = @(x) A * x - b + x.^3;
h4 = @(x) A + diag(3 * x.^2);
x04 = 2 * ones(n, 1);
[x4n, ~, o4n] = newton_method(f4, g4, h4, x04, struct('tol', 1e-12));
x4s = x4n;
fprintf('%-24s %7s %9s %9s %12s\n', 'Phuong phap (tol 1e-8)', 'buoc', 'tinh f', 'tinh grad', '||x - x*||');
o4 = cell(1, 2);
for m = 1:2
  [x, ~, o] = quasi_newton(f4, g4, x04, struct('method', methods{m}, 'tol', 1e-8));
  fprintf('%-24s %7d %9d %9d %12.1e\n', [upper(methods{m}) ', Wolfe'], o.iterations, o.nfev, o.ngev, norm(x - x4s));
  o4{m} = o;
end
[x, ~, o] = newton_method(f4, g4, h4, x04, struct('tol', 1e-8));
fprintf('%-24s %7d %9d %9d %12.1e  (+ %d Hessian)\n', 'Newton (cau 11)', o.iterations, o.iterations + 1, ...
  o.iterations + 1, norm(x - x4s), o.iterations);
[x, ~, o] = gradient_descent(f4, g4, x04, struct('step', 'armijo', 'tol', 1e-8));
fprintf('%-24s %7d %9d %9d %12.1e\n', 'Gradient Armijo', o.iterations, o.nfev, o.ngev, norm(x - x4s));
e4 = sqrt(sum((o4{1}.X - x4s).^2, 1));
fprintf('BFGS: e_k = %s\n', mat2str(e4, 3));
fprintf('BFGS: e_{k+1}/e_k = %s\n', mat2str(e4(2:end) ./ e4(1:end-1), 3));
fprintf('BFGS: ||H_k - (Hessian tai x*)^-1|| cuoi = %.2e (so voi ||(Hessian)^-1|| = %.2f)\n\n', ...
  norm(o4{1}.H - inv(h4(x4s))), norm(inv(h4(x4s))));

%% Hinh 1: vi du 1 - ket thuc huu han va ||H_k - A^{-1}||
fig = figure('visible', 'off', 'Position', [100 100 1100 440]);
subplot(1, 2, 1);
semilogy(0:o1{2, 1}.iterations, o1{2, 1}.G, 'r.-', 'LineWidth', 1.3, 'MarkerSize', 13); hold on; grid on;
semilogy(0:o1{2, 2}.iterations, o1{2, 2}.G, 'bo--', 'LineWidth', 1.1, 'MarkerSize', 6);
semilogy(0:o1{1, 1}.iterations, o1{1, 1}.G, 'k.-', 'LineWidth', 1.1, 'MarkerSize', 13);
xlim([0 11]); ylim([1e-16 1e14]); set(gca, 'YTick', 10.^(-16:4:12));
xlabel('buoc lap k'); ylabel('||\nabla f(x_k)||');
legend('BFGS, b = (1,...,10)', 'DFP, b = (1,...,10)', 'BFGS, b = (1,...,1)', 'Location', 'northeast');
title('Vi du 1: tim kiem chinh xac, n = 10');
subplot(1, 2, 2);
semilogy(0:n, Herr(1, :), 'r.-', 'LineWidth', 1.3, 'MarkerSize', 13); hold on; grid on;
semilogy(0:n, Herr(2, :), 'bo--', 'LineWidth', 1.1, 'MarkerSize', 6);
xlim([0 10]);
xlabel('buoc lap k'); ylabel('||H_k - A^{-1}||');
legend('BFGS', 'DFP', 'Location', 'southwest');
title('Vi du 1: H_k \rightarrow A^{-1} sau n buoc (H_0 = I)');
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_quadratic.png'));
close(fig);

%% Hinh 2: vi du 2 - Rosenbrock
fig = figure('visible', 'off', 'Position', [100 100 1150 460]);
subplot(1, 2, 1); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-1.6, 1.6, 400), linspace(-0.6, 2.6, 400));
Z = 100 * (X2 - X1.^2).^2 + (1 - X1).^2;
contour(X1, X2, Z, [0.5 2 5 10 25 50 100 200 400], 'LineColor', [0.7 0.7 0.7]);
hb = plot(o2b.X(1, :), o2b.X(2, :), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 11);
hd = plot(o2d.X(1, :), o2d.X(2, :), 'b-', 'LineWidth', 0.8);
plot(x02(1), x02(2), 'ko', 'MarkerFaceColor', 'k');
plot(1, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(-1.5, 1.2, 'x_0', 'FontSize', 10);
text(1.08, 0.9, 'x^*', 'FontSize', 10);
axis([-1.6 1.6 -0.6 2.6]);
xlabel('x_1'); ylabel('x_2');
legend([hb, hd], {sprintf('BFGS (%d buoc)', o2b.iterations), sprintf('DFP (%d buoc)', o2d.iterations)}, ...
  'Location', 'northeast');
title('Vi du 2: ham Rosenbrock, tim kiem Wolfe');
subplot(1, 2, 2);
semilogy(0:o2b.iterations, max(e2b, 1e-17), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 11); hold on; grid on;
semilogy(0:o2d.iterations, max(e2d, 1e-17), 'b-', 'LineWidth', 1.1);
semilogy(0:o2n.iterations, max(sqrt(sum((o2n.X - 1).^2, 1)), 1e-17), 'k.--', 'LineWidth', 1.1, 'MarkerSize', 11);
semilogy(0:o2g.iterations, sqrt(sum((o2g.X - 1).^2, 1)), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.1);
xlim([0 150]); ylim([1e-17 1e14]); set(gca, 'YTick', 10.^(-16:4:12));
xlabel('buoc lap k'); ylabel('||x_k - x^*||');
legend('BFGS', sprintf('DFP (%d buoc)', o2d.iterations), 'Newton giam buoc (cau 11)', ...
  sprintf('Gradient Armijo (%d buoc)', o2g.iterations), 'Location', 'northeast');
title('Sai so (thang log)');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_rosenbrock.png'));
close(fig);

%% Hinh 3: vi du 3 va 4 - toc do hoi tu
fig = figure('visible', 'off', 'Position', [100 100 1100 440]);
subplot(1, 2, 2);
semilogy(0:o4n.iterations, max(sqrt(sum((o4n.X - x4s).^2, 1)), 1e-17), 'k.--', 'LineWidth', 1.1, 'MarkerSize', 12);
hold on; grid on;
semilogy(0:o4{1}.iterations, max(e4, 1e-17), 'r.-', 'LineWidth', 1.3, 'MarkerSize', 12);
semilogy(0:o4{2}.iterations, max(sqrt(sum((o4{2}.X - x4s).^2, 1)), 1e-17), 'bo-', 'LineWidth', 1.1, 'MarkerSize', 5);
ylim([1e-17 1e2]);
xlabel('buoc lap k'); ylabel('||x_k - x^*||');
legend('Newton', 'BFGS', 'DFP', 'Location', 'northeast');
title('Vi du 4: n = 10, khong bac hai');
subplot(1, 2, 1);
e3 = sqrt(sum((o3b.X - [2; 1]).^2, 1));
[~, ~, o3n] = newton_method(F, G, @(x) [12*(x(1) - 2)^2 + 2, -4; -4, 8], [0; 3], struct('tol', 1e-8));
semilogy(0:o3n.iterations, sqrt(sum((o3n.X - [2; 1]).^2, 1)), 'k.--', 'LineWidth', 1.1, 'MarkerSize', 12);
hold on; grid on;
semilogy(0:o3b.iterations, e3, 'r.-', 'LineWidth', 1.3, 'MarkerSize', 12);
xlabel('buoc lap k'); ylabel('||x_k - x^*||');
legend('Newton (cau 11)', 'BFGS', 'Location', 'northeast');
title('Vi du 3: ham Bazaraa');
print(fig, '-dpng', '-r150', fullfile('figures', 'examples3_4_convergence.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
