% RUN_EXAMPLES  Cac vi du minh hoa giai thuat Frank-Wolfe (cau 14).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q05-lagrange-duality'));   % solve_lp.m, solve_qp.m
addpath(fullfile('..', 'q06-bisection-search'));   % bisection_min.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: nghiem nam tren mot canh - hien tuong zic-zac
%   min (x1 - 2.5)^2 + 4 (x2 - 1.25)^2   s.t.  x1 + 2 x2 <= 4,  x >= 0
fprintf('==================== VI DU 1 ====================\n');
f1 = @(x) (x(1) - 2.5)^2 + 4 * (x(2) - 1.25)^2;
g1 = @(x) [2 * (x(1) - 2.5); 8 * (x(2) - 1.25)];
A1 = [1 2];  b1 = 4;  lb1 = [0; 0];  x01 = [1; 0];
xs1 = [2; 1];  fs1 = f1(xs1);
[xq, fq] = solve_qp(diag([2 8]), [-5; -10], A1, b1, [], [], lb1);
fprintf('Solver qp: x* = %s, f* = %.6g (giai tay: x* = (2, 1), f* = 0.5, lambda = 1)\n', mat2str(xq', 6), fq + 12.5);
fprintf('-- Frank-Wolfe, tim kiem chinh xac, 6 buoc dau --\n');
frank_wolfe(f1, g1, A1, b1, [], [], lb1, x01, struct('display', true, 'maxIter', 6));
names = {'FW, tim kiem chinh xac', 'FW, gamma = 2/(k+2)', 'FW co buoc lui (away)'};
optsList = {struct('step', 'exact'), struct('step', 'standard'), struct('variant', 'away')};
o1 = cell(1, 3);
fprintf('%-26s %10s %10s %12s %12s\n', 'Phuong phap', 'gap<=1e-3', 'gap<=1e-6', 'f - f*', '||x - x*||');
for j = 1:3
  op = optsList{j};  op.maxIter = 20000;
  op.tol = 1e-3;  [~, ~, oa] = frank_wolfe(f1, g1, A1, b1, [], [], lb1, x01, op);
  op.tol = 1e-6;  [x, fx, ob] = frank_wolfe(f1, g1, A1, b1, [], [], lb1, x01, op);
  if ob.exitflag == 1, sb = sprintf('%d', ob.iterations); else, sb = sprintf('>%d', op.maxIter); end
  fprintf('%-26s %10d %10s %12.1e %12.1e\n', names{j}, oa.iterations, sb, fx - fs1, norm(x - xs1));
  o1{j} = ob;
end
% Doi chieu voi danh gia ly thuyet f(x_k) - f* <= 2 L diam^2 / (k + 2)
Lc = 8;  diam = norm([4; 0] - [0; 2]);
for k = [10 100 1000 10000 20000]
  fprintf('k = %5d: f - f* = %.3e (chinh xac), %.3e (2/(k+2)); can tren 2*L*diam^2/(k+2) = %.3e\n', ...
    k, o1{1}.F(k + 1) - fs1, o1{2}.F(k + 1) - fs1, 2 * Lc * diam^2 / (k + 2));
end
kk = [1e2 1e3 1e4 2e4];
fprintf('(k+2)*(f - f*) voi buoc chinh xac tai k = 1e2, 1e3, 1e4, 2e4: %s\n', ...
  mat2str((kk + 2) .* (o1{1}.F(kk + 1) - fs1), 4));
fprintf('(k+2)^2*(f - f*) voi buoc 2/(k+2) tai k = 1e2, 1e3, 1e4, 2e4: %s\n', ...
  mat2str((kk + 2).^2 .* (o1{2}.F(kk + 1) - fs1), 4));
fprintf('-- FW co buoc lui, tung buoc --\n');
frank_wolfe(f1, g1, A1, b1, [], [], lb1, x01, struct('variant', 'away', 'display', true));
fprintf('\n');

%% Vi du 2: nghiem nam trong mien - hoi tu tuyen tinh
%   min (x1 - 1)^2 + 4 (x2 - 0.5)^2  tren cung tam giac
fprintf('==================== VI DU 2 ====================\n');
f2 = @(x) (x(1) - 1)^2 + 4 * (x(2) - 0.5)^2;
g2 = @(x) [2 * (x(1) - 1); 8 * (x(2) - 0.5)];
[x2, f2v, o2] = frank_wolfe(f2, g2, A1, b1, [], [], lb1, x01, struct('tol', 1e-10));
e2 = sqrt(sum((o2.X - [1; 0.5]).^2, 1));
fprintf('FW chinh xac: %d buoc, x = %s, ||x - x*|| = %.1e\n', o2.iterations, mat2str(x2', 8), e2(end));
fprintf('Ti so f(x_{k+1}) - f* / f(x_k) - f*, 10 buoc dau: %s\n', mat2str(o2.F(2:11) ./ o2.F(1:10), 3));
fprintf('\n');

%% Vi du 3: danh muc dau tu Markowitz, n = 50, rang buoc don hinh
%   min x'Sigma x - tau mu'x   s.t.  sum x = 1,  x >= 0
fprintf('==================== VI DU 3 ====================\n');
n = 50;  i = (1:n)';
mu  = 0.05 + 0.1 * sin(i).^2;                          % loi suat ky vong
Fm  = cos(i * (1:3)) / sqrt(n);                        % 3 nhan to rui ro
Sig = Fm * Fm' + diag(0.01 + 0.01 * (1 + cos(3 * i)));  % ma tran hiep phuong sai
tau = 1;
f3 = @(x) x' * Sig * x - tau * mu' * x;
g3 = @(x) 2 * Sig * x - tau * mu;
[xq3, fq3] = solve_qp(2 * Sig, -tau * mu, [], [], ones(1, n), 1, zeros(n, 1));
fprintf('Solver qp: f* = %.10g, so tai san co ty trong > 1e-8: %d\n', fq3, nnz(xq3 > 1e-8));
o3 = cell(1, 3);
for j = 1:3
  op = optsList{j};  op.tol = 1e-8;  op.maxIter = 20000;
  [x, fx, o] = frank_wolfe(f3, g3, [], [], ones(1, n), 1, zeros(n, 1), [], op);
  fprintf('%-26s %6d buoc (exitflag %d), f - f* = %.1e, ||x - x_qp|| = %.1e, %d tai san, %d buoc lui, %d lan bo dinh\n', ...
    names{j}, o.iterations, o.exitflag, fx - fq3, norm(x - xq3), nnz(x > 1e-10), o.awaySteps, o.dropSteps);
  o3{j} = o;
end
fprintf('FW chinh xac: so tai san trong x_k tai k = 5, 10, 20: %s (x_k la to hop loi cua <= k+1 dinh)\n', ...
  mat2str([nnz(o3{1}.X(:, 6) > 0), nnz(o3{1}.X(:, 11) > 0), nnz(o3{1}.X(:, 21) > 0)]));
fprintf('Khe FW la chung chi: gap_k >= f(x_k) - f* tai moi buoc? %d\n', ...
  all(o3{1}.gap(1:end) >= o3{1}.F - fq3 - 1e-12));
fprintf('\n');

%% Hinh 1: vi du 1 - duong di
fig = figure('visible', 'off', 'Position', [100 100 1150 460]);
[X1, X2] = meshgrid(linspace(-0.2, 4.2, 300), linspace(-0.2, 2.2, 300));
Z = (X1 - 2.5).^2 + 4 * (X2 - 1.25).^2;
ttl = {sprintf('FW, tim kiem chinh xac (30 buoc dau)'), sprintf('FW co buoc lui: %d buoc', o1{3}.iterations)};
for p = 1:2
  subplot(1, 2, p); hold on; grid on; box on;
  fill([0 4 0], [0 0 2], [0.85 0.92 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 1.2);
  contour(X1, X2, Z, [0.5 0.6 0.8 1.2 2 3.5 5.5 8], 'LineColor', [0.6 0.6 0.6]);
  if p == 1
    P = o1{1}.X(:, 1:31);  col = 'r';
  else
    P = o1{3}.X;  col = [0 0.5 0];
  end
  plot(P(1, :), P(2, :), '.-', 'Color', col, 'LineWidth', 1.1, 'MarkerSize', 11);
  plot(x01(1), x01(2), 'ko', 'MarkerFaceColor', 'k');
  plot(2, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
  text(1.05, 0.12, 'x_0', 'FontSize', 10);
  text(2.15, 1.12, 'x^* = (2, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
  text(4.0, -0.12, '(4, 0)', 'FontSize', 9);  text(0.05, 2.08, '(0, 2)', 'FontSize', 9);
  axis equal; axis([-0.2 4.4 -0.2 2.3]);
  xlabel('x_1'); ylabel('x_2');
  title(ttl{p});
end
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_paths.png'));
close(fig);

%% Hinh 2: toc do hoi tu
fig = figure('visible', 'off', 'Position', [100 100 1150 440]);
subplot(1, 2, 1);
k1 = 1:o1{1}.iterations;  k2 = 1:o1{2}.iterations;
loglog(k1, o1{1}.F(k1 + 1) - fs1, 'r-', 'LineWidth', 1.3); hold on; grid on;
loglog(k2, o1{2}.F(k2 + 1) - fs1, 'b-', 'LineWidth', 1.3);
loglog(k1, 2 * Lc * diam^2 ./ (k1 + 2), 'k--', 'LineWidth', 1);
xlabel('buoc lap k'); ylabel('f(x_k) - f^*');
ylim([1e-9 1e7]); set(gca, 'YTick', 10.^(-8:2:6));
legend('FW, tim kiem chinh xac', 'FW, \gamma_k = 2/(k+2)', '2 L diam^2/(k+2)', 'Location', 'northeast');
title('Vi du 1: nghiem tren canh');
subplot(1, 2, 2);
for j = 1:3
  K = 0:o3{j}.iterations;
  semilogy(K, max(o3{j}.F - fq3, 1e-16), '-', 'LineWidth', 1.3); hold on;
end
grid on; xlim([0 2000]); ylim([1e-14 1]);
xlabel('buoc lap k'); ylabel('f(x_k) - f^*');
legend(names, 'Location', 'northeast');
title('Vi du 3: Markowitz, n = 50');
print(fig, '-dpng', '-r150', fullfile('figures', 'examples1_3_convergence.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
