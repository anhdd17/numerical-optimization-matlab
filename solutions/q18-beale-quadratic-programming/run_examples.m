% RUN_EXAMPLES  Cac vi du minh hoa giai thuat don hinh Beale cho quy hoach toan phuong (cau 18).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q01-simplex'));              % simplex.m (cho pha 1)
addpath(fullfile('..', 'q02-two-phase-simplex'));    % two_phase_simplex.m
addpath(fullfile('..', 'q05-lagrange-duality'));     % solve_qp.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: quy hoach toan phuong cua cau 15 (dang chinh tac, x3, x4 la bien bu)
%   min -4x1 - 6x2 + 2x1^2 + 2x2^2 - 2x1x2
%   s.t. x1 + x2 + x3 = 2,  x1 + 5x2 + x4 = 5,  x >= 0
fprintf('==================== VI DU 1 ====================\n');
c1 = [-4; -6; 0; 0];
D1 = zeros(4);  D1(1:2, 1:2) = [4 -2; -2 4];
A1 = [1 1 1 0; 1 5 0 1];  b1 = [2; 5];
[x1, f1, fl1, o1] = beale_qp(c1, D1, A1, b1, struct('B0', [3 4], 'display', true));
fprintf('x* = %s, f* = %.10g (= -222/31 = %.10g), exitflag = %d\n', mat2str(x1', 8), f1, -222/31, fl1);
fprintf('df/dx_j tai nghiem (nhan tu cua x_j >= 0): %s; 32/31 = %.6f\n', mat2str(o1.lambda', 5), 32/31);
[~, ~, ~, o1b] = beale_qp(c1, D1, A1, b1);
fprintf('Co so xuat phat tu pha 1 (cau 2): %d buoc, cung nghiem\n', o1b.iterations);
[xq1, fq1] = solve_qp([4 -2; -2 4], [-4; -6], [1 1; 1 5], [2; 5], [], [], [0; 0]);   % dang bat dang thuc 2 bien
fprintf('Solver qp (dang 2 bien): x* = %s, f* = %.10g\n\n', mat2str(xq1', 8), fq1);

%% Vi du 2: bien tu do bi loai (bai toan cua cau 14, vi du 1)
%   min (x1 - 2.5)^2 + 4 (x2 - 1.25)^2 = -5x1 - 10x2 + x1^2 + 4x2^2 + 12.5
%   s.t. x1 + 2x2 + x3 = 4,  x >= 0
fprintf('==================== VI DU 2 ====================\n');
c2 = [-5; -10; 0];  D2 = diag([2 8 0]);  A2 = [1 2 1];  b2 = 4;
[x2, f2, fl2, o2] = beale_qp(c2, D2, A2, b2, struct('B0', 3, 'display', true));
fprintf('x* = %s, f* = %.10g (+ 12.5 = %.10g), lambda = %s\n\n', mat2str(x2', 8), f2, f2 + 12.5, mat2str(o2.lambda', 5));

%% Vi du 3: ham muc tieu khong bi chan duoi (D suy bien)
%   min -x1 + 1/2 (x1 - x2)^2   s.t.  x1 - x2 + x3 = 1,  x >= 0
fprintf('==================== VI DU 3 ====================\n');
c3 = [-1; 0; 0];  D3 = [1 -1 0; -1 1 0; 0 0 0];  A3 = [1 -1 1];  b3 = 1;
[x3, f3, fl3, o3] = beale_qp(c3, D3, A3, b3, struct('B0', 3, 'display', true));
fprintf('exitflag = %d; f(t, t, 0) = -t -> -inf doc theo huong (1, 1, 0)\n\n', fl3);

%% Vi du 4: danh muc Markowitz n = 50 (cau 14, 15, 17)
fprintf('==================== VI DU 4 ====================\n');
n = 50;  i = (1:n)';
mu  = 0.05 + 0.1 * sin(i).^2;
Fm  = cos(i * (1:3)) / sqrt(n);
Sig = Fm * Fm' + diag(0.01 + 0.01 * (1 + cos(3 * i)));
[xq4, fq4] = solve_qp(2 * Sig, -mu, [], [], ones(1, n), 1, zeros(n, 1));
[x4, f4, fl4, o4] = beale_qp(-mu, 2 * Sig, ones(1, n), 1);
nPiv = sum(strncmp(o4.events, 'xoay', 4));
nDrop = sum(~cellfun(@isempty, strfind(o4.events, 'loai')));
fprintf('Beale: exitflag = %d, %d buoc (%d lan xoay, %d bien tu do tao ra, %d bi loai)\n', ...
  fl4, o4.iterations, nPiv, o4.nFree, nDrop);
fprintf('f = %.15g, f_qp = %.15g, f - f_qp = %.1e, ||x - x_qp|| = %.1e, %d tai san\n', f4, fq4, f4 - fq4, ...
  norm(x4 - xq4), nnz(x4 > 1e-10));
nb = setdiff(1:n, o4.basis);
fprintf('Diem xuat phat (pha 1): tai san %d; min df/dx_j tren %d bien phi co so = %.3e (>= 0)\n', ...
  find(o4.X(:, 1) > 0), numel(nb), min(o4.lambda(nb)));
fprintf('So tai san co ty trong duong theo buoc: %s\n\n', mat2str(sum(o4.X > 1e-12, 1)));

%% Hinh: duong di cua vi du 1 va 2
fig = figure('visible', 'off', 'Position', [100 100 1150 470]);
subplot(1, 2, 1); hold on; grid on; box on;
fill([0 2 1.25 0], [0 0 0.75 1], [0.85 0.92 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 1.2);
[X1, X2] = meshgrid(linspace(-0.3, 2.3, 300), linspace(-0.3, 1.6, 300));
contour(X1, X2, 2 * X1.^2 + 2 * X2.^2 - 2 * X1 .* X2 - 4 * X1 - 6 * X2, [-7.1 -6.8 -6 -5 -4 -3 -2 -1], ...
  'LineColor', [0.65 0.65 0.65]);
plot(o1.X(1, :), o1.X(2, :), 'r.-', 'LineWidth', 1.6, 'MarkerSize', 18);
text(0.05, -0.12, 'x_0: co so \{x_3, x_4\}', 'FontSize', 9);
text(0.05, 1.08, 'xoay: x_2 vao, x_4 ra', 'FontSize', 9, 'BackgroundColor', 'w');
text(1.15, 0.92, 'tao u_1: \partial f/\partial x_1 = 0', 'FontSize', 9, 'BackgroundColor', 'w');
plot(35/31, 24/31, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
axis equal; axis([-0.3 2.3 -0.3 1.5]);
xlabel('x_1'); ylabel('x_2');
title(sprintf('Vi du 1: %d buoc', o1.iterations));
subplot(1, 2, 2); hold on; grid on; box on;
fill([0 4 0], [0 0 2], [0.85 0.92 1], 'EdgeColor', [0.2 0.4 0.8], 'LineWidth', 1.2);
[X1, X2] = meshgrid(linspace(-0.2, 4.2, 300), linspace(-0.2, 2.2, 300));
contour(X1, X2, (X1 - 2.5).^2 + 4 * (X2 - 1.25).^2, [0.5 1 2 3.5 5.5 8], 'LineColor', [0.65 0.65 0.65]);
plot(o2.X(1, :), o2.X(2, :), 'r.-', 'LineWidth', 1.6, 'MarkerSize', 18);
text(0.08, -0.12, 'x_0', 'FontSize', 9);
text(-0.15, 1.42, 'tao u_1', 'FontSize', 9, 'BackgroundColor', 'w');
text(1.15, 1.48, 'xoay: x_1 vao, x_3 ra', 'FontSize', 9, 'BackgroundColor', 'w');
text(2.15, 0.85, 'tao u_2, loai u_1', 'FontSize', 9, 'BackgroundColor', 'w');
plot(2, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
axis equal; axis([-0.3 4.3 -0.3 2.3]);
xlabel('x_1'); ylabel('x_2');
title(sprintf('Vi du 2: %d buoc', o2.iterations));
print(fig, '-dpng', '-r150', fullfile('figures', 'examples1_2_paths.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
